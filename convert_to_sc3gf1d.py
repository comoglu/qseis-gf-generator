#!/usr/bin/env python3
"""
Convert Fomosto/Pyrocko Green's function store to SeisComP sc3gf1d format

Usage:
    python3 convert_to_sc3gf1d.py <fomosto_store_path> <output_dir> [options]

Options:
    --validate-only    Only validate the store, don't convert
    --skip-zeros       Skip empty/zero traces
    --verbose          Show detailed output
    --force            Overwrite existing output directory

Example:
    python3 convert_to_sc3gf1d.py regional_16hz gf_sc3gf1d/regional_16hz
    python3 convert_to_sc3gf1d.py regional_16hz gf_sc3gf1d/regional_16hz --verbose
"""

import sys
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime
import time

try:
    from pyrocko import gf, trace
    from obspy import Trace
    from obspy.core import UTCDateTime
except ImportError as e:
    print(f"Error: {e}")
    print("Please install required packages:")
    print("  pip install pyrocko obspy")
    sys.exit(1)


def validate_store(store):
    """
    Validate that store is ready for conversion

    Returns:
    --------
    tuple: (is_valid, warnings, errors)
    """
    warnings = []
    errors = []

    config = store.config

    # Check if store has been built
    try:
        # Try to access a sample GF to check if built
        depths_m = config.coords[0]
        distances_m = config.coords[1]

        if len(depths_m) == 0 or len(distances_m) == 0:
            errors.append("Store appears empty (no depth/distance coordinates)")
            return False, warnings, errors

        # Try to get one GF to verify store is built
        test_depth = depths_m[len(depths_m)//2]
        test_dist = distances_m[len(distances_m)//2]
        test_gf = store.get((test_depth, test_dist, 0), interpolation='nearest_neighbor')

        if test_gf is None:
            errors.append("Store has not been built yet (GFs are missing)")
            errors.append("Run: cd <store_path> && fomosto build --nworkers=8")
            return False, warnings, errors

    except Exception as e:
        errors.append(f"Error accessing store data: {e}")
        return False, warnings, errors

    # Check for reasonable time window
    if hasattr(config, 'deltat'):
        dt = config.deltat
    else:
        dt = 1.0 / config.sample_rate

    # Estimate time window from a sample trace
    if test_gf is not None and not test_gf.is_zero:
        duration = len(test_gf.data) * dt
        if duration < 100.0:
            warnings.append(f"⚠ Short time window: {duration}s (may truncate waveforms)")
        elif duration > 7200.0:
            warnings.append(f"⚠ Very long time window: {duration}s (large storage requirements)")

    # Check component count
    if config.ncomponents < 8:
        warnings.append(f"⚠ Only {config.ncomponents} components (expected 8-10 for full MT)")

    return True, warnings, errors


def convert_store_to_sc3gf1d(store_path, output_dir, skip_zeros=False, verbose=False, force=False):
    """
    Convert Fomosto store to sc3gf1d format

    Parameters:
    -----------
    store_path : str
        Path to Fomosto store directory
    output_dir : str
        Output directory for sc3gf1d files
    skip_zeros : bool
        Skip empty/zero traces
    verbose : bool
        Show detailed output
    force : bool
        Overwrite existing output directory
    """
    print("=" * 70)
    print("Converting Fomosto store to sc3gf1d format")
    print("=" * 70)
    print(f"Input store: {store_path}")
    print(f"Output directory: {output_dir}")
    print()

    # Load store
    try:
        store = gf.store.Store(str(store_path), 'r')
    except Exception as e:
        print(f"❌ Error loading store: {e}")
        print("Make sure the path points to a valid Fomosto store")
        sys.exit(1)

    # Validate store
    print("Validating store...")
    is_valid, warnings, errors = validate_store(store)

    if errors:
        print("\n❌ VALIDATION ERRORS:")
        for error in errors:
            print(f"  {error}")
        store.close()
        sys.exit(1)

    if warnings:
        print("\n⚠️  VALIDATION WARNINGS:")
        for warning in warnings:
            print(f"  {warning}")
        print()

    # Check output directory
    output_path = Path(output_dir)
    if output_path.exists() and not force:
        print(f"❌ Output directory already exists: {output_path}")
        print("Use --force to overwrite")
        store.close()
        sys.exit(1)

    output_path.mkdir(exist_ok=True, parents=True)

    config = store.config

    # Get depth and distance ranges
    depths_m = store.config.coords[0]  # Source depths in meters
    distances_m = store.config.coords[1]  # Distances in meters

    depths_km = depths_m / 1000.0
    distances_km = distances_m / 1000.0

    # Calculate expected duration
    if hasattr(config, 'deltat'):
        dt = config.deltat
    else:
        dt = 1.0 / config.sample_rate

    # Get sample trace to determine duration
    test_gf = store.get((depths_m[0], distances_m[0], 0), interpolation='nearest_neighbor')
    duration = len(test_gf.data) * dt if test_gf else 0.0

    print("Store configuration:")
    print(f"  Store ID: {config.id}")
    print(f"  Sample rate: {config.sample_rate} Hz")
    print(f"  Time window: {duration:.1f}s ({len(test_gf.data) if test_gf else 0} samples)")
    print(f"  Depths: {depths_km.min():.1f} - {depths_km.max():.1f} km ({len(depths_km)} values)")
    print(f"  Distances: {distances_km.min():.1f} - {distances_km.max():.1f} km ({len(distances_km)} values)")
    print(f"  Components: {config.ncomponents}")
    print("=" * 70)

    # Create description file
    create_description_file(output_path, depths_km, distances_km, config)

    # Component mapping for moment tensor sources
    # These are the standard sc3gf1d component names
    component_names = ['ZSS', 'ZDD', 'ZDS', 'RSS', 'RDD', 'RDS', 'TSS', 'TDS']

    if config.ncomponents >= 10:
        component_names.extend(['ZEP', 'REP'])

    # Convert each depth/distance pair
    total = len(depths_km) * len(distances_km)
    count = 0
    success = 0
    skipped = 0
    failed = 0
    failed_pairs = []

    start_time = time.time()

    print("\nConverting Green's Functions...")
    print(f"Total: {total} depth/distance pairs × {len(component_names)} components")
    print()

    for depth_km in depths_km:
        for distance_km in distances_km:
            count += 1

            # Progress with ETA
            if count > 1:
                elapsed = time.time() - start_time
                rate = count / elapsed
                remaining = (total - count) / rate if rate > 0 else 0
                eta_str = f"ETA: {int(remaining//60)}m {int(remaining%60)}s"
            else:
                eta_str = "calculating..."

            if not verbose:
                print(f"\r[{count}/{total}] depth={depth_km:.1f}km, dist={distance_km:.1f}km | {eta_str}",
                      end="", flush=True)
            else:
                print(f"[{count}/{total}] Processing depth={depth_km:.1f}km, distance={distance_km:.1f}km")

            try:
                n_written = convert_gf_pair(
                    store, output_path,
                    depth_km, distance_km,
                    component_names,
                    skip_zeros=skip_zeros,
                    verbose=verbose
                )
                if n_written > 0:
                    success += 1
                else:
                    skipped += 1
                    if verbose:
                        print(f"  Skipped (all zeros)")
            except Exception as e:
                failed += 1
                failed_pairs.append((depth_km, distance_km, str(e)))
                if verbose:
                    print(f"  ❌ Failed: {e}")

    elapsed_total = time.time() - start_time

    print(f"\n" + "=" * 70)
    print("Conversion Summary:")
    print(f"  Total pairs: {total}")
    print(f"  ✓ Successful: {success}")
    if skipped > 0:
        print(f"  ○ Skipped (zeros): {skipped}")
    if failed > 0:
        print(f"  ✗ Failed: {failed}")
    print(f"  Time elapsed: {int(elapsed_total//60)}m {int(elapsed_total%60)}s")
    print(f"  Output: {output_path.absolute()}")

    if failed_pairs and verbose:
        print("\nFailed pairs:")
        for depth, dist, err in failed_pairs[:10]:  # Show first 10
            print(f"  depth={depth:.1f}km, dist={dist:.1f}km: {err}")
        if len(failed_pairs) > 10:
            print(f"  ... and {len(failed_pairs)-10} more")

    print("=" * 70)

    store.close()

    return success, skipped, failed


def create_description_file(output_path, depths, distances, config):
    """Create sc3gf1d description file"""
    store_name = output_path.name
    desc_file = output_path.parent / f"{store_name}.desc"

    with open(desc_file, 'w') as f:
        f.write(f"# sc3gf1d Green's Functions: {store_name}\n")
        f.write(f"# Created: {datetime.now().isoformat()}\n")
        f.write(f"# Source: Fomosto/Pyrocko Store\n")
        f.write(f"# Sample rate: {config.sample_rate} Hz\n")
        f.write("#\n")

        # Depth ranges
        if len(depths) > 1:
            depth_step = depths[1] - depths[0]
            f.write(f"depth {depths.min():.1f} {depths.max():.1f} {depth_step:.1f}\n")
        else:
            f.write(f"depth {depths[0]:.1f} {depths[0]:.1f} 1.0\n")

        # Distance ranges
        if len(distances) > 1:
            dist_step = distances[1] - distances[0]
            f.write(f"distance {distances.min():.1f} {distances.max():.1f} {dist_step:.1f}\n")
        else:
            f.write(f"distance {distances[0]:.1f} {distances[0]:.1f} 1.0\n")

        # Travel time interface
        f.write("#\n")
        f.write("# Travel time tables\n")
        f.write("times LOCSAT iasp91\n")

    print(f"Created description file: {desc_file}")


def convert_gf_pair(store, output_path, depth_km, distance_km, component_names,
                   skip_zeros=False, verbose=False):
    """
    Convert single depth/distance pair to SAC files

    Parameters:
    -----------
    store : pyrocko.gf.Store
        Fomosto store
    output_path : Path
        Output directory
    depth_km : float
        Source depth in km
    distance_km : float
        Epicentral distance in km
    component_names : list
        List of component names to extract
    skip_zeros : bool
        Skip empty/zero traces
    verbose : bool
        Show detailed output

    Returns:
    --------
    int : Number of traces written
    """
    # Create output directory structure
    # Format: <name>/<depth_in_100m>/<dist_in_km>/
    depth_dir = output_path / f"{int(depth_km * 10):04d}"
    dist_dir = depth_dir / f"{int(distance_km):05d}"
    dist_dir.mkdir(exist_ok=True, parents=True)

    depth_m = depth_km * 1000.0
    distance_m = distance_km * 1000.0

    # Reference time (1970-01-01 00:00:00)
    reftime = UTCDateTime(1970, 1, 1, 0, 0, 0)

    # Extract GF traces for this source-receiver pair
    try:
        # Get sampling parameters
        dt = 1.0 / store.config.sample_rate

        # Write each component as SAC file
        ncomponents = min(len(component_names), store.config.ncomponents)
        traces_written = 0

        for idx in range(ncomponents):
            comp_name = component_names[idx]

            # Extract GF trace for this component using store.get()
            # Args: (source_depth, distance, component)
            gf_trace = store.get((depth_m, distance_m, idx), interpolation='nearest_neighbor')

            if gf_trace is None:
                if verbose:
                    print(f"    Component {comp_name}: None (missing)")
                continue

            if gf_trace.is_zero:
                if skip_zeros:
                    if verbose:
                        print(f"    Component {comp_name}: skipped (zero)")
                    continue
                elif verbose:
                    print(f"    Component {comp_name}: zero trace")

            # Get trace data
            data = gf_trace.data
            nsamples = len(data)

            if verbose:
                max_amp = np.abs(data).max()
                print(f"    Component {comp_name}: {nsamples} samples, max amp={max_amp:.2e}")

            # Calculate time offset
            tmin = gf_trace.itmin * dt

            # Convert to ObsPy trace
            obs_tr = Trace(data=data)
            obs_tr.stats.sampling_rate = store.config.sample_rate
            obs_tr.stats.starttime = reftime + tmin
            obs_tr.stats.delta = dt
            obs_tr.stats.npts = nsamples

            # Station/channel info
            obs_tr.stats.station = f"D{int(distance_km):05d}"
            obs_tr.stats.channel = comp_name[:3]
            obs_tr.stats.network = "SY"
            obs_tr.stats.location = ""

            # SAC headers
            tmax = tmin + (nsamples - 1) * dt
            obs_tr.stats.sac = {
                'stla': 0.0,
                'stlo': distance_km / 111.195,  # Approximate deg
                'stel': 0.0,
                'evla': 0.0,
                'evlo': 0.0,
                'evdp': depth_km,
                'dist': distance_km,
                'az': 0.0,
                'baz': 180.0,
                'b': float(tmin),
                'e': float(tmax),
                'iztype': 9,  # IB - begin time
            }

            # Write SAC file
            # Format: <depth_in_100m>.<dist_in_km>.<component>
            sac_filename = f"{int(depth_km * 10):04d}.{int(distance_km):05d}.{comp_name}"
            sac_path = dist_dir / sac_filename
            obs_tr.write(str(sac_path), format='SAC')
            traces_written += 1

        return traces_written

    except Exception as e:
        raise RuntimeError(f"Failed to extract GF: {e}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Convert Fomosto store to SeisComP sc3gf1d format',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('store_path', help='Path to Fomosto store directory')
    parser.add_argument('output_dir', help='Output directory for sc3gf1d files')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate store, do not convert')
    parser.add_argument('--skip-zeros', action='store_true',
                       help='Skip empty/zero traces')
    parser.add_argument('--verbose', action='store_true',
                       help='Show detailed output for each trace')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing output directory')

    args = parser.parse_args()

    if not Path(args.store_path).exists():
        print(f"❌ Error: Store path does not exist: {args.store_path}")
        sys.exit(1)

    # Validate-only mode
    if args.validate_only:
        print("=" * 70)
        print("Validating Fomosto store (validation only)")
        print("=" * 70)
        print(f"Store path: {args.store_path}\n")

        try:
            store = gf.store.Store(str(args.store_path), 'r')
        except Exception as e:
            print(f"❌ Error loading store: {e}")
            sys.exit(1)

        is_valid, warnings, errors = validate_store(store)

        config = store.config
        depths_m = config.coords[0]
        distances_m = config.coords[1]
        depths_km = depths_m / 1000.0
        distances_km = distances_m / 1000.0

        print("Store Information:")
        print(f"  Store ID: {config.id}")
        print(f"  Sample rate: {config.sample_rate} Hz")
        print(f"  Depths: {depths_km.min():.1f} - {depths_km.max():.1f} km ({len(depths_km)} values)")
        print(f"  Distances: {distances_km.min():.1f} - {distances_km.max():.1f} km ({len(distances_km)} values)")
        print(f"  Components: {config.ncomponents}")

        if warnings:
            print("\n⚠️  WARNINGS:")
            for warning in warnings:
                print(f"  {warning}")

        if errors:
            print("\n❌ ERRORS:")
            for error in errors:
                print(f"  {error}")
            store.close()
            sys.exit(1)

        print("\n✓ Store validation passed")
        store.close()
        sys.exit(0)

    # Normal conversion mode
    convert_store_to_sc3gf1d(
        args.store_path,
        args.output_dir,
        skip_zeros=args.skip_zeros,
        verbose=args.verbose,
        force=args.force
    )


if __name__ == "__main__":
    main()
