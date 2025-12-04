#!/usr/bin/env python3
"""
Convert Fomosto/Pyrocko Green's function store to SeisComP sc3gf1d format

Usage:
    python3 convert_to_sc3gf1d.py <fomosto_store_path> <output_dir>

Example:
    python3 convert_to_sc3gf1d.py regional_16hz gf_sc3gf1d/regional_16hz
"""

import sys
import numpy as np
from pathlib import Path
from datetime import datetime

try:
    from pyrocko import gf, trace
    from obspy import Trace
    from obspy.core import UTCDateTime
except ImportError as e:
    print(f"Error: {e}")
    print("Please install required packages:")
    print("  pip install pyrocko obspy")
    sys.exit(1)


def convert_store_to_sc3gf1d(store_path, output_dir):
    """
    Convert Fomosto store to sc3gf1d format

    Parameters:
    -----------
    store_path : str
        Path to Fomosto store directory
    output_dir : str
        Output directory for sc3gf1d files
    """
    print(f"Converting Fomosto store to sc3gf1d format")
    print(f"Input store: {store_path}")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

    # Load store
    try:
        store = gf.store.Store(str(store_path), 'r')
    except Exception as e:
        print(f"Error loading store: {e}")
        print("Make sure the path points to a valid Fomosto store")
        sys.exit(1)

    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    config = store.config

    # Get depth and distance ranges
    depths_m = store.config.coords[0]  # Source depths in meters
    distances_m = store.config.coords[1]  # Distances in meters

    depths_km = depths_m / 1000.0
    distances_km = distances_m / 1000.0

    print(f"Store configuration:")
    print(f"  Sample rate: {config.sample_rate} Hz")
    print(f"  Depths: {depths_km.min():.1f} - {depths_km.max():.1f} km ({len(depths_km)} values)")
    print(f"  Distances: {distances_km.min():.1f} - {distances_km.max():.1f} km ({len(distances_km)} values)")
    print(f"  Number of components: {config.ncomponents}")
    print("=" * 60)

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

    for depth_km in depths_km:
        for distance_km in distances_km:
            count += 1
            print(f"\r[{count}/{total}] Processing: depth={depth_km:.1f}km, distance={distance_km:.1f}km",
                  end="", flush=True)

            try:
                convert_gf_pair(
                    store, output_path,
                    depth_km, distance_km,
                    component_names
                )
                success += 1
            except Exception as e:
                print(f"\nWarning: Failed depth={depth_km:.1f}km, dist={distance_km:.1f}km: {e}")

    print(f"\n" + "=" * 60)
    print(f"Conversion complete: {success}/{total} successful")
    print(f"Output directory: {output_path.absolute()}")

    store.close()


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


def convert_gf_pair(store, output_path, depth_km, distance_km, component_names):
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

        for idx in range(ncomponents):
            comp_name = component_names[idx]

            # Extract GF trace for this component using store.get()
            # Args: (source_depth, distance, component)
            gf_trace = store.get((depth_m, distance_m, idx), interpolation='nearest_neighbor')

            if gf_trace is None or gf_trace.is_zero:
                continue

            # Get trace data
            data = gf_trace.data
            nsamples = len(data)

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

    except Exception as e:
        raise RuntimeError(f"Failed to extract GF: {e}")


def main():
    """Main entry point"""
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    store_path = sys.argv[1]
    output_dir = sys.argv[2]

    if not Path(store_path).exists():
        print(f"Error: Store path does not exist: {store_path}")
        sys.exit(1)

    convert_store_to_sc3gf1d(store_path, output_dir)


if __name__ == "__main__":
    main()
