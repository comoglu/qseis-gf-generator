#!/usr/bin/env python3
"""
Generate Green's Functions for Any Region using Fomosto/QSEIS

Flexible GF generator that supports:
- Custom velocity models (load from file)
- Built-in models (ak135, iasp91, prem)
- Any distance/depth range
- Configurable sampling rate

Usage:
    python3 generate_gf_custom.py --config myconfig.yaml
    python3 generate_gf_custom.py --model iasp91 --sample-rate 10 --distance-max 1000

Author: Generated for regional seismology
"""

import sys
import argparse
import yaml
from pathlib import Path

try:
    from pyrocko import gf
    from pyrocko.fomosto import qseis
    from pyrocko import cake
except ImportError:
    print("Error: Pyrocko is not installed")
    print("Install with: pip install pyrocko")
    sys.exit(1)


def load_velocity_model(model_spec):
    """
    Load velocity model from various sources

    Parameters:
    -----------
    model_spec : str
        Can be:
        - Built-in model name: 'iasp91', 'ak135-f-average', 'prem-no-ocean'
        - Path to .tvel file (travel time format)
        - Path to ND model file (.nd format)
        - Path to simple ASCII file with layers

    Returns:
    --------
    cake.LayeredModel
    """
    model_path = Path(model_spec)

    # Try built-in models first (use cake.builtin_models() to get proper paths)
    builtin_models = {
        'iasp91': 'iasp91-compat.m',
        'ak135-f-continental': 'ak135-f-continental.m',
        'ak135-f-average': 'ak135-f-average.m',
        'ak135-f-oceanic': 'ak135-f-oceanic.m',
        'prem-no-ocean': 'prem-no-ocean.m',
        'prem-no-crust': 'prem-no-crust.m'
    }

    if model_spec in builtin_models:
        print(f"Loading built-in Pyrocko model: {model_spec}")
        try:
            # Get the actual model file from Pyrocko's data directory
            from pyrocko import config as pconfig
            model_file = builtin_models[model_spec]
            # Try to load using Pyrocko's builtin function
            try:
                return cake.load_model(model_file)
            except:
                # Fallback: construct path manually
                import pyrocko
                pyrocko_dir = Path(pyrocko.__file__).parent
                model_path = pyrocko_dir / 'data' / model_file
                if model_path.exists():
                    return cake.load_model(str(model_path))
                raise FileNotFoundError(f"Built-in model not found: {model_file}")
        except Exception as e:
            print(f"Error loading built-in model: {e}")
            print(f"Falling back to simple model for {model_spec}")
            return create_simple_model_for_builtin(model_spec)

    # Try loading from file
    if model_path.exists():
        print(f"Loading velocity model from file: {model_spec}")

        # Try .tvel format (travel time velocity format)
        if model_spec.endswith('.tvel'):
            try:
                return load_tvel_model(model_path)
            except Exception as e:
                print(f"Error loading .tvel format: {e}")
                raise

        # Try ND format
        if model_spec.endswith('.nd') or model_spec.endswith('.m'):
            try:
                return cake.load_model(model_spec)
            except Exception as e:
                print(f"Error loading ND format: {e}")
                raise

        # Try simple ASCII format: depth(km) vp(km/s) vs(km/s) rho(g/cm3) qp qs
        try:
            return load_simple_velocity_model(model_path)
        except Exception as e:
            print(f"Error loading velocity model: {e}")
            raise

    raise ValueError(f"Unknown model specification: {model_spec}")


def load_tvel_model(filepath):
    """
    Load velocity model from .tvel format (travel time velocity format)

    Format (.tvel files typically used with NonLinLoc/TauP):
    # Comments start with #
    depth(km) vp(km/s) vs(km/s) rho(g/cm3)
    OR
    depth(km) vp(km/s) vs(km/s)  (rho calculated from vp)
    """
    model = cake.LayeredModel()

    depths = []
    materials = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            # Skip non-numeric lines (headers)
            try:
                depth_km = float(parts[0])
                vp_km_s = float(parts[1])
                vs_km_s = float(parts[2])
            except ValueError:
                continue  # Skip header lines

            # Density: either provided or estimated from Vp using Nafe-Drake
            if len(parts) >= 4:
                rho_g_cm3 = float(parts[3])
            else:
                # Nafe-Drake relation for density from Vp
                rho_g_cm3 = 1.6612 * vp_km_s - 0.4721 * vp_km_s**2 + \
                           0.0671 * vp_km_s**3 - 0.0043 * vp_km_s**4 + \
                           0.000106 * vp_km_s**5
                rho_g_cm3 = max(rho_g_cm3, 1.02)  # Minimum density

            # Quality factors - use typical crustal/mantle values
            if depth_km < 35:  # Crust
                qp, qs = 600., 300.
            else:  # Mantle
                qp, qs = 1400., 600.

            # Skip if this is a duplicate depth (velocity discontinuity)
            # We'll use the material from below the discontinuity
            if depths and abs(depth_km - depths[-1]) < 0.01:
                # Update the last material with the new values (below discontinuity)
                materials[-1] = cake.Material(
                    vp=vp_km_s * 1000.,
                    vs=vs_km_s * 1000.,
                    rho=rho_g_cm3 * 1000.,
                    qp=qp,
                    qs=qs
                )
            else:
                depths.append(depth_km)
                materials.append(cake.Material(
                    vp=vp_km_s * 1000.,    # Convert to m/s
                    vs=vs_km_s * 1000.,
                    rho=rho_g_cm3 * 1000., # Convert to kg/m3
                    qp=qp,
                    qs=qs
                ))

    # Build layered model from discontinuities
    # Limit to upper 60 km for efficiency (we don't need whole mantle for regional GFs)
    max_depth_km = 60
    for i in range(len(depths) - 1):
        if depths[i] >= max_depth_km:
            break
        zbot = min(depths[i+1] * 1000., max_depth_km * 1000.)
        ztop = depths[i] * 1000.

        # Ensure layer has some thickness (avoid zero thickness)
        if abs(zbot - ztop) < 1.0:  # Less than 1 meter
            continue

        model.append(cake.HomogeneousLayer(
            ztop=ztop,
            zbot=zbot,
            m=materials[i]
        ))

    # Last layer extends to 80 km depth
    if depths and len(depths) > 0:
        last_depth_m = min(depths[-1] * 1000., max_depth_km * 1000.)
        if last_depth_m < 80000:
            model.append(cake.HomogeneousLayer(
                ztop=last_depth_m,
                zbot=80000.,  # 80 km
                m=materials[-1]
            ))

    return model


def create_simple_model_for_builtin(model_name):
    """Create simplified version of standard Earth models"""
    model = cake.LayeredModel()

    if model_name == 'iasp91':
        # Simplified IASP91
        model.append(cake.HomogeneousLayer(
            ztop=0., zbot=20000.,
            m=cake.Material(vp=5800., vs=3360., rho=2720., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=20000., zbot=35000.,
            m=cake.Material(vp=6500., vs=3750., rho=2920., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=35000., zbot=80000.,
            m=cake.Material(vp=8040., vs=4470., rho=3320., qp=1400., qs=600.)
        ))
    elif model_name.startswith('ak135'):
        # Simplified AK135
        model.append(cake.HomogeneousLayer(
            ztop=0., zbot=20000.,
            m=cake.Material(vp=5800., vs=3460., rho=2720., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=20000., zbot=35000.,
            m=cake.Material(vp=6800., vs=3900., rho=2920., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=35000., zbot=80000.,
            m=cake.Material(vp=8110., vs=4491., rho=3371., qp=1400., qs=600.)
        ))
    else:  # PREM-like
        model.append(cake.HomogeneousLayer(
            ztop=0., zbot=15000.,
            m=cake.Material(vp=6800., vs=3900., rho=2900., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=15000., zbot=24400.,
            m=cake.Material(vp=6800., vs=3900., rho=2900., qp=600., qs=300.)
        ))
        model.append(cake.HomogeneousLayer(
            ztop=24400., zbot=80000.,
            m=cake.Material(vp=8110., vs=4500., rho=3380., qp=1400., qs=600.)
        ))

    return model


def load_simple_velocity_model(filepath):
    """
    Load velocity model from simple ASCII format

    Format (one layer per line):
    # Comments start with #
    depth_top(km) depth_bot(km) vp(km/s) vs(km/s) rho(g/cm3) qp qs

    Use depth_bot=0 for half-space

    Example:
    0.0   2.0  5.5  3.2  2.6  600  300
    2.0   10.0 6.2  3.6  2.8  600  300
    10.0  30.0 6.8  3.9  2.9  600  300
    30.0  0    8.1  4.5  3.3  1400 600  # Half-space
    """
    model = cake.LayeredModel()

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            ztop, zbot, vp, vs, rho, qp, qs = map(float, parts[:7])

            # Convert to meters and kg/m^3
            ztop_m = ztop * 1000.
            zbot_m = zbot * 1000. if zbot > 0 else None

            mat = cake.Material(
                vp=vp * 1000.,    # km/s -> m/s
                vs=vs * 1000.,
                rho=rho * 1000.,  # g/cm3 -> kg/m3
                qp=qp,
                qs=qs
            )

            if zbot_m is None or zbot_m > 100000.:
                # Half-space (limit to 100 km depth)
                zbot_m = 100000.

            model.append(cake.HomogeneousLayer(
                ztop=ztop_m,
                zbot=zbot_m,
                m=mat
            ))

    return model


def create_config_template(filename='gf_config_template.yaml'):
    """Create a template configuration file"""
    template = {
        'store_name': 'my_region_8hz',
        'earth_model': 'iasp91',  # or path to model file
        'sample_rate': 8.0,
        'distance': {
            'min': 10.0,
            'max': 1000.0,
            'delta': 10.0
        },
        'depth': {
            'min': 1.0,
            'max': 50.0,
            'delta': 1.0
        },
        'duration': 300.0,
        'modelling_code': 'qseis.2006a',
        'qseis': {
            'sw_algorithm': 1,
            'sw_flat_earth_transform': 1,
            'slowness_window': [0.0, 0.0, 0.5, 1.5]
        }
    }

    with open(filename, 'w') as f:
        yaml.dump(template, f, default_flow_style=False, sort_keys=False)

    print(f"Configuration template created: {filename}")
    print("\nEdit this file and run:")
    print(f"  python3 generate_gf_custom.py --config {filename}")


def create_gf_store_from_config(config, base_dir='gf_stores'):
    """
    Create GF store from configuration dictionary

    Parameters:
    -----------
    config : dict
        Configuration dictionary
    base_dir : str
        Base directory for stores
    """
    store_name = config['store_name']
    print("=" * 70)
    print(f"Creating Green's Function Store: {store_name}")
    print("=" * 70)

    # Load earth model
    earth_model = load_velocity_model(config['earth_model'])

    # Extract parameters
    sample_rate = config['sample_rate']
    dist_min = config['distance']['min'] * 1000.  # Convert to meters
    dist_max = config['distance']['max'] * 1000.
    dist_delta = config['distance']['delta'] * 1000.
    depth_min = config['depth']['min'] * 1000.
    depth_max = config['depth']['max'] * 1000.
    depth_delta = config['depth']['delta'] * 1000.
    duration = config.get('duration', 300.0)

    # Calculate statistics
    n_dist = int((dist_max - dist_min) / dist_delta) + 1
    n_depth = int((depth_max - depth_min) / depth_delta) + 1
    total_gfs = n_dist * n_depth

    print(f"Earth model: {config['earth_model']}")
    print(f"Sample rate: {sample_rate} Hz")
    print(f"Distance: {config['distance']['min']}-{config['distance']['max']} km " +
          f"(step: {config['distance']['delta']} km, n={n_dist})")
    print(f"Depth: {config['depth']['min']}-{config['depth']['max']} km " +
          f"(step: {config['depth']['delta']} km, n={n_depth})")
    print(f"Duration: {duration} s")
    print(f"Total GFs to compute: {total_gfs:,}")
    print("=" * 70)

    # Create store path
    store_path = Path(base_dir) / store_name
    store_path.parent.mkdir(exist_ok=True, parents=True)

    if store_path.exists():
        print(f"\nWarning: Store exists at {store_path}")
        response = input("Delete and recreate? (yes/no): ")
        if response.lower() != 'yes':
            print("Aborted.")
            return None
        import shutil
        shutil.rmtree(store_path)

    # Create QSEIS configuration
    qsconf = qseis.QSeisConfig()
    qsconf.qseis_version = '2006b'  # Using 2006b (most commonly available)
    qsconf.time_region = (gf.Timing('0'), gf.Timing(f'{duration}'))
    qsconf.cut = (gf.Timing('0'), gf.Timing(f'{duration}'))
    qsconf.wavelet_duration_samples = 0.001

    # Apply QSEIS-specific settings
    qseis_config = config.get('qseis', {})
    qsconf.sw_algorithm = qseis_config.get('sw_algorithm', 1)
    qsconf.sw_flat_earth_transform = qseis_config.get('sw_flat_earth_transform', 1)

    slowness = qseis_config.get('slowness_window', [0.0, 0.0, 0.5, 1.5])
    qsconf.slowness_window = tuple(slowness)

    # Create store configuration
    modelling_code = config.get('modelling_code', 'qseis.2006b')

    store_config = gf.ConfigTypeA(
        id=store_name,
        ncomponents=10,  # Full moment tensor
        sample_rate=sample_rate,
        receiver_depth=0.0,
        source_depth_min=depth_min,
        source_depth_max=depth_max,
        source_depth_delta=depth_delta,
        distance_min=dist_min,
        distance_max=dist_max,
        distance_delta=dist_delta,
        earthmodel_1d=earth_model,
        modelling_code_id=modelling_code,
    )

    # Validate and create
    print("\nValidating configuration...")
    store_config.validate()

    print(f"Creating store at: {store_path}")
    gf.store.Store.create(str(store_path), config=store_config, extra={'qseis': qsconf})

    print("Creating travel time tables...")
    store = gf.store.Store(str(store_path), 'r')
    store.make_ttt()
    store.close()

    print("\n" + "=" * 70)
    print("Store created successfully!")
    print(f"Location: {store_path.absolute()}")
    print("=" * 70)
    print("\nNext steps:")
    print(f"1. Build GFs: cd {store_path} && fomosto build --nworkers=8")
    print(f"2. Convert: python3 convert_to_sc3gf1d.py {store_path} gf_sc3gf1d/{store_name}")
    print("=" * 70)

    return store_path


def main():
    parser = argparse.ArgumentParser(
        description='Generate Green\'s Functions with custom configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--config', help='Path to YAML configuration file')
    parser.add_argument('--create-template', metavar='FILENAME',
                       help='Create a configuration template file')
    parser.add_argument('--model', default='iasp91',
                       help='Earth model (iasp91, ak135-f-average, etc.) or path to model file')
    parser.add_argument('--store-name', help='Name for the GF store')
    parser.add_argument('--sample-rate', type=float, default=8.0,
                       help='Sampling rate in Hz')
    parser.add_argument('--distance-max', type=float, default=1000.0,
                       help='Maximum distance in km')
    parser.add_argument('--depth-max', type=float, default=50.0,
                       help='Maximum depth in km')

    args = parser.parse_args()

    # Create template
    if args.create_template:
        create_config_template(args.create_template)
        return 0

    # Load or create configuration
    if args.config:
        print(f"Loading configuration from: {args.config}")
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    else:
        # Create config from command line arguments
        store_name = args.store_name or f"custom_{int(args.sample_rate)}hz"
        config = {
            'store_name': store_name,
            'earth_model': args.model,
            'sample_rate': args.sample_rate,
            'distance': {'min': 10.0, 'max': args.distance_max, 'delta': 10.0},
            'depth': {'min': 1.0, 'max': args.depth_max, 'delta': 1.0},
            'duration': 300.0,
        }

    # Create store
    store_path = create_gf_store_from_config(config)

    return 0 if store_path else 1


if __name__ == "__main__":
    sys.exit(main())
