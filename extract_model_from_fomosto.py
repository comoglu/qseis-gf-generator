#!/usr/bin/env python3
"""
Extract velocity model from Fomosto config and save as .tvel format
"""
import sys
from pathlib import Path

try:
    from pyrocko import gf
except ImportError:
    print("Error: Pyrocko not installed")
    sys.exit(1)

def extract_model_to_tvel(store_path, output_file):
    """
    Extract velocity model from Fomosto store and save as .tvel format

    Parameters:
    -----------
    store_path : str
        Path to Fomosto store
    output_file : str
        Output .tvel file path
    """
    # Load store config
    store = gf.store.Store(str(store_path), 'r')
    config = store.config

    # Get the earth model
    earthmodel = config.earthmodel_1d

    print(f"Extracting model from: {store_path}")
    print(f"Store ID: {config.id}")
    print(f"Sample rate: {config.sample_rate} Hz")

    # Write to .tvel format
    with open(output_file, 'w') as f:
        f.write(f"# Velocity model extracted from {config.id}\n")
        f.write(f"# Format: depth(km) vp(km/s) vs(km/s) rho(g/cm3)\n")

        # Iterate through layers
        from pyrocko import cake
        for layer in earthmodel.layers():
            ztop_km = layer.ztop / 1000.0
            zbot_km = layer.zbot / 1000.0

            # Handle both HomogeneousLayer and GradientLayer
            if isinstance(layer, cake.HomogeneousLayer):
                m = layer.m
                vp_km_s = m.vp / 1000.0
                vs_km_s = m.vs / 1000.0
                rho_g_cm3 = m.rho / 1000.0

                # Write top of layer
                f.write(f"{ztop_km:8.3f}  {vp_km_s:8.4f}  {vs_km_s:8.4f}  {rho_g_cm3:8.4f}\n")

                # Write bottom with same properties
                if zbot_km < 200:  # Write up to 200 km for subduction zones
                    f.write(f"{zbot_km:8.3f}  {vp_km_s:8.4f}  {vs_km_s:8.4f}  {rho_g_cm3:8.4f}\n")

            elif isinstance(layer, cake.GradientLayer):
                # For gradient layer, sample at top, middle, and bottom
                depths = [layer.ztop, (layer.ztop + layer.zbot) / 2, layer.zbot]
                for z in depths:
                    if z / 1000.0 > 200:  # Skip depths > 200 km
                        continue
                    m = layer.material(z)
                    z_km = z / 1000.0
                    vp_km_s = m.vp / 1000.0
                    vs_km_s = m.vs / 1000.0
                    rho_g_cm3 = m.rho / 1000.0
                    f.write(f"{z_km:8.3f}  {vp_km_s:8.4f}  {vs_km_s:8.4f}  {rho_g_cm3:8.4f}\n")

    store.close()

    print(f"\nModel saved to: {output_file}")
    print(f"\nYou can now use it with:")
    print(f"  python3 generate_gf_custom.py --model {output_file} --store-name <name> --sample-rate 8")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 extract_model_from_fomosto.py <store_path> <output.tvel>")
        print("\nExample:")
        print("  python3 extract_model_from_fomosto.py gf_koeri_qseis koeri_model.tvel")
        sys.exit(1)

    store_path = sys.argv[1]
    output_file = sys.argv[2]

    if not Path(store_path).exists():
        print(f"Error: Store not found: {store_path}")
        sys.exit(1)

    extract_model_to_tvel(store_path, output_file)

if __name__ == "__main__":
    main()
