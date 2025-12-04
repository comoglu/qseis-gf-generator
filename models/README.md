# Velocity Models

This directory contains velocity models in `.tvel` format for use with the GF generator.

## Included Models

### iasp91.tvel (in parent directory)
- **Region**: Global model, suitable for Australia
- **Description**: IASPEI91 velocity model
- **Depth coverage**: 0-60 km (crust and upper mantle)
- **Use case**: Regional networks, continental crust

### turkey_ak135_deep.tvel
- **Region**: Turkey and Mediterranean
- **Description**: AK135 model with deep coverage for subduction zones
- **Depth coverage**: 0-187.5 km (includes intermediate depth earthquakes)
- **Use case**: Subduction zones (Aegean, Mediterranean), deep earthquakes

## File Format

The `.tvel` format is simple ASCII:

```
# Comments start with #
# Format: depth(km) vp(km/s) vs(km/s) rho(g/cm3)
    0.000    5.8000    3.3600    2.7200
   20.000    5.8000    3.3600    2.7200
   20.000    6.5000    3.7500    2.9200  # Velocity discontinuity
   35.000    6.5000    3.7500    2.9200
   35.000    8.0400    4.4700    3.3198  # Moho
```

**Features:**
- Lines starting with `#` are comments
- Duplicate depths represent velocity discontinuities (e.g., Moho, Conrad)
- Density (4th column) is optional - will be calculated from Vp if not provided
- Model is automatically extended to 80 km depth if shallower

## Creating Custom Models

### Option 1: Manual Creation

Create a `.tvel` file with your velocity structure:

```bash
cat > my_model.tvel <<EOF
# My custom velocity model
# depth(km) vp(km/s) vs(km/s) rho(g/cm3)
0.0    5.5    3.2    2.6
10.0   5.5    3.2    2.6
10.0   6.0    3.5    2.7
30.0   6.0    3.5    2.7
30.0   8.0    4.5    3.3
60.0   8.0    4.5    3.3
EOF
```

### Option 2: Extract from Existing Fomosto Store

If you have an existing Fomosto store:

```bash
python3 ../extract_model_from_fomosto.py \
    /path/to/fomosto_store \
    my_extracted_model.tvel
```

### Option 3: Convert from Other Formats

For ND format (Pyrocko/cake format):
```bash
# Use Pyrocko's built-in models directly in generate_gf_custom.py:
python3 ../generate_gf_custom.py --model ak135-f-average ...
```

## Usage

Use any `.tvel` model with the generator:

```bash
python3 ../generate_gf_custom.py \
    --model models/turkey_ak135_deep.tvel \
    --store-name my_store \
    ...
```

## Model Selection Guidelines

**For shallow crustal earthquakes (< 40 km):**
- Use IASP91 or AK135 continental/average models
- Coverage: 0-60 km is sufficient

**For subduction zones (40-200 km):**
- Use deep models like `turkey_ak135_deep.tvel`
- Coverage: 0-200 km recommended
- Include mantle velocities

**For local studies:**
- Create custom models from local velocity studies
- Use tomography results if available
- Ensure model extends at least 20 km deeper than max source depth

## Quality Factors

The generator uses these Q values:
- **Crust** (< 35 km): Qp = 600, Qs = 300
- **Mantle** (> 35 km): Qp = 1400, Qs = 600

These are typical values and work well for most applications.

## References

- IASP91: Kennett & Engdahl (1991), Geophys. J. Int.
- AK135: Kennett et al. (1995), Geophys. J. Int.
- PREM: Dziewonski & Anderson (1981), Phys. Earth Planet. Inter.
