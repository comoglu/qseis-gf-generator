# QSEIS Green's Function Generator

Generate Green's Functions (GF) using QSEIS via Pyrocko/Fomosto and convert them to SeisComP sc3gf1d format for moment tensor inversion.

## Overview

This toolkit provides:
- **Custom GF generation** with any velocity model (.tvel format)
- **Model extraction** from existing Fomosto stores
- **Format conversion** from Fomosto to SeisComP sc3gf1d (SAC files)

## Requirements

```bash
pip install pyrocko obspy numpy
```

You also need QSEIS 2006b installed via Fomosto:
```bash
# Install Fomosto QSEIS backend
fomosto init qseis.2006b
```

## Quick Start

### 1. Generate Green's Functions with Custom Model

```bash
python3 generate_gf_custom.py \
    --model iasp91.tvel \
    --store-name australia_8hz \
    --sample-rate 8 \
    --depth-min 1 --depth-max 50 --depth-delta 1 \
    --distance-min 10 --distance-max 2000 --distance-delta 10 \
    --nworkers 8
```

**Parameters:**
- `--model`: Velocity model file (.tvel format) or built-in name (iasp91, ak135-f-average, prem-no-ocean)
- `--store-name`: Output store name
- `--sample-rate`: Sampling frequency in Hz (e.g., 8 for 8 Hz)
- `--depth-min/max/delta`: Source depth range in km
- `--distance-min/max/delta`: Epicentral distance range in km
- `--nworkers`: Number of parallel workers (default: 4)

### 2. Convert to SeisComP sc3gf1d Format

```bash
python3 convert_to_sc3gf1d.py \
    gf_stores/australia_8hz \
    gf_sc3gf1d/australia_8hz
```

This creates:
- `gf_sc3gf1d/australia_8hz.desc` - Description file for SeisComP
- `gf_sc3gf1d/australia_8hz/<depth>/<dist>/<files>` - SAC files with GF components

**Output structure:**
```
gf_sc3gf1d/
├── australia_8hz.desc              # Description file
└── australia_8hz/
    ├── 0010/                        # Depth: 1.0 km (in 100m units)
    │   ├── 00010/                   # Distance: 10 km
    │   │   ├── 0010.00010.ZSS
    │   │   ├── 0010.00010.ZDD
    │   │   ├── 0010.00010.ZDS
    │   │   ├── 0010.00010.RSS
    │   │   ├── 0010.00010.RDD
    │   │   ├── 0010.00010.RDS
    │   │   ├── 0010.00010.TSS
    │   │   ├── 0010.00010.TDS
    │   │   ├── 0010.00010.ZEP      # For full moment tensor
    │   │   └── 0010.00010.REP      # For full moment tensor
    │   └── ...
    └── ...
```

### 3. Extract Velocity Model from Existing Store

```bash
python3 extract_model_from_fomosto.py \
    gf_stores/existing_store \
    my_model.tvel
```

## Velocity Model Format (.tvel)

The `.tvel` format is simple ASCII:

```
# Comments start with #
# Format: depth(km) vp(km/s) vs(km/s) rho(g/cm3)
    0.000    5.8000    3.3600    2.7200
   20.000    5.8000    3.3600    2.7200
   20.000    6.5000    3.7500    2.9200  # Discontinuity (duplicate depth)
   35.000    6.5000    3.7500    2.9200
   35.000    8.0400    4.4700    3.3198  # Moho
   60.000    8.0450    4.4850    3.3500
```

**Features:**
- Lines starting with `#` are comments
- Duplicate depths represent velocity discontinuities
- Density (4th column) is optional - calculated from Vp using Nafe-Drake relation if not provided
- Model is automatically extended to 80 km depth

**Included model:**
- `iasp91.tvel` - IASP91 global model (suitable for Australia)

## SeisComP Integration

### Using GFs in SeisComP

1. Copy the converted GF directory to SeisComP:
```bash
cp -r gf_sc3gf1d/australia_8hz* ~/.seiscomp/share/gf/ "or your custom GF folder"
```

2. Configure in `scmtv` or moment tensor tools:
```ini
# In scmtv.cfg or global.cfg
Or use scconfig to adjust related configuration
```

3. The `.desc` file format:
```
# depth [from] [to] [step]
depth 1 50 1
#
# distance [from] [to] [step]
distance 10 2000 10
#
# travel times
times LOCSAT iasp91
```

## Examples

### Regional Network (Australia)

Generate 8 Hz GFs for regional moment tensor inversion:

```bash
# Generate Fomosto store
python3 generate_gf_custom.py \
    --model iasp91.tvel \
    --store-name australia_8hz_iasp91 \
    --sample-rate 8 \
    --depth-min 1 --depth-max 50 --depth-delta 1 \
    --distance-min 10 --distance-max 2000 --distance-delta 10 \
    --nworkers 8

# Convert to sc3gf1d
python3 convert_to_sc3gf1d.py \
    gf_stores/australia_8hz_iasp91 \
    gf_sc3gf1d/australia_8hz_iasp91
```

### Subduction Zone (Turkey/Mediterranean)

Generate GFs with deep sources for subduction zone earthquakes:

```bash
# First extract and modify the AK135 model
python3 extract_model_from_fomosto.py \
    existing_store \
    turkey_ak135.tvel

# Generate with deeper sources
python3 generate_gf_custom.py \
    --model turkey_ak135.tvel \
    --store-name turkey_deep_8hz \
    --sample-rate 8 \
    --depth-min 1 --depth-max 200 --depth-delta 5 \
    --distance-min 10 --distance-max 1000 --distance-delta 10 \
    --nworkers 8

# Convert to sc3gf1d
python3 convert_to_sc3gf1d.py \
    gf_stores/turkey_deep_8hz \
    gf_sc3gf1d/turkey_deep_8hz
```

### Test/Development

Generate a small test store:

```bash
python3 generate_gf_custom.py \
    --model iasp91.tvel \
    --store-name test_8hz \
    --sample-rate 8 \
    --depth-min 1 --depth-max 10 --depth-delta 1 \
    --distance-min 10 --distance-max 50 --distance-delta 10 \
    --nworkers 2
```

## Technical Details

### GF Components

The tools generate 10 components for full moment tensor:
- **ZSS, ZDD, ZDS** - Vertical component
- **RSS, RDD, RDS** - Radial component
- **TSS, TDS** - Transverse component
- **ZEP, REP** - Explosion/isotropic components

For deviatoric moment tensor, only the first 8 components are needed.

### Performance

- **Parallel computation**: Use `--nworkers` to speed up generation
- **Memory usage**: ~100-200 MB per worker
- **Storage**: ~1-5 MB per depth/distance pair (depends on duration and sample rate)

**Typical generation times:**
- Test store (10 depths × 5 distances = 50 GFs): ~2-5 minutes with 2 workers
- Regional (50 depths × 200 distances = 10,000 GFs): ~2-8 hours with 8 workers
- Full coverage (200 depths × 100 distances = 20,000 GFs): ~8-24 hours with 8 workers

### QSEIS Configuration

The generator uses these QSEIS settings:
- **Algorithm**: Reflectivity method
- **Flat Earth transform**: Enabled
- **Quality factors**: Crustal (Qp=600, Qs=300), Mantle (Qp=1400, Qs=600)
- **Duration**: 300s (default) - automatically adjusted based on distance

## Troubleshooting

### "qseis executable not found"

Install QSEIS 2006b:
```bash
fomosto init qseis.2006b
```

### "Store validation failed"

Check that:
- Depth/distance ranges are positive
- Step sizes are not zero
- Model extends deep enough (at least to max source depth + 20 km)

### "Zero-thickness layer error"

The `.tvel` file has duplicate depths that create zero-thickness layers. The tool automatically handles velocity discontinuities - just ensure depths are properly formatted.

### Conversion fails with "No GF data"

Make sure the store was fully built:
```bash
cd gf_stores/your_store
fomosto check
```

## References

- **Pyrocko/Fomosto**: https://pyrocko.org/
- **QSEIS**: Wang, R. (2006) - Reflectivity method for computing synthetic seismograms
- **SeisComP**: https://www.seiscomp.de/
- **sc3gf1d format**: SeisComP Green's function format specification

## License

MIT License - Feel free to use and modify for your research.

## Author

Created for seismological research and moment tensor inversion workflows.

## Citation

If you use this toolkit in your research, please cite:
- Pyrocko: Heimann et al. (2017) - https://pyrocko.org
- QSEIS: Wang (2006) - The reflectivity method for computing synthetic seismograms
