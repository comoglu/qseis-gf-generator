#!/bin/bash
# Example: Generate a small test Green's Function store

# Generate test store (1-10 km depth, 10-50 km distance)
python3 ../generate_gf_custom.py \
    --model ../iasp91.tvel \
    --store-name test_8hz_iasp91 \
    --sample-rate 8 \
    --depth-min 1 --depth-max 10 --depth-delta 1 \
    --distance-min 10 --distance-max 50 --distance-delta 10 \
    --nworkers 2

# Convert to sc3gf1d format
python3 ../convert_to_sc3gf1d.py \
    ../gf_stores/test_8hz_iasp91 \
    ../gf_sc3gf1d/test_8hz_iasp91

echo "Test store generated successfully!"
echo "Location: ../gf_sc3gf1d/test_8hz_iasp91/"
