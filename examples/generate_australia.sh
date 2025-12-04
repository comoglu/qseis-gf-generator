#!/bin/bash
# Example: Generate Green's Functions for Australia regional network

# Configuration
STORE_NAME="australia_8hz_iasp91"
SAMPLE_RATE=8
DEPTH_MIN=1
DEPTH_MAX=50
DEPTH_DELTA=1
DIST_MIN=10
DIST_MAX=2000
DIST_DELTA=10
NWORKERS=8

echo "Generating Green's Function store for Australia..."
echo "Store name: $STORE_NAME"
echo "Coverage: ${DEPTH_MIN}-${DEPTH_MAX} km depth, ${DIST_MIN}-${DIST_MAX} km distance"
echo "Workers: $NWORKERS"
echo ""

# Generate Fomosto store
python3 ../generate_gf_custom.py \
    --model ../iasp91.tvel \
    --store-name "$STORE_NAME" \
    --sample-rate "$SAMPLE_RATE" \
    --depth-min "$DEPTH_MIN" \
    --depth-max "$DEPTH_MAX" --depth-delta "$DEPTH_DELTA" \
    --distance-min "$DIST_MIN" --distance-max "$DIST_MAX" --distance-delta "$DIST_DELTA" \
    --nworkers "$NWORKERS"

if [ $? -eq 0 ]; then
    echo ""
    echo "Store generation complete. Converting to sc3gf1d format..."

    # Convert to sc3gf1d format
    python3 ../convert_to_sc3gf1d.py \
        "../gf_stores/$STORE_NAME" \
        "../gf_sc3gf1d/$STORE_NAME"

    if [ $? -eq 0 ]; then
        echo ""
        echo "SUCCESS! Green's Functions are ready:"
        echo "  Fomosto store: ../gf_stores/$STORE_NAME"
        echo "  sc3gf1d store: ../gf_sc3gf1d/$STORE_NAME"
        echo ""
        echo "To use in SeisComP:"
        echo "  cp -r ../gf_sc3gf1d/$STORE_NAME* ~/.seiscomp/share/gf/"
    else
        echo "ERROR: Conversion to sc3gf1d failed"
        exit 1
    fi
else
    echo "ERROR: Store generation failed"
    exit 1
fi
