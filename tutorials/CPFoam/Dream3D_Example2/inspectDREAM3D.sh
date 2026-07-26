#!/bin/bash
# Inspect DREAM3D file structure

if [ $# -eq 0 ]; then
    echo "Usage: $0 <dream3d_file>"
    echo "Example: $0 constant/microstructure.dream3d"
    exit 1
fi

DREAM3D_FILE="$1"

if [ ! -f "$DREAM3D_FILE" ]; then
    echo "Error: File not found: $DREAM3D_FILE"
    exit 1
fi

echo "========================================="
echo "DREAM3D File Structure Inspector"
echo "========================================="
echo "File: $DREAM3D_FILE"
echo ""

echo "Full file structure:"
echo "-------------------------------------"
h5dump -n "$DREAM3D_FILE"
echo ""

echo "========================================="
echo "Looking for data containers..."
echo "-------------------------------------"
h5ls -r "$DREAM3D_FILE" | grep "DataContainers"
echo ""

echo "========================================="
echo "Looking for geometry information..."
echo "-------------------------------------"
h5ls -r "$DREAM3D_FILE" | grep -i "geometry\|dimension\|origin\|spacing\|resolution"
echo ""

echo "========================================="
echo "Looking for cell data..."
echo "-------------------------------------"
h5ls -r "$DREAM3D_FILE" | grep -i "celldata\|featureid\|euler\|phase"
echo ""

echo "========================================="
echo "Data Container Details:"
echo "-------------------------------------"
# Try to find the actual data container
CONTAINERS=$(h5ls "$DREAM3D_FILE/DataContainers" 2>/dev/null || echo "")
if [ -n "$CONTAINERS" ]; then
    echo "Found containers:"
    echo "$CONTAINERS"
    echo ""
    
    # For each container, show its structure
    for container in $CONTAINERS; do
        container=${container%/}  # Remove trailing slash
        echo "Container: $container"
        h5ls -r "$DREAM3D_FILE/DataContainers/$container" 2>/dev/null | head -30
        echo ""
    done
else
    echo "No DataContainers group found"
fi

echo "========================================="
echo "Recommendations:"
echo "-------------------------------------"
echo "1. Check the 'Data Container Details' section above"
echo "2. Find where DIMENSIONS, ORIGIN, SPACING/RESOLUTION are located"
echo "3. Update the utility with correct paths"
echo ""
echo "Common alternatives:"
echo "  - RESOLUTION vs SPACING"
echo "  - _SIMPL_GEOMETRY vs GEOMETRY"
echo "  - Different container names"
echo ""
