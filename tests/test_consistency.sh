#!/bin/bash
# test_consistency.sh
# Automates the cross-checking of different k2tree formats

set -e

BUILD_DIR=${1:-"."}
ROOT_DIR=${2:-".."}

# Define paths to executables
# On Windows/MinGW, suffix might be .x.exe
K2SPARSE="$BUILD_DIR/k2sparse.x"
K2CPDF="$BUILD_DIR/k2cpdf.x"
B128SPARSE="$BUILD_DIR/b128sparse.x"
MATRIXCMP="$BUILD_DIR/test_matrixcmp.x"

# Compatibility check for Windows (.exe suffix)
if [[ ! -f "$K2SPARSE" && -f "$K2SPARSE.exe" ]]; then
    K2SPARSE="$K2SPARSE.exe"
    K2CPDF="$K2CPDF.exe"
    B128SPARSE="$B128SPARSE.exe"
    MATRIXCMP="$MATRIXCMP.exe"
fi

# Create a temporary directory for test files
TEMP_DIR="tests/temp_consistency"
mkdir -p "$TEMP_DIR"

# Matrix size for B128
DIM=16

echo "--- Generating test matrix ($DIM x $DIM sparse) ---"
cat <<EOF > "$INPUT_TXT"
0 0
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
9 9
10 10
11 11
12 12
13 13
14 14
15 15
0 15
15 0
5 10
3 3
3 4
4 3
4 4
10 10
10 11
11 10
11 11
EOF

echo -e "\n--- Testing Standard K2-TREE (.k2) ---"
"$K2SPARSE" -o "$TEST_K2" "$INPUT_TXT"
"$K2SPARSE" -d -o "$DECODED_TXT" "$TEST_K2"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: K2-TREE matches input."

echo -e "\n--- Testing Compressed K2-DFS (.ck2) ---"
"$K2CPDF" -o "$TEST_CK2" -t 16 "$TEST_K2"
BACKP="${TEST_CK2}.p"
"$K2SPARSE" -d -I "$BACKP" -o "$DECODED_TXT" "$TEST_CK2"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: K2-DFS (Compressed) matches input."

echo -e "\n--- Testing B128 Baseline (.b128) ---"
"$B128SPARSE" -o "$TEST_B128" -s "$DIM" "$INPUT_TXT"
"$B128SPARSE" -d -o "$DECODED_TXT" "$TEST_B128"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: B128 matches input."

echo -e "\nALL FORMATS ARE CONSISTENT!"

# Cleanup
# rm -rf "$TEMP_DIR"
