# test_consistency.ps1
# Automates the cross-checking of different k2tree formats

param (
    [string]$BuildDir = ".",
    [string]$RootDir = ".."
)

$ErrorActionPreference = "Stop"

# Define paths to executables
$k2sparse = Join-Path $BuildDir "k2sparse.x.exe"
$k2cpdf = Join-Path $BuildDir "k2cpdf.x.exe"
$b128sparse = Join-Path $BuildDir "b128sparse.x.exe"
$matrixcmp = Join-Path $BuildDir "matrixcmp.x.exe"

# Create a temporary directory for test files
$tempDir = Join-Path $PSScriptRoot "temp_consistency"
if (Test-Path $tempDir) { Remove-Item -Recurse -Force $tempDir }
New-Item -ItemType Directory -Path $tempDir | Out-Null

$inputTxt = Join-Path $tempDir "input.txt"
$testK2 = Join-Path $tempDir "test.k2"
$testCK2 = Join-Path $tempDir "test.ck2"
$testB128 = Join-Path $tempDir "test.b128"
$decodedTxt = Join-Path $tempDir "decoded.txt"

# Matrix size for B128
$dim = 16

Write-Host "--- Generating test matrix ($dim x $dim sparse) ---"
# Create a $dim x $dim sparse matrix as text
@"
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
"@ | Out-File -FilePath $inputTxt -Encoding ascii

Write-Host "`n--- Testing Standard K2-TREE (.k2) ---"
& $k2sparse -o $testK2 $inputTxt
& $k2sparse -d -o $decodedTxt $testK2
& $matrixcmp $inputTxt $decodedTxt
Write-Host "SUCCESS: K2-TREE matches input."

Write-Host "`n--- Testing Compressed K2-DFS (.ck2) ---"
# 1. Start from .k2
# 2. Compress to .ck2 (using k2cpdf)
# 3. Decompress .ck2 (using k2sparse -d -I)
# 4. Compare
& $k2cpdf -o $testCK2 -t 16 $testK2
$backp = $testCK2 + ".p"
& $k2sparse -d -I $backp -o $decodedTxt $testCK2
& $matrixcmp $inputTxt $decodedTxt
Write-Host "SUCCESS: K2-DFS (Compressed) matches input."

Write-Host "`n--- Testing B128 Baseline (.b128) ---"
& $b128sparse -o $testB128 -s $dim $inputTxt
& $b128sparse -d -o $decodedTxt $testB128
& $matrixcmp $inputTxt $decodedTxt
Write-Host "SUCCESS: B128 matches input."

Write-Host "`nALL FORMATS ARE CONSISTENT!"

# Cleanup
# Remove-Item -Recurse -Force $tempDir
