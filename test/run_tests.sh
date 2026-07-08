#!/bin/bash
# Test script for evalAdmix: validates old method, new corrected method,
# and -useInds functionality. Assumes the binary is already built.
#
# Usage: bash test/run_tests.sh   (run from the evalAdmix root directory)
#   or:  bash run_tests.sh        (run from inside test/)

set -e

# Resolve paths
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT_DIR/evalAdmix"
DATA="$ROOT_DIR/data"
TMP="$SCRIPT_DIR/tmp"
mkdir -p "$TMP"

PASS=0
FAIL=0

check_pass() { echo "[PASS] $1"; PASS=$((PASS+1)); }
check_fail() { echo "[FAIL] $1"; FAIL=$((FAIL+1)); }

# --- Helper: check matrix dimensions, NA diagonal ---
validate_matrix() {
  local f="$1"
  local expected_n="$2"
  local label="$3"

  local ncols
  ncols=$(head -1 "$f" | awk '{print NF}')
  if [ "$ncols" -ne "$expected_n" ]; then
    check_fail "$label: expected $expected_n columns, got $ncols"
    return 1
  fi
  local nrows
  nrows=$(wc -l < "$f" | tr -d ' ')
  if [ "$nrows" -ne "$expected_n" ]; then
    check_fail "$label: expected $expected_n rows, got $nrows"
    return 1
  fi
  # Check diagonal is NA
  local diag_ok
  diag_ok=$(awk '{if($NR!="NA"){print "BAD"; exit}}' "$f")
  if [ -n "$diag_ok" ]; then
    check_fail "$label: diagonal is not all NA"
    return 1
  fi
  check_pass "$label: ${nrows}x${ncols} matrix, NA diagonal OK"
}

# --- Helper: compare two matrices (exact) ---
compare_exact() {
  local a="$1"; local b="$2"; local label="$3"
  if diff -q "$a" "$b" > /dev/null 2>&1; then
    check_pass "$label: output matches reference exactly"
  else
    check_fail "$label: output differs from reference"
  fi
}

echo "=== evalAdmix test suite ==="
echo ""

# ----------------------------------------------------------
# Test 1: Old method (default) — regression vs reference
# ----------------------------------------------------------
echo "--- Test 1: Old method (default), -P 1 ---"
OUT1="$TMP/test1_old.corres"
"$BIN" -plink "$DATA/admixTjeck2" \
  -fname "$DATA/admixTjeck2.3.P" \
  -qname "$DATA/admixTjeck2.3.Q" \
  -P 1 -o "$OUT1" 2>"$TMP/test1.log"
validate_matrix "$OUT1" 374 "Test1 old method"
compare_exact "$OUT1" "$DATA/old_freq_correction.corres" "Test1 old method"
echo ""

# ----------------------------------------------------------
# Test 2: New method (corrected) — regression vs reference
# ----------------------------------------------------------
echo "--- Test 2: New method (-method corrected) ---"
OUT2="$TMP/test2_corrected.corres"
"$BIN" -plink "$DATA/admixTjeck2" \
  -fname "$DATA/admixTjeck2.3.P" \
  -qname "$DATA/admixTjeck2.3.Q" \
  -method corrected -o "$OUT2" 2>"$TMP/test2.log"
validate_matrix "$OUT2" 374 "Test2 corrected method"
compare_exact "$OUT2" "$DATA/new_freq_correction_v2.corres" "Test2 corrected method"
echo ""

# ----------------------------------------------------------
# Generate useInds file: exclude all YRI (108 excluded → 266 used)
# ----------------------------------------------------------
USEINDS_FILE="$TMP/useInds_excludeYRI.txt"
awk '{print $1"\t"($2=="YRI"?0:1)}' "$DATA/admixTjeck2.fam" > "$USEINDS_FILE"
N_INCLUDED=$(awk -F'\t' '{s+=$2}END{print s}' "$USEINDS_FILE")
echo "useInds file: $USEINDS_FILE ($N_INCLUDED individuals included, YRI excluded)"
echo ""

# ----------------------------------------------------------
# Test 3: Old method + useInds
# ----------------------------------------------------------
echo "--- Test 3: Old method + useInds (exclude YRI) ---"
OUT3="$TMP/test3_old_useInds.corres"
"$BIN" -plink "$DATA/admixTjeck2" \
  -fname "$DATA/admixTjeck2.3.P" \
  -qname "$DATA/admixTjeck2.3.Q" \
  -P 1 -useInds "$USEINDS_FILE" -o "$OUT3" 2>"$TMP/test3.log"
validate_matrix "$OUT3" "$N_INCLUDED" "Test3 old+useInds"
echo ""

# ----------------------------------------------------------
# Test 4: New method + useInds
# ----------------------------------------------------------
echo "--- Test 4: New method + useInds (exclude YRI) ---"
OUT4="$TMP/test4_corrected_useInds.corres"
"$BIN" -plink "$DATA/admixTjeck2" \
  -fname "$DATA/admixTjeck2.3.P" \
  -qname "$DATA/admixTjeck2.3.Q" \
  -method corrected -useInds "$USEINDS_FILE" -o "$OUT4" 2>"$TMP/test4.log"
validate_matrix "$OUT4" "$N_INCLUDED" "Test4 corrected+useInds"
echo ""

# ----------------------------------------------------------
# Test 5: Cross-check old vs new with useInds
# (Expect near-agreement, similar to full-data comparison)
# ----------------------------------------------------------
echo "--- Test 5: Cross-check old vs new (useInds) ---"
if python3 - "$OUT3" "$OUT4" <<'PYEOF' 2>/dev/null; then
  import math, sys
  a=[[float(x) if x!="NA" else float("nan") for x in l.split()] for l in open(sys.argv[1])]
  b=[[float(x) if x!="NA" else float("nan") for x in l.split()] for l in open(sys.argv[2])]
  n=len(a)
  ds=[b[i][j]-a[i][j] for i in range(n) for j in range(n)
      if i!=j and math.isfinite(a[i][j]) and math.isfinite(b[i][j])]
  if not ds:
      print("ERROR: no comparable pairs"); sys.exit(1)
  mx=max(abs(x) for x in ds)
  mean=abs(sum(abs(x) for x in ds)/len(ds))
  print(f"  pairs={len(ds)} max_abs_diff={mx:.6f} mean_abs_diff={mean:.6f}")
  if mx < 0.01:
      sys.exit(0)
  else:
      print("  max abs diff > 0.01 — methods disagree more than expected")
      sys.exit(1)
PYEOF
  check_pass "Test5 cross-check: old vs new agree (max abs diff < 0.01)"
else
  check_fail "Test5 cross-check: old vs new disagree or python3 unavailable"
fi
echo ""

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------
echo "=== Summary: $PASS passed, $FAIL failed ==="
if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
exit 0