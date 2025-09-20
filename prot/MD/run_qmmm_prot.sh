#!/usr/bin/env bash

set -euo pipefail
SANDER_BIN="${SANDER:-sander}"

die(){ echo "ERROR: $*" >&2; exit 1; }
have(){ command -v "$1" >/dev/null 2>&1; }
have "$SANDER_BIN" || die "sander not found (set SANDER=/path/to/sander)"

# Resolve paths relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT_DIR="$(dirname "$SCRIPT_DIR")"

PRMTOP="$PROT_DIR/solvated_complex_prot.prmtop"
START_MINLOCAL="$PROT_DIR/min_local_solv.rst"   
MD_DIR="$SCRIPT_DIR"

[[ -f "$PRMTOP" ]] || die "Missing $PRMTOP"
[[ -f "$START_MINLOCAL" ]] || die "Missing $START_MINLOCAL (run your local solvent min first)"

run_step () {
  local start="$1" mdin="$2" outbase="$3"
  [[ -f "$start" ]] || die "Missing coords: $start"
  [[ -f "$mdin"  ]] || die "Missing mdin:   $mdin"
  local out="${outbase}.out" rst="${outbase}.rst" nc="${outbase}.nc"
  if [[ -s "$rst" ]]; then echo "[SKIP] $mdin → $rst exists"; return; fi
  local ref_arg=""; [[ "$mdin" == *min1* ]] && ref_arg="-ref $start"
  local x_arg="";   [[ "$mdin" != *min1* && "$mdin" != *min2* ]] && x_arg="-x $nc"
  echo "[RUN ] $(basename "$mdin")"
  "$SANDER_BIN" -O -i "$mdin" -o "$out" -p "$PRMTOP" -c "$start" -r "$rst" $ref_arg $x_arg
  echo "[DONE] $outbase"
}

MIN1="$MD_DIR/tc_qmmm_prot_min1.in"
MIN2="$MD_DIR/tc_qmmm_prot_min2.in"
EQ1="$MD_DIR/tc_qmmm_eq1_nvt_heat_prot.in"
EQ2="$MD_DIR/tc_qmmm_eq2_npt_dense_prot.in"
EQ3="$MD_DIR/tc_qmmm_eq3_npt_hold_prot.in"
PRN="$MD_DIR/tc_qmmm_prod_nvt_prot.in"
PRP="$MD_DIR/tc_qmmm_prod_npt_prot.in"

run_step "$START_MINLOCAL"                   "$MIN1" "$MD_DIR/tc_qmmm_prot_min1"
run_step "$MD_DIR/tc_qmmm_prot_min1.rst"    "$MIN2" "$MD_DIR/tc_qmmm_prot_min2"
run_step "$MD_DIR/tc_qmmm_prot_min2.rst"    "$EQ1"  "$MD_DIR/tc_qmmm_eq1_nvt_heat_prot"
run_step "$MD_DIR/tc_qmmm_eq1_nvt_heat_prot.rst" "$EQ2" "$MD_DIR/tc_qmmm_eq2_npt_dense_prot"
run_step "$MD_DIR/tc_qmmm_eq2_npt_dense_prot.rst" "$EQ3" "$MD_DIR/tc_qmmm_eq3_npt_hold_prot"
run_step "$MD_DIR/tc_qmmm_eq3_npt_hold_prot.rst"  "$PRN" "$MD_DIR/tc_qmmm_prod_nvt_prot"
run_step "$MD_DIR/tc_qmmm_eq3_npt_hold_prot.rst"  "$PRP" "$MD_DIR/tc_qmmm_prod_npt_prot"
