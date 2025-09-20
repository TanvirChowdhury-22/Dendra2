#!/usr/bin/env bash

set -euo pipefail
SANDER_BIN="${SANDER:-sander}"

die(){ echo "ERROR: $*" >&2; exit 1; }
have(){ command -v "$1" >/dev/null 2>&1; }
have "$SANDER_BIN" || die "sander not found (set SANDER=/path/to/sander)"

# Resolve paths relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPROT_DIR="$(dirname "$SCRIPT_DIR")"

PRMTOP="$DEPROT_DIR/solvated_complex_deprot.prmtop"
START_MINLOCAL="$DEPROT_DIR/min_local_solv_deprot.rst"
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

MIN1="$MD_DIR/tc_qmmm_deprot_min1.in"
MIN2="$MD_DIR/tc_qmmm_deprot_min2.in"
EQ1="$MD_DIR/tc_qmmm_eq1_nvt_heat_deprot.in"
EQ2="$MD_DIR/tc_qmmm_eq2_npt_dense_deprot.in"
EQ3="$MD_DIR/tc_qmmm_eq3_npt_hold_deprot.in"
PRN="$MD_DIR/tc_qmmm_prod_nvt_deprot.in"
PRP="$MD_DIR/tc_qmmm_prod_npt_deprot.in"

run_step "$START_MINLOCAL"                      "$MIN1" "$MD_DIR/tc_qmmm_deprot_min1"
run_step "$MD_DIR/tc_qmmm_deprot_min1.rst"     "$MIN2" "$MD_DIR/tc_qmmm_deprot_min2"
run_step "$MD_DIR/tc_qmmm_deprot_min2.rst"     "$EQ1"  "$MD_DIR/tc_qmmm_eq1_nvt_heat_deprot"
run_step "$MD_DIR/tc_qmmm_eq1_nvt_heat_deprot.rst" "$EQ2" "$MD_DIR/tc_qmmm_eq2_npt_dense_deprot"
run_step "$MD_DIR/tc_qmmm_eq2_npt_dense_deprot.rst" "$EQ3" "$MD_DIR/tc_qmmm_eq3_npt_hold_deprot"
run_step "$MD_DIR/tc_qmmm_eq3_npt_hold_deprot.rst"  "$PRN" "$MD_DIR/tc_qmmm_prod_nvt_deprot"
run_step "$MD_DIR/tc_qmmm_eq3_npt_hold_deprot.rst"  "$PRP" "$MD_DIR/tc_qmmm_prod_npt_deprot"
