#!/bin/bash
# with_env.sh TOOL COMMAND... -- run COMMAND inside TOOL's configured environment.
#
#   TOOL is QUPKAKE, DB2C, PREP, DOCK or PKA (see config/hpc.env). Examples:
#
#   common/with_env.sh PREP  python run_qupkake/prepare_parents.py samples/H*.smi -o parents
#   common/with_env.sh PREP  python run_qupkake/assign_names.py --shards-dir protomers ...
#   common/with_env.sh DB2C  python run_qupkake/check_db2_stereo.py library/library.smi
#   common/with_env.sh DB2C  python run_qupkake/verify_db2_build.py library/library.tsv ...
#   common/with_env.sh DOCK  python3 run_qupkake/db2_to_mol2.py build/ -o library.mol2 ...
#
# Also exported when set: DB2C_SRC (read by check_db2_stereo.py),
# BUILD_LIGAND_EXE, DB2TOOL_EXE (read by db2_to_mol2.py) and MOLGPKA_DIR.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
tool="$(tr '[:lower:]' '[:upper:]' <<<"$1")"; shift
case "$tool" in QUPKAKE|DB2C|PREP|DOCK|PKA) ;; *) echo "unknown TOOL '$tool' (QUPKAKE, DB2C, PREP, DOCK, PKA)" >&2; exit 64 ;; esac

LLB_ROOT="${LLB_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config
llb_activate "$tool"
[[ -n "${DB2C_SRC:-}" ]] && export DB2C_SRC
[[ -n "${BUILD_LIGAND_EXE:-}" ]] && export BUILD_LIGAND_EXE
[[ -n "${DB2TOOL_EXE:-}" ]] && export DB2TOOL_EXE
[[ -n "${MOLGPKA_DIR:-}" ]] && export MOLGPKA_DIR
exec "$@"
