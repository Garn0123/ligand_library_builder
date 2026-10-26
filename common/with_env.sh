#!/bin/bash
# with_env.sh TOOL COMMAND... -- run COMMAND inside TOOL's configured environment.
#
#   TOOL is PREP or DOCK (see config/hpc.env.example). Examples:
#
#   common/with_env.sh PREP  python docking/mol2_library.py library.mol2 --library library.tsv -o mol2_lib
#   common/with_env.sh DOCK  python3 docking/db2_to_mol2.py build/ -o library.mol2 ...
#
# Also exported when set: DB2TOOL_EXE (read by db2_to_mol2.py).

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
tool="$(tr '[:lower:]' '[:upper:]' <<<"$1")"; shift
case "$tool" in PREP|DOCK) ;; *) echo "unknown TOOL '$tool' (PREP, DOCK)" >&2; exit 64 ;; esac

LLB_ROOT="${LLB_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config
[[ -n "${DB2TOOL_EXE:-}" ]] && export DB2TOOL_EXE
llb_activate "$tool"
exec "$@"
