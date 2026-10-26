#!/bin/bash
# smoke.sh -- does the toolchain image hold what llb.def says it does?
#
#   apptainer   exec llb.sif /opt/llb/smoke.sh [--allow-missing-bcl]
#   singularity exec llb.sif /opt/llb/smoke.sh [--allow-missing-bcl]
#
# Runs at the end of %post (build fails if it fails) and in CI under each
# runtime. Checks the image only: no LLB code, no ZBH tools. One line per check;
# exit 1 if any FAILs. Stands in for `llb doctor --no-zbh` until that exists.
#
# --allow-missing-bcl  a missing /opt/bcl is a WARN, not a FAIL (only for images
#                      built with ALLOW_MISSING_BCL; manifest.json says which)
set -uo pipefail

ALLOW_MISSING_BCL=0
[[ "${1:-}" == --allow-missing-bcl ]] && ALLOW_MISSING_BCL=1

export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-/opt/conda}
export XTBPATH=${XTBPATH:-/opt/xtb/bin/xtb}
export PYTHONNOUSERSITE=1
MM=/usr/local/bin/micromamba
nfail=0

emit() {  # STATUS NAME DETAIL
    printf '%-4s  %-22s %s\n' "$1" "$2" "$3"
    [[ "$1" == FAIL ]] && nfail=$((nfail + 1))
    return 0
}

check() {  # NAME CMD... : PASS with the last output line, or FAIL with the tail
    local name=$1 out rc
    shift
    out=$("$@" 2>&1); rc=$?
    if [[ $rc -eq 0 ]]; then
        emit PASS "$name" "$(tail -n1 <<<"$out")"
    else
        emit FAIL "$name" "rc=$rc: $(tail -n3 <<<"$out" | tr '\n' ' ')"
    fi
}

inenv() { local e=$1; shift; "$MM" run -n "$e" "$@"; }

# --- the three envs -------------------------------------------------------------
check env:qupkake inenv qupkake python -c '
import os, torch, rdkit, torch_geometric, pytorch_lightning, qupkake
assert qupkake.XTB_LOCATION == os.environ["XTBPATH"], qupkake.XTB_LOCATION
print(f"qupkake {qupkake.__version__}, torch {torch.__version__}, rdkit {rdkit.__version__}, "
      f"PyG {torch_geometric.__version__}, lightning {pytorch_lightning.__version__}, "
      f"xtb from {qupkake.XTB_LOCATION}")'

check env:db2_converter inenv db2_converter python -c '
import shutil, rdkit, posebusters, tabulate, openbabel
from db2_converter.config import config
assert rdkit.__version__ == "2024.03.5", rdkit.__version__
assert shutil.which("build_ligand"), "build_ligand not on PATH"
for k in ("BABEL_EXE", "ANTECHAMBER"):
    exe = config["all"][k]
    assert shutil.which(exe), f"{k}={exe} not on PATH"
print(f"rdkit {rdkit.__version__}, posebusters {posebusters.__version__}, build_ligand ok")'

check env:llb_pka inenv llb_pka python -c '
import torch, rdkit, torch_geometric, dimorphite_dl, numpy, pandas
assert rdkit.__version__ == "2024.03.5", rdkit.__version__
print(f"torch {torch.__version__}, rdkit {rdkit.__version__}, PyG {torch_geometric.__version__}")'

# --- conda stub (QupKake runs `conda list xtb` at import) --------------------------
check conda-stub inenv qupkake conda list xtb

# --- xtb ---------------------------------------------------------------------------
out=$("$XTBPATH" --version 2>&1 | grep -m1 'xtb version')
if [[ "$out" =~ xtb\ version\ 6\.4\.1 ]]; then
    emit PASS xtb "$XTBPATH: $(sed 's/^ *\* *//' <<<"$out")"
else
    emit FAIL xtb "$XTBPATH --version: '${out:-nothing}', want 6.4.1"
fi

# --- AMSOL (links libgfortran.so.5) -----------------------------------------------
amsol=$(inenv db2_converter python -c 'from db2_converter.config import config; print(config["all"]["AMSOLEXE"])' 2>/dev/null)
if [[ ! -x "$amsol" ]]; then
    emit FAIL amsol "AMSOLEXE '${amsol:-?}' is not executable"
elif ldd "$amsol" | grep -q 'not found'; then
    emit FAIL amsol "$amsol: $(ldd "$amsol" | grep 'not found' | tr -s ' \n' ' ')"
else
    emit PASS amsol "$amsol: every library found"
fi

# --- BCL ---------------------------------------------------------------------------
BCL_EXE=/opt/bcl/bin/bcl-apps-static.exe
if [[ ! -e /opt/bcl ]]; then
    if [[ $ALLOW_MISSING_BCL -eq 1 ]]; then
        emit WARN bcl "/opt/bcl absent (image built with ALLOW_MISSING_BCL); conformator_bcl will not run"
    else
        emit FAIL bcl "/opt/bcl absent: this image cannot run -m conformator_bcl"
    fi
else
    if [[ ! -x $BCL_EXE ]]; then
        emit FAIL bcl "$BCL_EXE missing or not executable"
    elif ldd "$BCL_EXE" 2>&1 | grep -q 'not found'; then
        emit FAIL bcl "$BCL_EXE: $(ldd "$BCL_EXE" | grep 'not found' | tr -s ' \n' ' ')"
    else
        # BCL's own flag syntax is `-help`; accept either spelling.
        if out=$("$BCL_EXE" --help 2>&1) || out=$("$BCL_EXE" -help 2>&1); then
            emit PASS bcl "$BCL_EXE runs ($(grep -m1 -io 'bcl[^ ]* *v*[0-9][0-9.]*' <<<"$out" || echo 'help printed'))"
        else
            emit FAIL bcl "$BCL_EXE --help / -help failed: $(tail -n2 <<<"$out" | tr '\n' ' ')"
        fi
    fi
    check bcl:db2_converter inenv db2_converter python -c '
import os
from db2_converter.config import config
base, cmd = config["bcl"]["BCLBASE"], config["bcl"]["BCL"]
exe = os.path.join(base, cmd.split()[0])
assert os.access(exe, os.X_OK), f"BCLBASE/BCL -> {exe} not executable"
print(f"BCLBASE={base} BCL={cmd}")'
fi
case ":${LD_LIBRARY_PATH:-}:" in
    *:/opt/bcl*) emit FAIL ld_library_path "contains /opt/bcl: it shadows the envs' libstdc++/libgomp" ;;
    *) emit PASS ld_library_path "no BCL directory on it" ;;
esac

# --- MolGpKa + manifest ---------------------------------------------------------------
MOLGPKA_DIR=${MOLGPKA_DIR:-/opt/MolGpKa}
# the sha256s pka_triage/molgpka.py checks (KNOWN_WEIGHTS)
declare -A WEIGHTS=(
    [weight_acid.pth]=e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39
    [weight_base.pth]=212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d
)
bad=""
for w in "${!WEIGHTS[@]}"; do
    got=$(sha256sum "$MOLGPKA_DIR/models/$w" 2>/dev/null | cut -d' ' -f1)
    [[ "$got" == "${WEIGHTS[$w]}" ]] || bad+=" $w(${got:-missing})"
done
[[ -f "$MOLGPKA_DIR/src/utils/smarts_pattern.tsv" ]] || bad+=" src/utils/smarts_pattern.tsv(missing)"
if [[ -z "$bad" ]]; then
    emit PASS molgpka "$MOLGPKA_DIR: weights' sha256 as molgpka.py expects, SMARTS present"
else
    emit FAIL molgpka "$MOLGPKA_DIR:$bad"
fi
check manifest inenv llb_pka python -c '
import json
m = json.load(open("/opt/llb/manifest.json"))
v, d, b = m["toolchain_version"], m["build_date"], m["bcl"]["installed"]
print(f"toolchain {v}, built {d}, BCL installed: {b}")'

echo
if [[ $nfail -gt 0 ]]; then
    echo "smoke: $nfail check(s) FAILED"
    exit 1
fi
echo "smoke: all checks passed"
