#!/usr/bin/env python3
"""
msani3d.py -- MolSanitizer's 3D + db2 step for a list of named states, in ONE process.

    python msani3d.py --todo todo.smi --chunk 20 -- <msani flags: -3d -f db2 --numconfs ...>

Run by build_shard.py (pass 2) in its work directory, with the env's python. It does
what `msani -i todo.smi <flags>` does for already-prepared states (msani/cli.py with no
2D flags: canonical SMILES, then conformers.gen_conf_chunk), with two differences:

* msani embeds every molecule in a fresh child process (conformers._run_embedding_worker,
  forkserver start method), and that child imports msani.conformers.conformers, which
  parses the torsion and ring libraries at import. Here that module is PRELOADED into the
  forkserver, so it is imported and parsed once per task, not once per molecule. With the
  env on Turbo (NFS) and 140 tasks at once, the per-molecule import is the suspected cause
  of the 8x slowdown on 2026-10-26.
* gen_conf_chunk is called on chunks of --chunk names, and "CHUNK a b" is printed after
  each, so build_shard.py can check and keep those db2s while the process keeps running.

Outputs are MolSanitizer's own, in the current directory: db2/<name>.db2,
msani_error.err, msani_timing.csv, todo.log. MolSanitizer itself is not modified.
"""
import argparse
import multiprocessing as mp
import sys
import time


def main() -> int:
    argv = sys.argv[1:]
    if "--" not in argv:
        raise SystemExit("usage: msani3d.py --todo FILE --chunk N -- <msani flags>")
    cut = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--todo", required=True)
    ap.add_argument("--chunk", type=int, default=20)
    ap.add_argument("--no-preload", action="store_true", help="msani's own behaviour, for A/B timing")
    own, msani_flags = ap.parse_args(argv[:cut]), argv[cut + 1:]

    if not own.no_preload:
        # must precede the first Process: the forkserver starts lazily on first use
        mp.set_forkserver_preload(["msani.conformers.conformers"])

    from pandas import DataFrame
    from rdkit import Chem, RDLogger
    from msani.conformers import conformers
    from msani.io import loggers, parsers
    RDLogger.DisableLog("rdApp.*")

    args = parsers.Sanitycheck(parsers.parseArguments(["-i", own.todo, *msani_flags]))
    loggers.setup_logger("todo.log")
    rows = [line.split()[:2] for line in open(own.todo) if line.strip()]
    t0 = time.perf_counter()
    for i in range(0, len(rows), own.chunk):
        part = rows[i:i + own.chunk]
        df = DataFrame({"smiles": [s for s, _ in part], "ids": [n for _, n in part]})
        df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)
        df["smiles"] = df["mol"].apply(Chem.MolToSmiles)      # as msani's CLI does with no 2D flags
        conformers.gen_conf_chunk(df, args, input_file=f"chunk{i // own.chunk:04d}")
        print(f"CHUNK {i} {i + len(part)} {time.perf_counter() - t0:.1f}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
