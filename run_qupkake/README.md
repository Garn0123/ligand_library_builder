# QupKake → protomer enumeration → db2

> **For a ZINC library, follow `../UPSTREAM.md`.** It runs this script per shard
> with throwaway names, then `assign_names.py` assigns the 16-character contract
> names (`NAMING_CONTRACT.md`) over the whole library. The `short` names described
> in §3c below are 15 characters and predate the contract; they fail
> `validate_names.py`.

Takes a 2-column `.smi` file, predicts micro-pKa with QupKake, enumerates the
protonation states populated at pH 6.4 / 7.4 / 8.4, dedupes them, labels them,
and writes a `.smi` ready for a DOCK db2 build — plus a statistics report.

## 1. Install

QupKake pulls in torch + pytorch-geometric + pytorch-lightning and pins
`xtb == 6.4.1`. Keep it in its own environment; do not try to share one with
your DOCK/db2 build tools.

```bash
git clone https://github.com/Shualdon/QupKake.git
cd QupKake
conda env create -f environment.yaml
conda activate qupkake
pip install .
```

The README warns that the conda build of xtb is buggy and should be built from
[source (v6.4.1)](https://github.com/grimme-lab/xtb/releases/tag/v6.4.1):

```bash
export XTBPATH=/path/to/xtb
```

Linux binaries ship with the package and are used if neither the conda package
nor `$XTBPATH` is set, but the from-source route is what the authors recommend.

Then, in the same environment:

```bash
pip install rdkit pandas          # already present via QupKake, but pin them
python qupkake_protomers.py --help
```

Smoke-test the enumeration logic without waiting on xtb:

```bash
python selftest_make_sdf.py                       # writes fake_qupkake_output.sdf
python qupkake_protomers.py test_input.smi -o out --from-sdf fake_qupkake_output.sdf
```

## 2. Run

```bash
python qupkake_protomers.py ligands.smi \
    -o protomers \
    --ph 6.4 7.4 8.4 \
    --margin 1.0 \
    --min-population 0.01 \
    --max-states 8 \
    -mp 8
```

Input is `SMILES name`, whitespace-separated, no header — the same format the
DOCK 3D pipeline expects.

Outputs in `protomers/`:

| file | contents |
|---|---|
| `protomers.smi` | 2-column `SMILES name` — **this is the db2 input** |
| `protomers.csv` | per-protomer: parent, charge, +/- atom counts, zwitterion flag, pH list, population, n changes |
| `name_map.tsv` | protomer → parent lookup for joining OUTDOCK results back |
| `stats.json` | the full statistics block |
| `rejected_input.smi` | SMILES RDKit could not parse |
| `qupkake_work/` | QupKake's own `raw/`, `processed/`, `output/`, `logs/` |

### Handing off to db2_converter

`protomers.smi` is already in db2_converter's input format — `SMILES LIG_NAME`,
optionally with a third `CONF_NUM` column (`--max-conf 600`) that overrides `-n`
per molecule.

```bash
conda activate db2_converter
workingpath=$(readlink -f build)
OMP_NUM_THREADS=1 build_ligand \
    -i protomers/protomers.smi \
    -n 600 -nr 30 --keep_max_conf \
    -m ccdc \
    --workingpath $workingpath --outputpath $workingpath \
    --checkstereo --reseth --rotateh --dock38
```

Or chain it from the enumeration step:

```bash
python qupkake_protomers.py ligands.smi -o protomers \
    --db2-cmd 'build_ligand -i {smi} -n 600 -nr 30 --keep_max_conf -m ccdc \
               --workingpath build --outputpath build --checkstereo --reseth --rotateh'
```

**Do not pass `--sampletp`.** It calls
`unicon -i ... -o ... -t single -p single`, which collapses each molecule to a
*single* tautomer and protomer — it would discard the protonation states you
just spent QupKake time computing. That flag is the alternative to this
pipeline, not a companion to it.

`--mergeiso` — **not for contract-named libraries** (see `../UPSTREAM.md` §5: after stereo closure it is a no-op, and if a split slipped through it would hide it). For ad-hoc runs of this script alone, it sets the db2 M-line name to
`name.split(".")[0]`, and `--checkstereo` is what appends the `.0`/`.1`
stereoisomer suffixes — so `--mergeiso` merges enumerated *stereoisomers* under
one db2 name while leaving your `_p1`/`_p2` protomer suffixes intact, which is
usually what you want. (There is a commented-out line in `pipeline.py` that also
splits on `_`; if anyone uncomments it, your protomers will silently merge too.)

Then verify what was actually built:

```bash
python verify_db2_build.py protomers/protomers.csv \
    --outputpath build \
    --faillist build/protomers.smi.ccdc.faillist
```

Re-running enumeration with different pH or margin settings, without paying for
prediction again:

```bash
python qupkake_protomers.py ligands.smi -o wider \
    --from-sdf protomers/qupkake_work/output/qupkake_output.sdf --margin 1.5
```

## 3. How states are chosen

QupKake writes **one SDF record per reaction site**, tagged:

- `idx` — heavy-atom index of the site
- `pka_type` — `acidic` = proton **loss** here, `basic` = proton **gain** here
  (so a `basic` pKa is the pKa of the conjugate acid)
- `pka` — the predicted micro-pKa

Per site, Henderson–Hasselbalch gives the fraction in the changed form:

```
acidic:  f_deprot = 1 / (1 + 10^(pKa − pH))
basic:   f_prot   = 1 / (1 + 10^(pH − pKa))
```

Sites further than `--margin` from the pH are locked to their dominant form.
Sites inside the margin are enumerated both ways, giving 2^n states. Each
state's population is the product of its per-site fractions across **all**
sites; states below `--min-population` are dropped, then the top
`--max-states` survive per molecule per pH.

Identical SMILES arising at more than one pH collapse into one entry whose
`ph_values` column lists every pH where it appeared, which is usually most of
them — for a typical drug-like molecule the three pH points differ only for
sites with pKa in roughly 5.4–9.4.

Naming, via `--name-style`:

- `index` (default) — `LIG_p1`, `LIG_p2` — shortest, safest
- `charge` — `LIG_qp1_2`, `LIG_qm1_3` — charge visible in OUTDOCK
- `ph` — `LIG_pH6p4-7p4_2` — most informative, longest

## 3b. Stereochemistry

Two axes, treated differently. A **protomer** is the same compound in another
state — enumerate it, because you're uncertain. A **stereoisomer** is a
different compound — specify it, because you shouldn't be.

`--stereo` picks the policy:

| mode | behaviour |
|---|---|
| `carry` (default) | keep whatever the input declares; db2_converter's `--checkstereo` enumerates anything undefined into `NAME.0`, `NAME.1` |
| `enumerate` | expand undefined centres here, run QupKake **once per parent**, fan the protomers across the isomers |
| `keep` | trust the SDF's 3D-derived tags (rarely what you want) |
| `strip` | flatten everything |

With `enumerate`, isomers are generated from the *re-parsed canonical SMILES* —
the exact string QupKake receives — so every isomer already shares QupKake's
atom ordering and the protonation change is applied straight to the isomer.
Nothing is transferred, which matters because copying chiral tags onto a rebuilt
skeleton does not carry E/Z double-bond stereo. `onlyUnassigned=True` means a
partially-defined molecule keeps what it declared and only the gaps are filled.

Protomer states are enumerated once per parent and indexed on the first isomer,
so state `1` means the same thing for every isomer of that molecule — the output
is a clean `n_states x n_isomers` grid. QupKake (the xTB-bound step) runs once
per parent either way; a molecule with 4 isomers and 2 protomers costs one
prediction and yields 8 db2 entries.

Before spending prediction time, see what you're committing to:

```python
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
for line in open("ligands.smi"):
    smi, name = line.split()[:2]
    n = len(tuple(EnumerateStereoisomers(Chem.MolFromSmiles(smi))))
    if n > 1:
        print(f"{name}\t{n} isomers\t{smi}")
```

Anything that prints is a molecule where you're about to make a choice.
`--max-stereo` caps the expansion (default 32, matching db2_converter's own
`2**5` limit); `--stereo-embed` discards configurations RDKit can't embed, which
removes impossible ring stereochemistry at some cost.

Note that ZINC SMILES usually arrive already protonated at reference pH. QupKake
wants neutral input, so the pipeline neutralizes on read — don't pass
`--no-neutralize` on ZINC tranches.

## 3c. Names that survive the db2 field

`--name-style short` (default) builds a fixed-width name:

```
C 000012345678 1 0
| |            | |
| |            | stereoisomer  (base62: 0-9 A-Z a-z)
| |            protomer state  (base62)
| the ZINC digits, kept verbatim
leading ZINC rewritten to C
```

15 characters, inside db2's 16-char `name[-16:]` field, so nothing truncates.
Base62 gives 62 slots per field — more than db2_converter's 32-isomer cap.

**There is deliberately no dot.** `--mergeiso` sets the db2 name to
`name.split(".")[0]`, so `C000012345678.10` would collapse every protomer *and*
every stereoisomer of that molecule into one db2 name. With no dot, `--mergeiso`
is a harmless no-op and you can leave it on or off. This is also why the
`ZINC000012345678.0` scheme mangles: at 18 characters the tail wins and you get
`NC000012345678.0`.

If a name won't fit — a long internal ID, a non-ZINC identifier — the whole set
switches to compact serials (`C00`, `C01`, ...) rather than mixing schemes, and
says so. Either way `name_map.tsv` carries the full parent name, protomer index,
stereo index, charge, and pH list. Keep it: db2_converter hardcodes
`longname="NO_LONG_NAME"`, so the original name is *not* recoverable from the
db2 itself.

## 4. Things that will bite you

**db2 truncates names to the last 16 characters.** The bundled `mol2db2`
(both `mol2db2/hierarchy.py` and `mol2db2_py3_strain/hierarchy.py`) writes the M
line as `'M %16s ...' % mol2data.name[-16:]`. Appending `_p1` to an already-long
name chops the front off — `benzene_no_sites_p1` is stored as
`zene_no_sites_p1` — and two different parents can collapse to the same db2
name. Both scripts audit this: `stats.json` reports `db2_naming` before you
build, and `verify_db2_build.py` reports post-truncation collisions in what was
actually written. Starting from 16-character ZINC IDs, use short internal IDs
for the build and join back through `name_map.tsv`.

**Ligand names must be `[A-Za-z0-9_-]`.** `db2_converter.py` uses the name as a
directory, as a filename, and — unquoted — inside
`subprocess.run(f"rm -r {zinc} {zinc}.smi", shell=True)`. A name containing a
space, a semicolon, or a glob character is at best a broken build. Parent names
are sanitized on read and any substitution is reported; dots are replaced too,
since `--mergeiso` splits the name on the first dot.

**`chemistrycheck` cannot police protonation.** It compares each conformer to
the input by standard InChI, which normalizes mobile hydrogens and zwitterions —
glycine's zwitterion and neutral form share an InChI. So a molecule that comes
back from the sampler/UNICON round trip in a different protonation state can
still pass. That is precisely the check you would want in a protomer workflow,
and it is the one that does not work. `verify_db2_build.py` compares the AMSOL
total charge written into the db2 (second M line, `solvdata.totalCharge`)
against the formal charge you requested, which does catch it.

**Independent-site approximation.** Micro-pKa values are quoted for the
reference microstate. Protonating one site shifts every other site's pKa, and
the script does not model that coupling. For a monoprotic or well-separated
polyprotic molecule this is fine; for something like a polyamine with adjacent
basic nitrogens the populations will be wrong. Widen `--margin` there rather
than trusting the ranking.

**QupKake was not trained on zwitterions** and returns values near 7 when it is
handed one, so the script neutralizes and takes the largest fragment before
prediction (`--no-neutralize` to disable). Zwitterionic *outputs* are fine and
counted separately in the stats — it's zwitterionic *input* that degrades the
prediction.

**Carbon acids are absent from the training set** — the paper states the model
cannot even identify the reaction site for aliphatic C–H deprotonation. Also
expect weaker accuracy on aliphatic alcohols, thiols, phenols, and thiophenols
than on amines, carboxylic acids, anilines, and N-heteroaromatics.

**The 3D-embedding stereo trap.** QupKake embeds every molecule in 3D before
predicting, so its output SDF carries chiral tags read off an arbitrary
conformer. Parse it naively and `CC(C)NCC(O)COc1cccc2ccccc12` comes back as
`CC(C)[NH2+]C[C@@H](O)COc1cccc2ccccc12` — an unrequested enantiomer commitment
across your whole library. `--stereo input` (the default) copies the input
SMILES' stereo back over the protomer; `--stereo strip` flattens; `--stereo
keep` accepts the SDF's version.

**`-t/--tautomerize` runs an xtb tautomer search** and replaces the molecule
with the lowest-energy tautomer, so your protomers may be built on a different
skeleton than you submitted. It's slow, and it breaks the index correspondence
the stereo fix relies on (the script detects this and falls back to stripping
stereo). Worth it when your input tautomers are unreliable; skip it otherwise.

**QupKake crashes are excluded, not passed through.** QupKake drops molecules and sites it fails on; the script audits its working tree and lists them in `qupkake_failed.tsv` (with `--from-sdf` on a bare SDF, the audit cannot run and says so).

**Molecules with no predicted sites** pass through unchanged as `NAME_p1` with
`note=no_sites_predicted`, so nothing silently disappears between the input and
the db2 build. `input.molecules_without_sites` in the stats should match your
expectation for how many of your ligands are genuinely non-ionizable.

## 5. Statistics reported

`stats.json` and the terminal summary cover: molecules read / predicted /
rejected; total sites split acidic vs basic with pKa min/median/max and sites
per molecule; unique protomers overall and per pH, how many are present at all
three pH values vs unique to one, and the expansion factor; charge breakdown
(cationic / anionic / neutral / zwitterionic plus the full net-charge
histogram); protomers per molecule (min/mean/median/max, histogram, count of
single-protomer molecules); heavy-atom size range; and the db2 naming audit.
