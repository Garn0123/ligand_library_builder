# docking/

Tools on either side of a docking run that take finished db2 files, from
MolSanitizer or any other builder. Chunking db2s into balanced work units is
`mol_compiler/` (`llb chunk`); the docking run itself is `rna_pipeline`
`pipeline/library_docking/`.

| `llb` | script | env |
| --- | --- | --- |
| `llb mol2` | `db2_to_mol2.py`: one mol2 per molecule, read back out of its db2 | DOCK |
| `llb mol2-library` | `mol2_library.py`: DOCK-header multi-mol2 shards | PREP |
| `llb collapse-invertomers` | `collapse_invertomers.py`: best score per N-H+ invertomer group | PREP |

These moved here from `run_qupkake/` on 2026-10-26 (tag `llb-full-2026-10` has
them at the old path).

## Mol2 from the db2

```bash
llb mol2 db2_dir/ -o library.mol2 [--library library.tsv]
llb mol2-library library.mol2 --library library.tsv -o mol2_lib --shard-size 1000
```

No second 3D tool such as CORINA. `db2tool tomol2` (dock6_claude `294d235`)
reads a conformer back out of the db2 with DOCK's own db2 reader and mol2
writer, so the mol2 is consistent with the db2 by construction:

| | mol2 | db2 |
| --- | --- | --- |
| coordinates | one conformer set from the db2, lowest internal energy | all sets |
| partial charges | the db2's charges | same |
| atom types | the db2's SYBYL types | same |
| title | the db2's molecule name | same |

Over plain `tomol2`, `db2_to_mol2.py` adds three things:

- **One conformer per molecule.** A db2 can hold one record per rigid
  fragment, all under the same name, and `tomol2 --best` would emit each.
- **Nothing silently missing.** A molecule whose every set is flagged broken is
  taken from its lowest broken set and flagged in the `.tsv`
  (`--no-broken-fallback` to drop them instead).
- **The charge check.** With `--library`, each mol2's summed charge must
  match the library's net charge, and every library name must come out. It
  exits 1 otherwise.

It accepts any db2 files, including mol_compiler chunks: run it per chunk to get
mol2 chunks that line up with the db2 chunks.

**Broken sets.** In a db2_converter db2 (63d6656, mol2db2/clash.py) a set is
flagged broken by one rule only: a rotatable hydrogen within 1.70 Å of another
hydrogen more than two bonds away. DOCK6 samples broken sets unless
`skip_broken yes` (default `no`, conf_gen_hdb.cpp); `llb env DOCK db2tool audit`
counts all_broken / some_broken per molecule.

`mol2_library.py` writes DOCK's header block (name, molecular weight, formal
charge, H-bond acceptors/donors and heavy atoms computed DOCK's way) plus `LLB_`
lines taken from `--library`, the `library.tsv` that `llb names` writes. It has the
columns LLB's own construction stages wrote before 2026-10-26, plus `ms_name`
(MolSanitizer's own name for the row), so old and new libraries read the same.

## Invertomer collapse, after docking

```bash
llb collapse-invertomers scores.parquet --library library.tsv -o scores_collapsed.parquet
```

A protonated amine's N-H+ configuration is fixed in a db2, so both invertomers
are named and docked: MolSanitizer enumerates them as stereoisomers, and its
db2 pass builds each named SMILES as its own db2. This keeps the best score per
group, per target, and writes which member won (`<out stem>.invertomers.tsv`). It needs the
`invertomer_of` column of `library.tsv` (`llb names`). DRAP's
`stage4/collapse_invertomers.py` is the streaming port used at matrix scale.
