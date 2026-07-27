# tests

Automated suite for db2pipe. Pure-Python fixtures — no committed binary data, no
network, no cluster. Every fixture tree is generated on a temp dir by
`db2gen.py`, so runs are deterministic and self-cleaning.

## Run

```bash
pip install -r ../requirements-dev.txt   # pytest only
python3 -m pytest                          # from the mol_compiler dir
python3 -m pytest -v                        # verbose
python3 -m pytest tests/test_parallel.py    # one file
```

`pytest` is a **dev-only** dependency; the pipeline itself remains stdlib-only.

## Layout

| File | Covers |
|------|--------|
| `db2gen.py` | synthetic record + fixture-tree generators (not a test module) |
| `conftest.py` | sys.path setup + the `run_script` subprocess fixture |
| `test_parsers.py` | text vs bytes parser equivalence; buffer boundaries; truncation; the `E`-line divergence (open question) |
| `test_parallel.py` | `make_shards`→`p1`→`p2`→`p3`; manifest positions match assembled bytes; lossless round-trip |
| `test_pipelines.py` | serial `01` and the parallel pipeline agree on the id-multiset; every chunk is whole records |
| `test_labeling.py` | `02`/`03` duplicate suffixing; stale-label rejection aborts without corrupting data |
| `test_units.py` | `db2common` pure functions: id extraction, relabel, weights, dedup, gzip sniff, tsv validation |
| `test_edge_cases.py` | truncated input, plaintext-with-`.gz`-name, hardlink dedup, empty shard |

## Notes

- Stages are driven as subprocesses (`run_script`) to exercise the real CLI and
  exit codes.
- `test_parsers.py::test_interior_e_line_divergence_documents_open_question`
  pins the current behavior of the `\nE\n` assumption (README "Open questions").
  If that assumption is ever hardened, update that test to match.
