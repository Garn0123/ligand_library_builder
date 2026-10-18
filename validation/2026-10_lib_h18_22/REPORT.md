# Library run lib_h18_22

Started 2026-10-16T12:38:41, finished 2026-10-16T20:43:06. LLB e03b65a6e6ec29f456202ab82cf64e2efc79d453.

## Attrition (molecules, by heavy-atom bin)

| heavy_atoms | inputs | parents | named | passed_checks | built |
| --- | --- | --- | --- | --- | --- |
| 18 | 10000 | 10000 | 9998 | 9997 | 9786 |
| 19 | 10000 | 10000 | 10000 | 10000 | 9762 |
| 20 | 10000 | 10000 | 10000 | 10000 | 9790 |
| 21 | 10000 | 10000 | 9999 | 9999 | 9791 |
| 22 | 10000 | 10000 | 10000 | 10000 | 9807 |
| all | 50000 | 50000 | 49997 | 49996 | 48936 |

## Stages

| stage | status | jobs | retries | result |
| --- | --- | --- | --- | --- |
| parents | done |  | 0 | {"lost": 0, "total": 50000, "fraction": 0.0, "threshold": 0.01, "what": "inputs rejected"} |
| triage | done |  | 0 | {"lost": 0, "total": 50000, "fraction": 0.0, "threshold": 0.01, "what": "parents MolGpKa failed on"} |
| qupkake | done | 64160613 64169053 64181918 | 1 | {"lost": 408, "total": 19064, "fraction": 0.0214, "threshold": 0.1, "what": "routed parents QupKake failed on (MolGpKa's states used instead)"} |
| merge | done |  | 0 | {"aromatic_n_swapped_parents": 5645} |
| names | done | 64190979 | 0 | {"lost": 3, "total": 50000, "fraction": 6e-05, "threshold": 0.02, "what": "parents with no name"} |
| checks | done |  | 0 | {"lost": 3, "total": 98091, "fraction": 3e-05, "threshold": 0.005, "what": "names excluded by check-stereo", "split": 2, "crashed": 1} |
| db2 | done | 64194666 64198684 | 1 | {} |
| verify | done |  | 0 | {"db2": {"lost": 1061, "total": 49997, "fraction": 0.02122, "threshold": 0.05, "what": "named parents with nothing built"}, "end_to_end": {"lost": 1064, "total": 50000, "fraction": 0.02128, "threshold": 0.05, "what": "input parents with no db2"}} |

## Settings

```json
{
  "validate": false,
  "shard_size": 10,
  "preflight": 10,
  "qupkake_time": "00:30:00",
  "names_cpus": 16,
  "names_mem": "32G",
  "names_time": "03:00:00",
  "db2_ncpus": 8,
  "db2_names_per_task": 25,
  "db2_time": "01:00:00",
  "db2_mmff": "restrained",
  "control_time": "02:00:00",
  "control_mem": "8G",
  "max_array_tasks": 4500,
  "stop": {
    "parents": 0.01,
    "triage": 0.01,
    "qupkake": 0.1,
    "names": 0.02,
    "checks": 0.005,
    "db2": 0.05,
    "end_to_end": 0.05
  }
}
```

## Versions

```json
{
  "llb_commit": "e03b65a6e6ec29f456202ab82cf64e2efc79d453",
  "prep": {
    "aromatic_n_version": "2",
    "coupling_version": "1",
    "molgpka_weights_sha256": {
      "weight_acid.pth": "e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39",
      "weight_base.pth": "212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d"
    },
    "numpy": "1.26.4",
    "pandas": "2.3.1",
    "python": "3.9.23",
    "qupkake": "0.1+74.gf7a294e",
    "rdkit": "2025.9.2",
    "site_filter_version": "6",
    "torch": "2.1.0"
  },
  "pka": {
    "aromatic_n_version": "2",
    "coupling_version": "1",
    "dimorphite_dl": "2.1.0",
    "molgpka_weights_sha256": {
      "weight_acid.pth": "e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39",
      "weight_base.pth": "212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d"
    },
    "numpy": "1.26.4",
    "pandas": "3.0.6",
    "python": "3.11.17",
    "rdkit": "2024.3.5",
    "site_filter_version": "6",
    "torch": "2.4.1"
  },
  "qupkake": {
    "aromatic_n_version": "2",
    "coupling_version": "1",
    "molgpka_weights_sha256": {
      "weight_acid.pth": "e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39",
      "weight_base.pth": "212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d"
    },
    "numpy": "1.26.4",
    "pandas": "2.3.1",
    "python": "3.9.23",
    "qupkake": "0.1+74.gf7a294e",
    "rdkit": "2025.9.2",
    "site_filter_version": "6",
    "torch": "2.1.0"
  },
  "db2c": {
    "error": "oftware.\n    For more information please use:\n        \n        $ module help Chemistry\n\n\n\n    Provides Chemistry software.\n    For more information please use:\n        \n        $ module help Chemistry"
  }
}
```

## Cost

```

== triage triage: 5000 finished shard(s)
  50,000 parents in 492.7 s: 9.9 ms per parent

== QupKake protomers_qk: 1907 finished shard(s)
  19,064 parents, 147.04 CPU-h: 27.8 s per parent overall; per shard median 26.8 s, p90 36.8 s
  slowest shard shard_00195: 10.8 min for 10 parents (size --time from this, with headroom)

== db2 build db2_run: 98,088 molecule(s) timed, from 3924 shard file(s)
  time per molecule: median 41.8 s, p90 77.8 s, max 497.8 s; total 1,246.29 CPU-h
  built 95,696 of 98,088 (97.6%); failures: 1generate 1341, 3chemistrycheck 1042, 5amsolfail 9
  conformers produced (n_out): 0 0, 1-24 1550, 25-99 5386, 100-299 17149, 300+ 71611  -> 1.6% under 25
  db2.gz size: median 445 KB, mean 469 KB, total 42.82 GB
  per 100k molecules at these rates: 1,271 CPU-h, 45 GB of db2
  by heavy atoms: H18: 36s (n=20457)  H19: 39s (n=19934)  H20: 43s (n=19958)  H21: 46s (n=18832)  H22: 48s (n=18907)

```

## Next: the mol2 library

```bash
llb mol2 db2_run/out --library library/library.tsv -o library/library.mol2
llb mol2-library library/library.mol2 --library library/library.tsv -o mol2_lib --shard-size 1000
```
