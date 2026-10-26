"""A stand-in for MolSanitizer 0.7.4, for testing build_shard.py without it.

Only the surface build_shard.py touches: msani.api.Msani(...).run(df),
msani.filtering.filters.Filters.remove_exotic_chem_to_db2(df), msani.__version__,
and the `msani` CLI (bin/msani) for the 3D step. The chemistry is canned (see
api.py); the id suffixes and the db2 header lines follow the real 0.7.4 formats
(api.expand_ids, stereoisomers.py, db2/layout.serialize_db2).
"""
__version__ = "0.7.4-fake"
