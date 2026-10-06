# ZINC22 tranche fetching — a user's manual

Getting a few thousand `.db2.tgz` tranches off `files.docking.org` onto a shared
cluster filesystem, verifying them, and reading them into a pipeline — without
being rude to a single academic web server or drowning GPFS in small files.

Written after a full run against `zinc-22a` / `zinc-22d` / `zinc-22u` on Great
Lakes. Everything in the Gotchas section actually happened.

---

## Quickstart

```bash
chmod +x *.sh                          # yes, really — see Gotcha 1

# 1. get a URL list out of the tranche browser's download script (curl, wget or PowerShell)
./extract_urls.sh download_all.sh > urls.txt

# 2. first pass
./run.sh urls.txt 4

# 3. integrity check, purge anything truncated
./verify.sh --purge

# 4. what's left?
./status.sh urls.txt

# 5. retry the gaps, gently
./run.sh retry.txt 2

# 6. repeat 3–5 until retry.txt is empty
# 7. confirm the leftovers are real absences, not a bad URL list
./tranche_pattern.sh permanent.tsv
```

## SMILES (2D) instead of db2

The same loop fetches 2D tranches (`zinc22/2d/H17/H17M100.smi.gz`): pick SMILES
in the CartBlanche22 tranche browser, build `urls.txt` as above, and run.
`status.sh` and `verify.sh` match `.smi.gz` as well as `.tgz`. `gzip -t` validates
both the same way. 2D filenames have no charge letter; `tranche_pattern.sh`
reports them as `none (2D file)`.

**Check the volume first.** 2D bins grow fast: H17 is 1.3 GB, H22 21 GB, H28
304 GB. To draw N molecules per bin, use `sample_2d.py`: it streams the front of
each tranche file and stops, a few MB per bin. `../UPSTREAM.md` covers the
whole path from here to db2.

---

## What's in here

| File | Purpose |
|---|---|
| `extract_urls.sh` | URLs out of a CartBlanche22 curl/wget/PowerShell script, cleaned for `run.sh`. |
| `fetch.sh` | Fetches one URL. Idempotent, purges partials, retries 5xx, bails on 404. |
| `run.sh` | Drives `fetch.sh` through GNU parallel at a polite concurrency. |
| `status.sh` | Diffs expected-vs-present on disk. **The authoritative completion check.** |
| `verify.sh` | `gzip -t` every archive; `--purge` deletes corrupt ones. |
| `parse_log.sh` | Classifies HTTP errors in `wget.log`. Diagnostic only. |
| `tranche_pattern.sh` | Sanity-checks whether 404s are chemically real or a bug. |
| `db2_sources.py` | Streams `.db2` records straight out of the tarballs. Don't extract. |
| `sample_2d.py` | N SMILES per heavy-atom bin from ZINC22 2D, streamed, without downloading the bins. |

Generated as you go: `urls.txt` (yours), `missing.txt`, `retry.txt`,
`permanent.tsv`, `failed.tsv`, `wget.log`, `joblog-*.tsv`, `corrupt.txt`.

`permanent.tsv` accumulates across passes on purpose — once a tranche 404s it
stays 404. Everything else is rewritten per pass.

---

## Doing this for a new set of tranches

**1. Get the URL list.** The tranche browser at `cartblanche22.docking.org`
generates a shell script of `wget` commands. Strip it down to bare URLs:

```bash
./extract_urls.sh download_all.sh > urls.txt     # curl, wget or PowerShell script
wc -l < urls.txt
```

Not `grep -o 'https://[^ ]*'`: it keeps a closing quote on quoted URLs and a
trailing `\r` on scripts with Windows line endings. With a `\r`, every URL
differs from the path on disk and `status.sh` reports every file missing
forever. `extract_urls.sh` strips both and rejects anything that isn't a file
on files.docking.org.

Don't run the generated script as-is. It ships `-r -l7 -np -A '<pattern>'`,
which are recursive-crawl flags copied from directory-download examples. Your
URLs point at individual files, so recursion buys nothing — but if the server
ever returns HTML instead of a tarball (error page, redirect), `-r` sends wget
crawling seven levels deep. Multiply by a few thousand and you have a problem.
`fetch.sh` uses `-x -nH` instead, which preserves the same directory tree with
no recursion.

**2. Pick a landing directory** on scratch, not home. Run everything from there;
`-x -nH` writes relative to the working directory, mirroring the URL path:

```
zinc22/zinc-22a/H05/H05M000/a/H05M000-O-aaaaaa.db2.tgz
```

**3. Run the loop** from Quickstart. Expect two or three passes. A run of ~2400
files at `-j 4` takes roughly 15–20 minutes plus retries.

**4. Stop when `status.sh` reports `-> retry.txt  0`.** Not when the log looks
clean — see Gotcha 4.

---

## Gotchas

### 1. `Exitval 126` on every job — script isn't executable

The joblog shows exit 126 across the board and every job finishes in
milliseconds. 126 is "found it, can't run it"; 127 is "command not found."
Nothing hit the network.

```bash
chmod +x *.sh
```

Then start a **fresh joblog**. A log containing both the failed and re-run
passes has duplicate sequence numbers, and `--resume-failed` matches jobs by
sequence number against the input file. Feeding it an ambiguous log will run the
wrong command for the wrong URL.

### 2. The `-c` append trap — the one that silently corrupts data

`wget -c` resumes by **appending to whatever is already on disk**. If a failed
attempt left a truncated file, or wrote an HTML error page under a `.tgz` name,
`-c` appends to that garbage. The result is a file that looks complete, fails
`gzip -t` forever, and never converges no matter how many times you retry.

`fetch.sh` handles this by `rm -f`-ing the target before every attempt. If you
modify the script, keep that. If you ever retry by hand:

```bash
while read -r url; do rm -f "${url#https://files.docking.org/}"; done < retry.txt
```

### 3. wget 1.19 doesn't retry HTTP 5xx

`--tries` only covers network-level failures. HTTP 500 responses fail on first
contact and are never retried. `--retry-on-http-error=500,502,503,504` fixes
this but needs **wget 1.20+**.

Check with `wget --version | head -1`. On 1.19 (which is what the Great Lakes
base environment ships), `fetch.sh`'s bash retry loop supplies the behaviour:
escalating backoff on 5xx and TLS drops, immediate bail on 404/403/410 so
permanent failures don't burn attempts.

If you'd rather have it natively: `conda install -c conda-forge wget`, then you
can simplify `fetch.sh` back to a single call with `--retry-on-http-error`.

### 4. The log undercounts failures — trust the filesystem

TLS-level failures appear in `wget.log` as bare lines with **no URL attached**:

```
GnuTLS: Error in the pull function.
Unable to establish SSL connection.
Read error (Success.) in headers.
```

There is no way to recover which file each one belonged to. In one real run the
log named 37 failures (24× 404, 13× 500) while the filesystem showed **52**
files missing. Building a retry list from the log would have silently dropped 39
files that exist and were wanted.

**Always use `status.sh` to decide what to retry.** `parse_log.sh` is for
understanding *why* things failed, never for deciding *what* to refetch.

### 5. Exit code 8 means "HTTP error," not "404"

wget collapses 404, 403, and every 5xx into exit 8. Bucketing exit-8 as
"permanently gone" throws away files that would succeed on retry. Evidence from
a real log:

```
10:29:28  ERROR 500   H10P050-N-daaaaa
10:29:30  SUCCESS     H10P050-M-daaaaa   ← same tranche, 2 seconds later
```

`fetch.sh` reads the actual status code out of wget's stderr and splits
`permanent.tsv` (404/403/410) from `failed.tsv` (everything else) so this
distinction survives.

### 6. Successful-looking files can be truncated

A dropped connection can leave a short file that wget reports as exit 0. `[[ -s
$f ]]` won't catch it. Only `gzip -t` will. Run `./verify.sh --purge` after
every pass and again before handing anything to the pipeline. Finding
corruption here costs seconds; finding it three hours into a DOCK run doesn't.

### 7. Bursty TLS errors mean back off, not push harder

Four `Read error in headers` followed by five consecutive `GnuTLS` failures in
one second is the server saying it's at capacity. Retrying at higher concurrency
has a *lower* success rate than retrying at `-j 2`. This isn't only etiquette —
it's what actually finishes faster.

Also worth checking when failures cluster in time: your own scratch quota.

```bash
df -h . ; lfs quota -h -u "$USER" /scratch 2>/dev/null
```

A quota filling mid-run produces the same signature as server trouble (wget exit
3 is local file I/O).

### 8. Some 404s are real — check the pattern before discarding

The tranche browser emits the full combinatorial grid of heavy-atom count ×
logP × charge whether or not a cell was ever populated. Real absences cluster:
in one run **every** 404 was in an `M` tranche (negative cLogP) at low heavy-atom
count, and not one `P` tranche 404'd. Charged, very polar, very small molecules
are a genuinely sparse corner of make-on-demand chemical space.

```bash
./tranche_pattern.sh permanent.tsv
```

If failures skew hard toward `M` and low HAC, discard with confidence. If
they're spread evenly across logP signs and charges, or they all share one
directory or one charge letter, the URL list is wrong — fix that rather than
writing off data you wanted.

---

## Reading the data

**Do not untar.** Each archive holds thousands of tiny `.db2` members; across a
few thousand archives that's millions of small files. On GPFS you'll exhaust the
inode allocation before the byte quota, metadata operations at that scale are
what shared filesystems are worst at, and you gain nothing — the archive is
already the right unit of parallelism.

`db2_sources.py` streams members in place:

```python
from db2_sources import iter_sources, find_archives, shard

archives = shard(find_archives("zinc22"), args.shard, args.nshards)

for fi, (label, fh) in enumerate(iter_sources(archives, on_error)):
    for src_idx, (rec, complete) in enumerate(iter_records_bytes(fh)):
        ...   # unchanged
```

Smoke-test any new tranche set before a full run:

```bash
python3 db2_sources.py zinc22/zinc-22a/H05/H05M000/a/H05M000-O-aaaaaa.db2.tgz
# -> "N db2 members, 0 unreadable archives"
```

Three things to know:

- **Streaming handles are one-shot.** Opened `r|gz`; each handle is valid only
  until the generator advances. Consume fully, never buffer into a list or hand
  to a thread pool.
- **Manifest paths become `archive::member`.** Anything downstream that reopens
  paths must split on `::` (`split_label()` does it). The upside is exact
  per-molecule provenance.
- **Keep archives whole within a shard.** Per-shard bin-balancing state
  (`writers`, `counts`, `heap`) assumes it. `shard()` respects this.

Decompression is CPU-bound and single-threaded per archive, so shard across
SLURM tasks. If you profile as I/O-bound rather than CPU-bound — likely on GPFS
with thousands of archives — staging a batch to node-local `/tmp` helps more
than adding ranks.

---

## Globus — use it instead, when it works

HTTP is the fallback, not the intended route for this volume. ZINC exposes
**"UCSF BKS Lab Guest Collection"**. As of this writing its certificate is
expired; when that's fixed, prefer it.

```bash
globus login
globus endpoint search "UCSF BKS Lab Guest Collection"
globus ls "$ZINC:/"                    # discover the real path layout first
```

Globus paths do **not** match the HTTPS URLs — the collection exposes the lab's
underlying NFS mounts (`/nfs/exd/zinc-22a/...`, `/nfs/exl/...`, `/nfs/exc/...`).
`globus ls` only lists directories, so to check existence in bulk, list each
parent once and diff against expected filenames. Globus also doesn't follow
symlinks, so a path that works over HTTPS can come back not-found; that's worth
asking the ZINC folks about rather than assuming the data is gone.

Once the paths are known, one batch transfer replaces the whole loop above:

```bash
globus transfer --batch batch.txt --sync-level checksum --preserve-mtime \
    --label "zinc22-H10" "$ZINC" "$MINE"
```

Server-to-server, survives logout, checksums included, and `--sync-level
checksum` makes re-runs skip what's already correct. All of Gotchas 2–7 simply
don't apply.

---

## Being a good citizen

`files.docking.org` is one server at UCSF on an NIH grant. Not a CDN.

- **4 concurrent connections max** on a first pass, 2 on retries. More doesn't
  go faster; per-connection bandwidth is the bottleneck.
- **Don't use a SLURM job array for downloads.** 2400 array tasks across dozens
  of nodes means dozens of simultaneous connections, and it burns scheduler
  slots on I/O-bound work. One job, or a `tmux` session on a data transfer node.
- **Check whether compute nodes have outbound internet** — on many clusters they
  don't. The DTN is where bulk transfers belong anyway.
- **Back off on failure** rather than looping tight. Escalating sleeps, three
  passes, then stop and diagnose.
- **Re-runs should be cheap.** `fetch.sh` skips anything already valid, so a
  full re-run costs local CPU and zero requests.
- **Report breakage.** The expired Globus certificate and persistent 500s are
  both worth a short email — they're the kind of thing that gets fixed within
  days once someone mentions it.
