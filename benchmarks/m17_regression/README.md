# M17 production regression — archived evidence

Text evidence from SLURM job 2833238 (`qcpu_18i`, 16 CPUs, 2026-09-29). The gzipped
FASTQ output written by these runs is **not** archived here; only summaries,
commands, parity reports and binary hashes are. Each run directory contains:

- `command.txt` — exact invocation, binary path preserved verbatim
- `binary.sha256` — hash of the M17 release binary that produced the result
- `node.txt` — host the run executed on
- `PASS` — the run script's own pass marker
- `python_parity.txt` — per-sample comparison against the lab Python reference
  (mismatch=0 runs only)
- `baseline_comparison.txt` — diff vs the frozen baseline (mismatch=0) or vs the
  previous mismatch=1 summary (delta=0)
- `i464.summary.tsv` / `i395.summary.tsv` — the SeqMux summary itself

The frozen comparison baselines the run script diffs against are archived under
`baselines/`, because two of the three lived only in `tmp/` and would have been
lost when `tmp/` was cleaned:

- `baselines/v03_i464.summary.tsv`, `baselines/v03_i395.summary.tsv` — the
  mismatch=0 production baselines
- `baselines/m13_i464_mm1.summary.tsv` — the pre-M17 mismatch=1 summary that the
  delta=0 run is checked against

## Result 1 — exact (mismatch=0) regression: PASS

`i464_exact/python_parity.txt`:

```
total=140771720 assigned=87477095 unassigned=53294625
orientation canonical=45393500 swapped=42083595
python_samples=35 seqmux_samples=35
OK: per-sample counts match
```

`i395_exact/python_parity.txt`:

```
total=60523000 assigned=39455648 unassigned=21067352
orientation canonical=18606371 swapped=20849277
python_samples=18 seqmux_samples=18
OK: per-sample counts match
```

Both reproduce the frozen v0.3 production numbers in `docs/REAL_DATA.md` exactly,
and both match the lab Python exact demux per sample. This is the gate the M17
plan required.

## Result 2 — the mismatch=1 delta sweep is inert on this data

`i464_mm1_delta{0,1,2}/i464.summary.tsv` are **byte-identical**:

| metric | delta=0 | delta=1 | delta=2 |
|--------|--------:|--------:|--------:|
| total_reads | 140,771,720 | 140,771,720 | 140,771,720 |
| assigned | 90,802,891 | 90,802,891 | 90,802,891 |
| unassigned | 49,968,829 | 49,968,829 | 49,968,829 |
| ambiguous | 0 | 0 | 0 |
| orientation_canonical | 47,441,842 | 47,441,842 | 47,441,842 |
| orientation_swapped | 43,361,049 | 43,361,049 | 43,361,049 |
| match_rate | 0.645036 | 0.645036 | 0.645036 |

This is not a bug, a no-op flag, or a config-plumbing failure. The parameter is
**provably incapable of binding** on this barcode table, which is why all three
deltas agree.

### Why it cannot bind

A read must be within threshold of two *different* samples for the margin test to
apply at all. With `--mismatches-1 1 --mismatches-2 1` and dual 8 bp barcodes,
samples A and B can both accept some read only if the read sits within 1 bp of
both of A's barcodes **and** within 1 bp of both of B's barcodes. By the triangle
inequality that requires

```
hamming(Barcode1_A, Barcode1_B) <= 2   AND   hamming(Barcode2_A, Barcode2_B) <= 2
```

Checking every sample pair in both production tables:

| Table | samples | pairs satisfying both conditions | min B1 distance | min B2 distance |
|-------|--------:|--------------------------------:|----------------:|----------------:|
| I464 (`tests/fixtures/I464-469erdai_barcode_and_name.csv`) | 35 | **0** | 2 | 3 |
| I395 (`I395erdai_gst_barcode_and_name.csv`) | 18 | **0** | 3 | 4 |

No pair satisfies both. The near misses are instructive: on I464 the closest pair
by Barcode1 (`I464469-B2` / `I464469-E6`) differs by 2 bp in Barcode1 but 6 bp in
Barcode2, so it is nowhere near co-eligible. On I395 the closest pair
(`I395-16-5` / `I395-16-6`) sits at 3 bp on Barcode1 and 4 bp on Barcode2.

So at `mismatches=1` every read has at most one eligible sample, `second_best` is
always `None`, and `finish()` short-circuits to `MatchResult::Match` before the
margin comparison is ever consulted. `ambiguous` stays 0 at every delta.

### What this means

- The M17 implementation is **not** validated as a behaviour change by these runs.
  It is validated as *backward compatible*: `Assigned(delta=0) == Assigned(delta=1)
  == Assigned(delta=2)`, which is the required subset property, trivially satisfied
  here.
- `--min-mismatch-delta` is a **latent** safety control for these datasets. It
  would start to matter only if a future barcode table placed two samples within
  1 bp on both barcodes, or if mismatches were raised above 1.
- A synthetic probe confirms the flag is live on tables where it can bind. With
  single barcodes `A=AAAA`, `B=AACA`, `C=CCCC` and read `AAAT` (d(A)=1, d(B)=2,
  d(C)=4, so margin 1):

  | `--min-mismatch-delta` | result |
  |-----------------------:|--------|
  | 0 | assigned → A |
  | 1 (margin 1 ≥ 1) | assigned → A |
  | 2 (margin 1 < 2) | ambiguous |
  | 3 | ambiguous |

  The `delta=2` threshold produces exactly the intended ambiguity, so the
  behaviour is correct where the margin can actually be evaluated.

Reporting the three identical summaries without this explanation would read as if
the sweep had failed to take effect. The finding — that the lab's barcode tables
are structurally safe from 1-mismatch cross-sample confusion — is itself the
useful result the M17 plan asked for.
