# Concurrent Polars GWAMA export — v0.6.2

Validated on 2026-09-30 using Python 3.10.21, Polars 0.20.31, pandas 2.0.3,
NumPy 1.23.5 and pigz 2.8 on local macOS. This release changes post-GWAMA export;
LDSC, PCA and the original GWAMA implementation were not modified.

## Validation

- **19 export tests passed** in `test_postprocess.py`, including native gzip/pigz
  writing, precise numeric strings, tiny p-values, leading-zero identifiers,
  literal NA, empty cells, quoting, tabs/newlines and extra columns.
- Sorting tests cover chromosome aliases, shuffled inputs, equal-position stable
  order, integer keys above 2^53, and the already-sorted fast path.
- Override tests verify that only the selected summary receives N_eff/INFO changes;
  absent columns require the corresponding explicit override.
- A synchronized concurrent-reader test proves two readers overlap, never exceed
  the requested two-reader budget, and retain input order when completion order
  is deliberately reversed. CSV parsing is limited to one thread per file.
- Fresh-process checks verify that both --n_cores forms configure the shared Polars
  pool before import and preserve an explicit POLARS_MAX_THREADS override.
- Invalid/missing rows, duplicate SNPIDs, malformed headers, corrupted gzip data,
  stale/missing/failed status records and archive collisions fail explicitly.
- Compressor, disk-write and publication-race failure tests confirm that new
  partial outputs are not published, temporary files are removed, and existing
  destinations are preserved. A 100,003-row output checks writer batch boundaries.
- Both GPCA wrappers forward the same --n_cores and --gzip_level settings.
- **15 Python interface tests passed**; the separate opt-in native-R interface
  test was skipped in that suite because LDSC_GPCA_TEST_RSCRIPT was not set.
- **10 thread/worker tests passed**, including actual R argument-parser checks
  using /usr/local/bin/Rscript and child-process thread-environment inheritance.

Expected valid outputs retained values, columns and stable ordering. Expected
invalid cases failed before publishing exports. Observed results matched those
expectations. An initial additional filename test exposed Polars 0.20 routing
bracket-containing paths to a compressed-CSV-incompatible scan path even with
`glob=False`; such paths now use compressed bytes as the literal input. The final
filename test passes. Malformed empty/duplicate source headers now fail rather
than being silently renamed by a parser.

## Synthetic performance comparison

All versions used the same deterministic **1,000,000 variants in 22 gzip files**,
180-character Direction strings, precise numerical strings and shuffled positions.
Overrides were N_eff=33000.0 and INFO=0.9. Timings include reads, validation, counts,
sorting, summary creation and both compressed output writes; fixture generation
and upstream GWAMA are excluded.

| Implementation | Read/combine/prepare | Write outputs | Total | Compressed bytes, both tables |
| --- | ---: | ---: | ---: | ---: |
| v0.6.0 pandas, serial gzip level 9 | 13.19 s | 344.06 s | 357.26 s | 143,382,233 |
| v0.6.1 pandas, pigz level 1, 4 workers | 8.05 s | 7.46 s | 15.51 s | 172,069,464 |
| v0.6.2 Polars, 1 reader/pool/compression worker | 4.34 s | 4.52 s | 8.86 s | 172,136,449 |
| v0.6.2 Polars, 4 readers/pool/compression workers | 1.21 s | 1.21 s | 2.42 s | 172,069,464 |

The four-reader Polars run was approximately **6.4 times faster than v0.6.1** on
this fixture. Concurrent reads plus parallel table processing/compression were
approximately 3.7 times faster than the single-worker Polars configuration.
The shared table pool and compression workers can run simultaneously during writes;
these labels are stage worker settings, not hard total-CPU limits.

Whole-process peak RSS was approximately 1.73 GB for the v0.6.1 pandas run and
1.10 GB for the four-reader Polars run (macOS resource.getrusage, bytes). These
are individual local measurements, not a repeated controlled benchmark. Brief
independent regression checks overlapped the first concurrent benchmark. Cache,
storage and machine differences affect timing. The user's 10-million-variant
analysis on the 100-CPU/900-GB server was not run; no server speedup is promised.

**Both complete decompressed outputs were byte-for-byte identical** across all
four implementations/configurations, verified with streamed SHA-256:

- Combined: `bcbf05bf52253d81cce506c36241cdad27fd00f0a666eda2dc25472bf3827ab5`
- Summary: `14b1456c7aff63a70bd9a9468d5e554231fb9dacab703cb2bb926681539ad2cc`

The v0.6.1 single-worker Python gzip fallback was separately benchmarked at 18.58 s
and also produced those same decompressed hashes. The current Polars fallback
writer passed the full-value fixture tests; that earlier fallback timing is not
a performance measurement of the Polars fallback.

## Reproducibility and limits

The benchmark generator, old exporter snapshot, JSON metrics, output hashes and
wheel are retained with the local v0.6.2 artifacts. The wheel was built and
installed to an isolated target for export tests, version checks and both GPCA
help routes. The user's active environment was not reinstalled, and a fresh full
Conda/Docker installation was not performed.

Source variant columns stay strings. Integer sort keys use Polars; decimal or
scientific position spellings use the previous pandas conversion rules for those
keys only. The small run-status CSV still uses pandas. Header names, output file
names and override scope are unchanged. Complete tables and concurrent reader
buffers still reside in memory. An already initialized Polars pool in an API
caller cannot be resized; the actual pool size is reported in the audit.

The original GWAMA SHA-256 remains
`99c041b1e31366e62f1650088e750615d2284047163a55edceeafdf00fcb34d2`.

APIs were verified using the installed version and official documentation:
[Polars read_csv](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.read_csv.html),
[Polars thread_pool_size](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.thread_pool_size.html),
[Python gzip](https://docs.python.org/3.10/library/gzip.html), and
[pigz](https://zlib.net/pigz/pigz.pdf).
