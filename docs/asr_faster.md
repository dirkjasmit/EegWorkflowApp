# A faster and deterministic implementation of Artifact Subspace Reconstruction

*EegWorkflowApp technical report, October 2026*

## Summary

The burst-correction step of `clean_rawdata` (ASR) was reimplemented as
`asr_process_fast` and `clean_asr_faster`. The algorithm is unchanged; the
only difference is how the local covariance is obtained. On a 126-channel,
190 s recording the processing step went from 24.3 s to 1.1 s (22×), the ASR
step in the application from 11.4 s to 4.7 s, and the peak intermediate
memory from 12.6 GB to a few megabytes. Given the same update grid the two
implementations agree to numerical precision (correlation 1.000000, RMS
difference 4 × 10⁻¹⁰ µV). Because the new version needs no memory blocks, its
output no longer depends on the memory available, whereas the original
differs from itself by 26 % in RMS (correlation 0.965) between its 64 MB
default and a 1 GB setting.

## Background

ASR calibrates on a clean stretch of the recording, obtaining a mixing matrix
**M** and a threshold vector **T**. It then moves through the data and, at
regular positions, compares the local covariance with that norm: the
eigenvectors of the local covariance are computed, directions carrying more
variance than calibration allows are marked for replacement, and a
reconstruction matrix **R** rebuilds them from the remaining directions.
Between two update positions the data are multiplied by a cosine-weighted
blend of the previous and the current **R**.

Throughout, **C** is the number of channels, **S** the number of samples,
**N** the window length in samples (`windowlen` × sampling rate, 0.5 s by
default) and *step* the update interval (`stepsize`, half a window by
default; 125 samples at 500 Hz).

The expensive quantity is the local covariance. The original forms the outer
product **x(t)·x(t)ᵀ** at every sample, flattens it to a column of length C²,
and smooths it over time with a trailing moving average of length N,

> Xcov(:, k) = (1/N) · Σ over the N samples ending at k of x(t)·x(t)ᵀ ,

which is an array of C² rows and S columns: 12.6 GB for the recording used
here. Because it does not fit in memory, the data are processed in blocks
sized by a `maxmem` argument (64 MB by default, about 5900 blocks for this
file), with the filter state, the covariance state and the last **R** carried
from block to block. Of the smoothed covariance, only every *step*-th column
is used, so 124 of every 125 computed columns are discarded.

## Implementation

The new version computes the same covariance only where it is used. At an
update position *u*,

> Xcov(u) = (1/N) · X(:, u−N+1 … u) · X(:, u−N+1 … u)ᵀ ,

which is exactly what the moving average yields at column *u*: a trailing
mean of outer products over the last N samples equals the product of those N
columns with their own transpose, divided by N. Before the start of the
recording the original pads with zeros, and so does this version.

Everything else follows the original: the IIR filter and its state, the
mirrored extension at the front by the lookahead P and the matching trim, the
eigendecomposition with eigenvalues sorted ascending, the rule for which
directions to keep (variance below the calibrated threshold, or among the
lowest-variance directions protected by `maxdims`), the reconstruction matrix
and its identity case, and the cosine blend since the last update. A `splits`
argument reproduces the original's per-block update grid, which exists only
to make a like-for-like comparison possible; the default is a single grid
over the recording.

Three details of the wrapper `clean_asr` had to be matched exactly, each of
which produced a clearly wrong result until it was: the calibration is called
through the plugin's disk cache *without* `maxmem` (inside `asr_calibrate`
that value sets the block size over which covariance matrices are averaged,
so a different value yields a different **M** and **T**); a mirrored tail is
appended before processing, without which the output is early by P samples
and the correlation with the original collapses to about zero; and the
leading carry columns are deleted afterwards, without which the output is P
samples too long. Calibration itself still runs through the plugin, so
`clean_rawdata` remains a requirement.

## Validation

All comparisons used one recording: 126 channels after channel removal, 190 s
at 500 Hz, with an identical calibration state supplied to both
implementations.

The covariance matrices themselves agree to **1 × 10⁻¹³** relative
difference, confirming that the direct product and the smoothed series are
the same quantity. With the update grid matched to the original's block
layout, the cleaned data agree at **correlation 1.000000 with an RMS
difference of 4 × 10⁻¹⁰ µV**, both at 95 blocks (its 1 GB layout) and at 5905
blocks (its 64 MB layout).

On a smaller test (16 channels, 20 s) the agreement at a matched grid was
correlation 0.999998 rather than exactly 1, with an RMS difference of 0.19 %
of the signal. The cause is floating point: the original accumulates a
cumulative sum over thousands of columns and drifts slightly, so where an
eigenvalue lies at the threshold the two implementations can make opposite
keep-or-drop decisions, changing a few short stretches. The direct
computation is the more accurate of the two, as each covariance is formed
from scratch.

Run with its own default settings — one update grid instead of the
memory-determined one — the new version differs from the plugin as much as
the plugin differs from itself across memory settings: correlation 0.956
against the plugin at 1 GB, where the plugin at 64 MB against the plugin at
1 GB gives 0.965. In other words, the deviation introduced by removing the
blocks is of the same size as the arbitrary variation the blocks were already
causing.

## Results

| Part of the algorithm | Original procedure or value | New procedure or value | Speed-up |
|---|---|---|---|
| Local covariance | Outer product at every sample, smoothed by a moving average over C² × S (cumulative sum, memory-bound) | One C × C matrix product per update position (BLAS, multithreaded) | — (the change itself) |
| Peak intermediate memory | 12.6 GB (C² × S × 8 bytes) | A few MB (C × C per position, C × N of state) | ~10⁴× less |
| Columns computed vs used | S computed, S/step used (124 of 125 discarded) | S/step computed, all used | 125× less work in that step |
| Memory blocks | 5905 at the 64 MB default; 100 at 1 GB | None | — |
| Update grid | Restarts in every block; depends on memory | One grid over the recording | — |
| Processing step (`asr_process`) | 24.3 s at 64 MB; 3.4 s at 1 GB | 1.1 s | **22×** vs default, 3.1× vs 1 GB |
| ASR step in the application | 11.4 s | 4.7 s | **2.4×** |
| Reproducibility | Output depends on memory: 26 % RMS difference between 64 MB and 1 GB (r = 0.965) | Identical on every run | — |
| Agreement, matched grid | — | r = 1.000000, RMS difference 4 × 10⁻¹⁰ µV | — |
| Calibration (`asr_calibrate`) | Unchanged, through the plugin | Unchanged, through the plugin | 1× |
| Filter, eigendecomposition, keep rule, blending, edges | Unchanged | Unchanged | 1× |

The remaining time in the application's ASR step is no longer ASR: channel
detection (`clean_channels`, RANSAC) takes 6.5 s of it, which is a separate
target.

## Limitations

The implementation depends on `clean_rawdata` for calibration and on the
plugin's behaviour staying as it is; the three wrapper details above are not
documented by the plugin and were established by comparison. The agreement
shown here is for one recording and one parameter set, the one the
application uses. The borderline keep-or-drop differences described above
mean that exact equality cannot be guaranteed in general, only agreement to
numerical precision in the quantities that enter the decision. Users
reproducing a published `clean_rawdata` analysis should switch the fast
option off: the DISCOVER-EEG preset does so by default, since that pipeline
calls `pop_clean_rawdata` itself.

## Use in the application

The ASR step offers **Fast ASR (deterministic)**, on by default, and
**Working memory (GB)**, default 1, which is passed to `clean_rawdata` for
the parts that still run through the plugin (0 leaves it the regular 64 MB). With
the fast option on, the step follows the order of `clean_artifacts`: flatline,
drift, channel and line-noise tests, then ASR, then the window criterion;
burst rejection marks the samples ASR changed, as `clean_artifacts` does.

---

*The original ASR is by Christian Kothe (`clean_rawdata`; Mullen et al.,
IEEE TBME 2015). This reimplementation keeps that algorithm and changes only
how the local covariance is obtained.*
