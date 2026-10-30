# A spurious comb in the spectrum after resampling to a non-dividing rate

*EegWorkflowApp technical note, October 2026*

## Summary

Resampling 1024 Hz recordings to 250 Hz with `pop_resample` adds a comb of
spurious peaks at 24, 48, 72, 96 and 120 Hz, 2 to 4 dB above the surrounding
spectrum. The peaks are not in the original data and do not appear when the
same recording is resampled to 256 Hz. They are images of the very strong
low-frequency content, folded down by the decimation and passed by the finite
stopband of the anti-aliasing filter. The comb sits at multiples of
`mod(fs_in, fs_out)`, which is 24 Hz for this pair, so the first tooth falls
inside the beta band. What is being imaged is mainly the DC offset of the
recording, 5488 uV per channel in these files: resampling to a rate that
divides the original, or simply taking the per-channel mean out first,
removes the artefact completely. The application's `spline` method avoids
the comb as well, but by having no anti-aliasing filter at all, which costs
a broadband bias of up to 3 dB towards the new Nyquist frequency.

## Background

`pop_resample` turns the rate change into an exact fraction

    [p, q] = rat(fs_out/fs_in, 1e-12)

and hands it to MATLAB's `resample`, which upsamples by `p`, low-pass filters
and decimates by `q`. The filter is not `resample`'s default: EEGLAB builds a
Kaiser-windowed sinc with `firws`, cutoff `0.9` and transition width `0.2` of
the output Nyquist frequency, designed for a passband/stopband ripple of
0.002, i.e. about 54 dB of stopband attenuation.

For the rates at issue:

| | p/q | filter length | passband | stopband |
|---|---|---|---|---|
| 1024 -> 250 | 125/512 | 16 400 taps | to 100 Hz | from 125 Hz |
| 2048 -> 250 | 125/1024 | 32 796 taps | to 100 Hz | from 125 Hz |
| 1024 -> 256 | 1/4 | 130 taps | to 102.4 Hz | from 128 Hz |
| 2048 -> 256 | 1/8 | 258 taps | to 102.4 Hz | from 128 Hz |

The filter length differs so much because the filter is designed at the
upsampled rate, `p` times the input rate.

## The observation

One recording of the Necker-cube set (`100_ec1.bdf`, 128 channels at 1024 Hz,
188 s, average reference, no other processing) was resampled to 250 and to
256 Hz, and the power spectrum averaged over channels was compared with the
original. A peak appears at 24 Hz in the 250 Hz version only. It is absent
from the 1024 Hz original and from the 256 Hz version.

## Mechanism

Upsampling by `p` replicates the spectrum at every multiple of the input
sampling rate: at 1024, 2048, 3072 Hz and so on. The interpolation filter is
meant to remove those images, but it only attenuates them by about 54 dB. In
EEG the power below 1 Hz stands roughly 35 dB above the level at 24 Hz (26.5
against -8 dB in this recording), so what survives the stopband is still
visible. Decimation then folds each image down into the new band, and the
image of the low-frequency peak at `k * fs_in` lands at

> f_k = | k * fs_in - n * fs_out |, taken into [0, fs_out/2]

For 1024 and 250 Hz: 1024 - 4*250 = 24, 2048 - 8*250 = 48, and so on. Since
`k * 1024 mod 250 = 24k mod 250`, the artefact is a comb at multiples of
`mod(fs_in, fs_out) = 24` Hz.

This is imaging, not ordinary aliasing: no content of the original recording
above 125 Hz is involved. The filter does its anti-aliasing job; what leaks
is the copy of the signal that the upsampling created.

## Evidence

**The comb is where the arithmetic puts it.** Peak height above the local
spectrum, 250 Hz output:

| on the comb | | off the comb | |
|---|---|---|---|
| 24 Hz | +2.10 dB | 36 Hz | -0.02 dB |
| 48 Hz | +4.38 dB | 60 Hz | -0.13 dB |
| 72 Hz | +4.31 dB | 84 Hz | +0.09 dB |
| 96 Hz | +1.83 dB | | |
| 120 Hz | +2.09 dB | | |

**The fold frequency follows the output rate.** Resampling the same recording
to other rates moves the peak exactly as the formula predicts:

| output | p/q | predicted | measured there | at 24 Hz |
|---|---|---|---|---|
| 250 Hz | 125/512 | 24 Hz | +2.10 dB | +2.10 dB |
| 240 Hz | 15/64 | 64 Hz | +5.24 dB | -0.65 dB |
| 200 Hz | 25/128 | 24 Hz | +0.34 dB | +0.34 dB |
| 256 Hz | 1/4 | none | - | -0.67 dB |

At 240 Hz the peak moves to 64 Hz and 24 Hz is clean, which rules out any
physiological origin.

**It is driven by the low-frequency power, and above all by the DC offset.**
The raw Biosemi files carry a mean offset of 5488 uV per channel. Taking the
per-channel mean out before resampling is enough to remove the comb
completely, and a 1 Hz high-pass does the same:

| | power 0.1-1 Hz | 24 Hz | 48 Hz |
|---|---|---|---|
| as recorded | 26.5 dB | +2.10 dB | +4.38 dB |
| per-channel mean removed | - | -0.72 dB | -0.17 dB |
| 1 Hz high-passed | 5.7 dB | -0.72 dB | -0.18 dB |

With the DC removed, `pop_resample` reproduces the original spectrum to
0.00 dB in every band below 100 Hz.

**It is stopband leakage.** Tightening the anti-aliasing filter,
`pop_resample(EEG, 250, 0.8, 0.1)`, roughly halves it: 24 Hz from +2.10 to
+1.12 dB, 48 Hz from +4.38 to +1.15 dB.

**It is not an edge or padding effect.** The peak is stationary over the
recording; discarding data at both ends does not reduce it:

| data dropped at each end | 24 Hz | 48 Hz |
|---|---|---|
| none (188 s) | +2.10 dB | +4.38 dB |
| 2 s (184 s) | +2.10 dB | +4.40 dB |
| 10 s (168 s) | +2.17 dB | +4.49 dB |
| 30 s (128 s) | +2.24 dB | +4.38 dB |

## Why 256 Hz is immune

For 1024 -> 256, p/q = 1/4 and p = 1: there is no upsampling, so no images
are created and nothing can leak. The same holds for 2048 -> 256 (1/8). This
is also why the filter is 130 instead of 16 400 taps and the step runs about
five times faster. A target rate that divides the original avoids the
artefact by construction, not by a better filter.

## The spline alternative in the application

The Resample step offers a `spline` method for non-integer ratios, which
interpolates the unfiltered data with `interp1(..., 'spline')`. Because it
never upsamples, it creates no images and shows no comb at all. It pays for
that with aliasing: nothing removes the content above the new Nyquist
frequency, which folds back as a positive bias that grows with frequency.

Error against the 1024 Hz original, same recording, 1024 -> 250 Hz:

| band | pop_resample | spline |
|---|---|---|
| 8-13 Hz | -0.01 dB | +0.04 dB |
| 18-28 Hz | -0.01 dB | +0.13 dB |
| 28-40 Hz | 0.00 dB | +0.40 dB |
| 40-60 Hz | +0.01 dB | +0.63 dB |
| 60-80 Hz | 0.00 dB | +1.08 dB |
| 80-100 Hz | 0.00 dB | +1.66 dB |
| 100-110 Hz | -1.08 dB | +2.23 dB |
| 110-124 Hz | -10 to -30 dB (the intended roll-off) | +2.6 to +2.95 dB |

| comb, dB above the local floor | 24 Hz | 48 Hz | 72 Hz |
|---|---|---|---|
| pop_resample, DC present | +2.10 | +4.38 | +4.31 |
| pop_resample, DC removed | -0.72 | -0.17 | +0.38 |
| spline, either way | -0.68 | -0.04 | +0.50 |

So the spline trades a comb of narrow peaks for a broadband bias, and the
bias is unaffected by removing the DC because it is intrinsic to having no
filter. By 80 Hz it adds more power than the comb teeth do. The ranking on
this recording is: `pop_resample` on DC-removed data, clearly first; then
the spline and `pop_resample` with the DC left in, each wrong in a different
place.

## On the padding

`pop_resample` pads both ends with a DC constant (the first and last sample,
repeated) before filtering, and removes the padding afterwards. This was
added to stop a filtered DC offset from ringing into the data for seconds
(EEGLAB bugs 1017 and 1757). A mirror or point-reflected padding would be
slightly gentler still, since it continues the waveform instead of flattening
it, and point reflection (`2*x(1) - x(n:-1:2)`, as `clean_asr_faster` uses)
keeps both the value and the slope continuous where DC padding keeps only the
value.

It would make no difference to the artefact described here. Padding only
affects the first and last half-length of the filter, about 128 ms, whereas
this comb is stationary over the whole recording, as the table above shows.
Better padding is a reasonable thing to want for the edges; it is not a fix
for imaging.

## Recommendations

1. **Resample to a rate that divides the original.** For Biosemi recordings
   at 1024 or 2048 Hz, 256 Hz costs nothing and removes the problem at the
   source. The application and the validation scripts use 256 for this
   reason, 250 Hz being the rounder number only to human eyes.
2. **If a published pipeline requires 250 Hz, take the DC out first.**
   Subtracting the per-channel mean, or detrending, or a 1 Hz high-pass: any
   of them removes the comb, because the offset is what is being imaged.
   This is the cheapest fix and it leaves the spectrum exact to 0.00 dB
   below 100 Hz. Note that DISCOVER-EEG resamples to 250 Hz before its
   high-pass, so its data still carry the offset at that point.
3. **Otherwise tighten the filter**, `pop_resample(EEG, 250, 0.8, 0.1)`,
   which halves the artefact at the cost of a longer filter.
4. **The spline method avoids the comb but is not more accurate.** It has no
   anti-aliasing at all, so it adds a broadband bias instead (see above).
   Use it only when the comb specifically matters more than accuracy near
   the top of the band, and prefer removing the DC and using `pop_resample`.

## Limitations

The measurements are from one recording, one montage and one reference, and
the absolute size of the comb depends on how much low-frequency power the
recording carries, so a drifting or DC-coupled recording will show more and a
high-passed one none. The mechanism, however, follows from the arithmetic of
`p/q` and the filter specification, and the predicted fold frequencies were
confirmed at four output rates.

The practical point is that the first tooth at 24 Hz falls inside the beta
band. Any beta measure taken from 1024 Hz data resampled to 250 Hz without a
prior high-pass carries it.

---

*`pop_resample` is by Arnaud Delorme and Andreas Widmann (EEGLAB). Nothing in
it is wrong: the artefact is the predictable consequence of a 54 dB stopband
and a rate ratio with a large numerator, and the note exists to record when
that matters.*
