# EEG preprocessing pipelines and what EegWorkflow needs to emulate them

Citation counts from OpenAlex (looked up 2026-09-22; Google Scholar blocked
automated lookup). OpenAlex counts typically run 20–50% below Google Scholar,
but the ranking is comparable. Sorted by citations.

| Pipeline (paper) | Cites | What's particular about it | What EegWorkflow needs to emulate it |
|---|---|---|---|
| **PREP** (Bigdely-Shamlo 2015, [doi](https://doi.org/10.3389/fninf.2015.00016)) | 1468 | Standardisation for large datasets: line-noise removal, then a **robust average reference**. It repeats detect bad channels → interpolate → re-reference until the set of bad channels stops changing. Detection uses four tests: robust deviation, correlation in 1 s windows, RANSAC predictability, high-frequency noise. | **RANSAC** and **high-frequency-noise** tests in Bad chans; an iterative **"Robust average (PREP)"** option in Rereference. Deviation and correlation tests already exist. |
| **FASTER** (Nolan 2010, [doi](https://doi.org/10.1016/j.jneumeth.2010.07.015)) | 1147 | Z-score > 3 on several features at four levels: channels (variance, correlation, Hurst exponent); epochs (amplitude range, variance, deviation); ICs (EOG correlation, spatial kurtosis, spectral slope, Hurst, median gradient); **channels within epochs**, which are interpolated per epoch. | **Epoch-level operations** (reject epochs, interpolate a channel inside one epoch); feature-based IC rejection as an alternative to ICLabel. |
| **Autoreject** (Jas 2017, [doi](https://doi.org/10.1016/j.neuroimage.2017.06.030)) | 672 | MNE/Python. **Learns a peak-to-peak threshold per channel by cross-validation**; per epoch it interpolates up to ρ bad channels and drops the epoch if more than κ are bad. | Epoching plus per-epoch interpolation (as for FASTER), plus the cross-validated threshold search. Largest effort of all. |
| **HAPPE** (Gabard-Durnam 2018, [doi](https://doi.org/10.3389/fnins.2018.00097)); HAPPE+ER (2022, 59); HAPPILEE (2022, 65) | 583 | Built for short, noisy, developmental recordings (HBN-like). **Wavelet thresholding** of the artefacts (in v2 on the channel data), ICA with MARA or ICLabel, and a **quality report per file** (% data kept, % channels, artefact-to-signal measures). Variants for ERPs and for low channel counts. | A **wavelet-thresholding step**; a **per-file quality table** written by Batch. |
| **Automagic** (Pedroni 2019, [doi](https://doi.org/10.1016/j.neuroimage.2019.06.046)) | 301 | A wrapper around PREP, clean_rawdata, **EOG regression** and MARA/ICLabel. Then an **automatic quality rating** of each file (Good/OK/Bad from high amplitude, time/channel variance and bad-channel ratio) and a viewer to check the rated files. | PREP (above); **EOG regression** (current EOG buttons are ICA-based); quality ratings in the batch output. |
| **MADE** (Debnath 2020, [doi](https://doi.org/10.1111/psyp.13580)) | 226 | Developmental data. FASTER bad channels; **ICA on a 1 Hz high-passed copy cut into 1 s epochs, with the weights copied back** to the 0.1 Hz data; **adjusted-ADJUST** for component selection; then per-epoch ±150 µV rejection and channel interpolation. | **ICA on a filtered copy with the weights transferred** (small, also useful on its own); adjusted-ADJUST as a classifier; epoch-level steps. |
| **BEAPP** (Levin 2018, [doi](https://doi.org/10.3389/fnins.2018.00513)) | 141 | Mostly infrastructure: batch processing of many formats and modules, including HAPPE and spectral output. | Already covered by batch mode and multi-format Open. |
| **RELAX** (Bailey 2023; pt 1: 94, pt 2: 32) | 126 | **Multi-channel Wiener filter (MWF)** for blinks, muscle and drift. Muscle is flagged by a **log-log spectral slope > −0.59**, close to the EMG slope test. Then **wavelet-ICA** on ICLabel artefact components. | **MWF step** (new); wavelet-ICA. The EMG slope test already matches the muscle criterion. |
| **DISCOVER-EEG** (Gil Avila 2023, [doi](https://doi.org/10.1038/s41597-023-02525-0)) | 76 | clean_rawdata channel removal, average reference, 10× ICA with the most typical run kept, ICLabel Muscle/Eye ≥ 0.8, ASR bad segments. | **Done**, apart from interpolating the removed channels back. |
| **APP** (da Cruz 2018, [doi](https://doi.org/10.1016/j.clinph.2018.04.600)) | 75 | Like FASTER, but on **robust statistics** (median/IQR instead of mean/SD) for channels, epochs and ICs. | Robust-statistics options for the z-score criteria; epoch-level steps. |
| **Lossless / EEG-IP-L** (Desjardins 2021, [doi](https://doi.org/10.1016/j.jneumeth.2020.108961)) | 71 | **Annotates instead of deleting**: channels, time and ICs are flagged iteratively (with AMICA), then checked in a **QC review**. The output is the full data plus annotations. | A **"flag only" mode** (mark as events or masks instead of removing); AMICA; the Review viewer covers the QC part. |
| **NEAR** (Kumaravel 2022, [doi](https://doi.org/10.1016/j.dcn.2022.101068)) | 55 | Newborns. Bad channels by **Local Outlier Factor (LOF)**, then ASR with tuned settings. | A **LOF bad-channel test**; the ASR button covers the rest. |

## Also relevant

- **Makoto's pipeline** (EEGLAB wiki, no paper, so no count) is probably the
  most copied recipe. EegWorkflow already covers it except AMICA and
  dipole-based component rejection.
- **ASR** itself (Mullen 2015: 847; Chang 2019: 573) is a method rather than a
  pipeline, and is in the ASR button.
- **MNE-Python** (Gramfort 2013: 4180) is a toolbox, not a pipeline.

## Missing building blocks, by how many pipelines need them

1. **Epoch-level operations:** per-epoch channel interpolation and epoch
   rejection. Needed by FASTER, Autoreject, MADE and APP; a real change in
   direction for a continuous-data app.
2. **Wavelet-based cleaning:** HAPPE and RELAX.
3. **Per-file quality metrics in Batch:** HAPPE and Automagic. Cheap, since the
   spectrum tracker already computes most of it.
4. **ICA on a filtered copy with weights transferred:** MADE; also recommended
   in Makoto's pipeline. Small.
5. **PREP robust reference with RANSAC:** PREP and Automagic.
6. **Alternative component classifiers** (MARA, adjusted-ADJUST), **EOG
   regression**, **MWF**, **LOF**, and a **flag-only mode**: one pipeline each.

Suggested order for HBN-type data: 3 and 4 first (both cheap), then wavelet
thresholding, then PREP.
