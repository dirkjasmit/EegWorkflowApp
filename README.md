# EegWorkflowApp

###Function
The EegWorkflowApp provides a simplified method for running EEG preprocessing and cleaning steps in a semi-automated or fully automated manner. Follow the buttons on the left from top to bottom for standard data reading and preprocessing, next to the middle button column for additional cleaning. The right column provides options for viewing, reviewing, epoching and saving.

Settings are saved between sessions

### Batch mode
The batch button allows automated batch processing using the settings as visible in the main window. 

### History

## version 0.9
The EegWorkflowApp is the AppDesigner version of the GUIDE based EegWorkflow. That version is obsolete and not supported since MATLAB 2025 or so.

## version 1.0 (Legacy)
Find the legacy version with the old layout here https://github.com/dirkjasmit/EegWorkflowApp/tree/Legacy
features:
- Manual or batch processing of EEG data files in various formats
- Flexible batch processing order
- Batch and manual processing share settings visible in the main window.
- Settings are saved on close (wait for the old settings to reappear after opening! Takes a while.)
- Tracking power values for each channel after each step in EEG.etc.<Power>
- loads of information in a "stdout" tracking listbox

The version before the revamp is on the [`Legacy` branch](https://github.com/dirkjasmit/EegWorkflowApp/tree/Legacy)
([download as zip](https://github.com/dirkjasmit/EegWorkflowApp/archive/refs/heads/Legacy.zip)),
or clone it with `git clone -b Legacy https://github.com/dirkjasmit/EegWorkflowApp.git`.

## version 1.1 (Revamp)
The current version with the novel layout is in the main branch. Added features are:
- Increased parametrization of the cleaning steps, including choice of ICA algorithm and number of 
- Added clustering of buttons
- Reorganized layout into main processing steps, additional processing / segment deletion, viewing & saving
- Added information window on main analysis steps
- Added cleaning step tracking with average PSD plots w/ SD around it
- Hugely optimized code for way faster preprocessing, like eeg_linenoise() which yields much better results besidess being way faster than Cleanline.


And, there's an undo button for backtracking analyses.

### Install
Install EEGLAB
From the EEGLAB window, choose manage extensions. Install:
- Biosig data import
- BDFimport (data import)
- neuroscanio (data import)
- ANTeepimport (data import)
- IClabel
- clean_rawdata
- PICARD
- firfilt
- AAR

For faster ICA decomposition, install and test binica. PICARD is 

Download the full github.

To run, start the appdesigner, open the guiEegWorkflow.mlapp and click the RUN button (guiEegAutoflow_App.mlapp in the legacy version)


## Copyright statements for included code

### REST reference

All copyrights of the software are reserved by the Key Laboratory for NeuroInformation of Ministry of Education, School of Life Science and Technology, University of Electronic Science and Technology of China. This software is for non-commercial use only. It is freeware but not in the public domain.

For more see http://www.neuro.uestc.edu.cn/rest/
Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
                      Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305

Written by Li Dong (Li_dong729@163.com) and Shiang Hu (hushiang@126.com)
Date: Aug. 19, 2019
% ---------------------------------------------------------------------------------------------
######################Please cite this toolbox as:###############################

Li Dong*, Fali Li, Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. MATLAB Toolboxes for Reference Electrode Standardization Technique (REST) of Scalp EEG. Frontiers in Neuroscience,  2017:11(601).
% --------------------------------------------------------------------------------------------

Log

REST_reference_v1.0_20170411
      [1] Retain unselected channels (e.g. EMG, EOG, band channels etc.), while saving re-reference EEG data.

REST_reference_v1.1_20190819
      [1] Calculate leadfield at once, based on canonical concentric-three-spheres head model.
      [2] remove the leadfield.exe.
      [3] add reference papers on the interface







