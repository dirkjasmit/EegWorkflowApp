#!/bin/bash
# Repair the TDBrain ses-1 restEC links in validation/datalinks.
#
# They were copied out of
#   /Volumes/FourTB_SSD/EEG/RuweData/TDBrain/eeg_ses1_restEC
# where the entries are themselves relative symlinks reading
#   ../TDBRAIN-dataset/sub-<n>/ses-1/eeg/sub-<n>_ses-1_task-restEC_eeg.vhdr
# so the copies inherited a target whose ".." no longer means anything here.
# Resolved against that source directory, ".." is TDBrain/, which is where
# TDBRAIN-dataset actually lives.
#
# Five files per subject: _channels.tsv, _eeg.eeg, _eeg.json, _eeg.vhdr,
# _eeg.vmrk. The pattern below takes all five and nothing else - the HBN
# sub-NDAR* links and the lnk_*.bdf links in the same folder are not matched.
#
# Dry run by default; pass -f to apply. Any other option is handed through.
#
# Usage: abs_symlinks_tdbrain.sh [-f] [DIR]      (DIR default: validation/datalinks)

set -uo pipefail

SOURCE="/Volumes/FourTB_SSD/EEG/RuweData/TDBrain/eeg_ses1_restEC"
PATTERN="sub-*_ses-1_task-restEC*"

here=$(cd "$(dirname "$0")" && pwd -P)
DIR="$here/../validation/datalinks"

args=()
for a in "$@"; do
    case $a in
        -*) args+=("$a") ;;
        *)  DIR=$a ;;
    esac
done

exec "$here/abs_symlinks.sh" -s "$SOURCE" -p "$PATTERN" ${args[@]+"${args[@]}"} "$DIR"
