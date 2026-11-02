#!/bin/bash
# Turn relative symlinks into absolute ones.
#
# Each link's target has its leading "../" components stripped and BASE put in
# front, so
#
#   sub-X_task-RestingState_eeg.set -> ../R07_ds005511/sub-X/eeg/sub-X_...set
#
# becomes
#
#   sub-X_task-RestingState_eeg.set -> BASE/R07_ds005511/sub-X/eeg/sub-X_...set
#
# The release folder (R01..R11) is whatever the link already said: nothing is
# hard-coded beyond BASE. Links that are already absolute are left alone.
#
# Only links whose NAME matches PATTERN are considered, "sub-*" by default, so
# a directory holding other links as well (lnk_*.bdf pointing somewhere else,
# say) is safe to run this in. Links that are already absolute are skipped too.
#
# With -s, the target is instead resolved against the directory the links were
# originally made in, which is the general rule and the one to prefer: what
# "../" meant is decided by that directory, not by where the links sit now.
# For links copied out of TDBrain/eeg_ses1_restEC, whose targets read
# "../TDBRAIN-dataset/...", -s puts them back onto TDBrain/TDBRAIN-dataset/...
#
# Dry run by default; pass -f to actually rewrite the links.
#
# Usage: abs_symlinks.sh [-f] [-r] [-b BASE | -s SOURCE] [-p PATTERN] [DIR]
#   -f          apply the changes (default: only show what would happen)
#   -r          recurse into subdirectories (default: DIR itself only)
#   -b BASE     strip leading "../" and put BASE in front
#               (default: the HBN resting-state set)
#   -s SOURCE   resolve the target against SOURCE, the directory the links
#               were created in. Overrides -b
#   -p PATTERN  only links whose name matches this (default: sub-*)
#   DIR         directory holding the links (default: the current directory)

set -uo pipefail

BASE="/Volumes/FourTB_SSD/EEG/RuweData/HBN/HBN_EEG_restingstate"
SOURCE=""
APPLY=0
RECURSE=0
PATTERN="sub-*"

while getopts ":frb:s:p:h" opt; do
    case $opt in
        f) APPLY=1 ;;
        r) RECURSE=1 ;;
        b) BASE=$OPTARG ;;
        s) SOURCE=$OPTARG ;;
        p) PATTERN=$OPTARG ;;
        h) sed -n '2,36p' "$0"; exit 0 ;;
        \?) echo "unknown option -$OPTARG" >&2; exit 2 ;;
        :)  echo "-$OPTARG needs an argument" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))

DIR=${1:-.}
BASE=${BASE%/}                      # no trailing slash

if [[ ! -d $DIR ]]; then
    echo "not a directory: $DIR" >&2; exit 1
fi
if [[ -n $SOURCE ]]; then
    SOURCE=${SOURCE%/}
    if [[ ! -d $SOURCE ]]; then
        echo "source directory does not exist: $SOURCE" >&2; exit 1
    fi
elif [[ ! -d $BASE ]]; then
    echo "*** warning *** base path does not exist: $BASE" >&2
fi

maxdepth=1                          # an empty array trips bash 3.2 under set -u
[[ $RECURSE -eq 1 ]] && maxdepth=999

undo="$DIR/.abs_symlinks_undo_$(date +%Y%m%d_%H%M%S).sh"
n_total=0; n_abs=0; n_done=0; n_broken=0; n_odd=0; n_same=0

[[ $APPLY -eq 1 ]] && printf '#!/bin/bash\n# restores the link targets as they were before abs_symlinks.sh ran\nset -e\n' > "$undo"

while IFS= read -r -d '' link; do
    n_total=$((n_total + 1))
    target=$(readlink "$link")

    if [[ $target == /* ]]; then
        n_abs=$((n_abs + 1))
        continue
    fi

    if [[ -n $SOURCE ]]; then
        # resolve against the directory the links were made in: cd into the
        # target's parent so the shell collapses "..", then keep the filename
        parent=$(cd "$SOURCE/$(dirname "$target")" 2>/dev/null && pwd -P)
        if [[ -z $parent ]]; then
            printf '*** no such directory *** %s\n                         ->  %s\n' \
                "$link" "$SOURCE/$target"
            n_broken=$((n_broken + 1))
            continue
        fi
        new="$parent/$(basename "$target")"
    else
        rest=$target
        while [[ $rest == ../* ]]; do rest=${rest#../}; done
        if [[ $rest == "$target" ]]; then
            # relative, but no leading ".." - not what this script is for
            printf 'skipped (no leading ..): %s -> %s\n' "$link" "$target"
            n_odd=$((n_odd + 1))
            continue
        fi
        new="$BASE/$rest"
    fi
    if [[ $new == "$target" ]]; then
        n_same=$((n_same + 1))
        continue
    fi

    if [[ ! -e $new ]]; then
        printf '*** missing target *** %s\n                   ->  %s\n' "$link" "$new"
        n_broken=$((n_broken + 1))
        continue
    fi

    if [[ $APPLY -eq 1 ]]; then
        if ln -sfn "$new" "$link"; then
            printf 'ln -sfn %q %q\n' "$target" "$link" >> "$undo"
            n_done=$((n_done + 1))
        else
            printf '*** failed *** %s\n' "$link" >&2
        fi
    else
        n_done=$((n_done + 1))
        [[ $n_done -le 3 ]] && printf 'would relink: %s\n          ->  %s\n' "$link" "$new"
    fi
done < <(find "$DIR" -maxdepth "$maxdepth" -type l -name "$PATTERN" -print0)

echo
printf '%-34s %6d\n' "symlinks matching '$PATTERN'" "$n_total"
printf '%-34s %6d\n' "already absolute, left alone" "$n_abs"
printf '%-34s %6d\n' "relative without .., skipped" "$n_odd"
printf '%-34s %6d\n' "already correct" "$n_same"
printf '%-34s %6d\n' "target does not exist, skipped" "$n_broken"
if [[ $APPLY -eq 1 ]]; then
    printf '%-34s %6d\n' "relinked" "$n_done"
    if [[ $n_done -gt 0 ]]; then
        chmod +x "$undo"
        echo
        echo "to undo: $undo"
    else
        rm -f "$undo"
    fi
else
    printf '%-34s %6d\n' "would be relinked" "$n_done"
    echo
    echo "dry run - nothing was changed. Re-run with -f to apply."
fi
