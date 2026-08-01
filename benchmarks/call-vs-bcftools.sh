#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 5 || $# -gt 7 ]]; then
    echo "usage: $0 RSOMICS_CALL BCFTOOLS REFERENCE BAM OUTPUT_DIR [REPEATS] [CPU]" >&2
    exit 2
fi

ours=$1
bcftools=$2
reference=$3
alignment=$4
output_dir=$5
repeats=${6:-5}
cpu=${7:-0}
annotations=FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/QM,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/NMBZ,INFO/NM,INFO/SCR
times=$output_dir/times.tsv

[[ $(uname -s) == Linux ]]
command -v taskset >/dev/null
/usr/bin/time --version 2>&1 | grep -q GNU
[[ $repeats =~ ^[1-9][0-9]*$ ]]
[[ $cpu =~ ^[0-9]+$ ]]
mkdir -p "$output_dir"

{
    uname -a
    lscpu
    "$ours" --version
    "$bcftools" --version
    sha256sum "$ours" "$bcftools" "$reference" "$reference.fai" "$alignment" "$alignment.bai"
} > "$output_dir/environment.txt"
printf 'round\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_kib\n' > "$times"

run_ours() {
    local round=$1
    local output=$output_dir/ours-$round.bcf
    /usr/bin/time -a -o "$times" -f "$round\tours\t%e\t%U\t%S\t%M" \
        taskset -c "$cpu" "$ours" run -f "$reference" -v -O u -o "$output" "$alignment"
}

run_bcftools() {
    local round=$1
    local output=$output_dir/bcftools-$round.bcf
    /usr/bin/time -a -o "$times" -f "$round\tbcftools\t%e\t%U\t%S\t%M" \
        taskset -c "$cpu" bash -o pipefail -c \
        '"$1" mpileup -f "$2" -a "$3" -Ou "$4" | "$1" call -mv -Ou -o "$5"' \
        _ "$bcftools" "$reference" "$annotations" "$alignment" "$output"
}

compare_calls() {
    local round=$1
    "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$output_dir/ours-$round.bcf" \
        | LC_ALL=C sort > "$output_dir/ours-$round.calls.tsv"
    "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$output_dir/bcftools-$round.bcf" \
        | LC_ALL=C sort > "$output_dir/bcftools-$round.calls.tsv"
    cmp "$output_dir/ours-$round.calls.tsv" "$output_dir/bcftools-$round.calls.tsv"
}

taskset -c "$cpu" "$ours" run -f "$reference" -v -O u \
    -o "$output_dir/ours-warmup.bcf" "$alignment"
taskset -c "$cpu" bash -o pipefail -c \
    '"$1" mpileup -f "$2" -a "$3" -Ou "$4" | "$1" call -mv -Ou -o "$5"' \
    _ "$bcftools" "$reference" "$annotations" "$alignment" "$output_dir/bcftools-warmup.bcf"
compare_calls warmup

for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        run_ours "$round"
        run_bcftools "$round"
    else
        run_bcftools "$round"
        run_ours "$round"
    fi
    compare_calls "$round"
done
