#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?Usage: $0 <root_dir> [out_dir] [jobs]}"
OUTDIR="${2:-motif_stats}"
JOBS="${3:-4}"

mkdir -p "$OUTDIR"

mapfile -t R2S < <(find "$ROOT" -type f -name "*_combined_R2.fastq.gz" | sort)

echo "Found ${#R2S[@]} R2 files under: $ROOT"
echo "Output dir: $OUTDIR"
echo "Parallel jobs: $JOBS"
echo

printf "%s\n" "${R2S[@]}" \
| xargs -I{} -P "$JOBS" bash -lc '
  fq="{}"
  sample=$(basename "$fq" | sed "s/_combined_R2\.fastq\.gz//")
  out="'"$OUTDIR"'/${sample}.motif.tsv"
  bash count_motifs_onepass.sh "$fq" "$out"
  echo "[done] $sample"
'

# 汇总成一个大表（方便后续 R/Excel）
summary="$OUTDIR/summary.tsv"
{
  echo -e "sample\tmotif\tmotif_n\tmotif_pct\tP7+motif_n\tP7+motif_pct"
  for f in "$OUTDIR"/*.motif.tsv; do
    s=$(basename "$f" .motif.tsv)
    awk -v s="$s" 'BEGIN{OFS="\t"} NR>=4 {print s,$1,$2,$3,$4,$5}' "$f"
  done
} > "$summary"

echo
echo "Summary written: $summary"
