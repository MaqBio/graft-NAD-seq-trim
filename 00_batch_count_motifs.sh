#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# batch_count_motifs.sh
#
# Usage:
#   bash batch_count_motifs.sh <root_dir> [out_dir]
#
# Optional env:
#   TOTAL_CPU=16   # total CPU budget you want to use (your interactive limit is 16)
#   JOBS=4         # number of samples to run in parallel
#
# It will auto set:
#   PIGZ_THREADS = floor(TOTAL_CPU / JOBS), at least 1
#
# Example (recommended for your 16-core interactive):
#   TOTAL_CPU=16 JOBS=4 bash batch_count_motifs.sh /labs/lilab/maq/graft/dxo1 motif_stats
# -----------------------------

ROOT="${1:?Usage: $0 <root_dir> [out_dir]}"
OUTDIR="${2:-motif_stats}"

# ===== 配置（可用环境变量覆盖）=====
JOBS="${JOBS:-4}"
TOTAL_CPU="${TOTAL_CPU:-16}"
# ==================================

mkdir -p "$OUTDIR"

mapfile -t R2S < <(find "$ROOT" -type f -name "*_combined_R2.fastq.gz" | sort)

if (( ${#R2S[@]} == 0 )); then
  echo "No *_combined_R2.fastq.gz found under: $ROOT" >&2
  exit 1
fi

# 自动给每个样本分配 pigz 线程：TOTAL_CPU / JOBS
PIGZ_THREADS=$(( TOTAL_CPU / JOBS ))
if (( PIGZ_THREADS < 1 )); then PIGZ_THREADS=1; fi

echo "Found ${#R2S[@]} R2 files under: $ROOT"
echo "Output dir: $OUTDIR"
echo "TOTAL_CPU=$TOTAL_CPU  JOBS=$JOBS  => PIGZ_THREADS(per sample)=$PIGZ_THREADS"
echo

# 并行处理每个 R2
printf "%s\n" "${R2S[@]}" \
| xargs -I{} -P "$JOBS" bash -lc '
  fq="{}"
  sample=$(basename "$fq" | sed "s/_combined_R2\.fastq\.gz//")
  out="'"$OUTDIR"'/${sample}.motif.tsv"
  PIGZ_THREADS="'"$PIGZ_THREADS"'" bash '"$(pwd)"'/count_motifs_onepass.sh "$fq" "$out"
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
