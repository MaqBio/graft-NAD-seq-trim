#!/usr/bin/env bash
set -euo pipefail

fq="${1:?Usage: $0 <*_combined_R2.fastq.gz> [out.tsv]}"
out="${2:-/dev/stdout}"

# ================== 配置区===================
P7CONST="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTC"

# 这里就是你的 motif 列表（单独的 b）
motifs=(
  "TCTTCTTGT"
  "TCTTCTTG"
  "TCTTCTT"
  "TCTTCT"
  "TCTTC"
)
# ============================================

# 更快解压：pigz > gzip
if command -v pigz >/dev/null 2>&1; then
  DECOMPRESS=(pigz -dc)     # 需要更快可改成：pigz -p 8 -dc
else
  DECOMPRESS=(gzip -dc)
fi

# 把 bash 数组传给 awk（逗号拼接）
motif_csv=$(IFS=,; echo "${motifs[*]}")

"${DECOMPRESS[@]}" "$fq" \
| awk -v motifs="$motif_csv" -v p7="$P7CONST" -v OFS="\t" '
BEGIN{
  n = split(motifs, m, ",")
  for(i=1;i<=n;i++){
    c[i]=0      # motif alone
    cp7[i]=0    # P7+motif
    p7m[i]=p7 m[i]
  }
  total=0
}
NR%4==2{
  total++
  seq=$0
  for(i=1;i<=n;i++){
    if(index(seq, m[i])>0)   c[i]++
    if(index(seq, p7m[i])>0) cp7[i]++
  }
}
END{
  print "fq", "'"$fq"'"
  print "total_reads", total
  print "motif", "motif_n", "motif_pct", "P7+motif_n", "P7+motif_pct"
  for(i=1;i<=n;i++){
    pct   = (total>0)? (100.0*c[i]/total)   : 0
    pctp7 = (total>0)? (100.0*cp7[i]/total) : 0
    printf "%s\t%d\t%.4f\t%d\t%.4f\n", m[i], c[i], pct, cp7[i], pctp7
  }
}' > "$out"
