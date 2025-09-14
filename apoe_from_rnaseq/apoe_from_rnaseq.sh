#!/usr/bin/env bash
# apoe_from_rnaseq_mac.sh
# Usage: apoe_from_rnaseq_mac.sh <bam> <fasta> <out_prefix>

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <bam> <fasta> <out_prefix>" >&2
  exit 1
fi

BAM="$1"
FA="$2"
OUT="$3"

for x in samtools bcftools; do
  command -v "$x" >/dev/null 2>&1 || { echo "Missing $x in PATH" >&2; exit 2; }
done

# Ensure indexes exist
[[ -f "${BAM}.bai" ]] || samtools index "$BAM"
[[ -f "${FA}.fai"  ]] || samtools faidx "$FA"

# Detect contig naming
bam_chr_style=$(samtools view -H "$BAM" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
fa_chr_style=$(grep -m1 '^>' "$FA" | sed 's/^>//;s/ .*//')

# Build both BED styles
bed_nochr="$(mktemp)"; cat > "$bed_nochr" <<EOF
19	44908683	44908684	rs429358
19	44908821	44908822	rs7412
EOF
bed_chr="$(mktemp)"; cat > "$bed_chr" <<EOF
chr19	44908683	44908684	rs429358
chr19	44908821	44908822	rs7412
EOF

choose_bed() {
  case "$1" in chr*) echo "$bed_chr" ;; *) echo "$bed_nochr" ;; esac
}

BED_BAM=$(choose_bed "$bam_chr_style")
BED_FA=$(choose_bed "$fa_chr_style")
if [[ "$BED_BAM" != "$BED_FA" ]]; then
  echo "Chromosome styles differ between BAM ($bam_chr_style) and FASTA ($fa_chr_style). Harmonize them." >&2
  exit 3
fi
BED="$BED_BAM"

# Variant calling at the two APOE ε sites (no downsampling; add MAPQ/BaseQ and strand-split AD)
VCF="${OUT}.apoe.ALL.vcf.gz"
bcftools mpileup -f "$FA" -q 20 -Q 20 -a AD,ADF,ADR,DP -d 100000 -R "$BED" "$BAM" \
| bcftools call -m -A -Oz -o "$VCF"
bcftools index "$VCF"

# Single tidy table with header
OUT_TSV="${OUT}.apoe.sites.withheader.tsv"
{
  echo -e "chrom\tpos\tref\talt\tgt\tdp\tad"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%DP]\t[%AD]\n' "$VCF" | sort -k2,2n
} > "$OUT_TSV"

echo "Wrote:"
echo "  $VCF"
echo "  ${VCF}.csi"
echo "  $OUT_TSV"
