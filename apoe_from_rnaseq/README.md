# APOE Genotyping from RNA-seq Data

## Overview

`apoe_from_rnaseq.sh` is a Bash script that determines APOE (Apolipoprotein E) genotype from RNA-seq BAM files by analyzing two critical SNP positions that define the ε2, ε3, and ε4 alleles. The script performs targeted variant calling at these specific genomic positions and outputs results in both VCF and TSV formats.

## Background

APOE genotype is determined by two SNPs on chromosome 19:
- **rs429358** (position 44,908,684 in GRCh37/hg19): C→T variant
- **rs7412** (position 44,908,822 in GRCh37/hg19): C→T variant

These combinations define the three major APOE alleles:
- **ε2**: rs429358=T, rs7412=T
- **ε3**: rs429358=T, rs7412=C (most common)
- **ε4**: rs429358=C, rs7412=C 

## Requirements

### Software Dependencies
- **samtools** (≥1.10 recommended)
- **bcftools** (≥1.10 recommended)
- Standard Unix utilities: `awk`, `grep`, `sed`, `sort`

### Input Files
1. **BAM file**: Aligned RNA-seq reads (coordinate-sorted)
2. **Reference FASTA**: Genome reference file (**must be identical to the reference used for mapping**)
3. **Output prefix**: Base name for output files

⚠️ **Critical**: The reference FASTA must be the exact same file used during the original read alignment. Using a different reference (even of the same genome build from a different source) will cause coordinate mismatches, incorrect variant calls, and failure of the quality control system.

## Reference Genome Downloads

### Common Reference Genomes



**GRCh38/hg38**:
```bash
# Ensembl GRCh38 (recommended)
wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Alternative: UCSC hg38
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

⚠️ **Critical**: Use the exact same reference source (Ensembl vs UCSC) and version that was used for your original RNA-seq alignment. The script's quality control system will detect reference mismatches.

## Usage

```bash
./apoe_from_rnaseq.sh <bam_file> <reference_fasta> <output_prefix>
```

### Parameters
- `<bam_file>`: Path to coordinate-sorted BAM file
- `<reference_fasta>`: Path to reference genome FASTA file
- `<output_prefix>`: Prefix for output files (no extension)

### Example
```bash
./apoe_from_rnaseq.sh sample1.bam hg19.fa sample1_apoe
```

## Features

### Automatic Index Generation
- Creates BAM index (`.bai`) if missing
- Creates FASTA index (`.fai`) if missing

### Chromosome Naming Compatibility
- Automatically detects chromosome naming convention (`chr19` vs `19`)
- Supports both styles in BAM and FASTA files
- Validates compatibility between input files
- Exits with error if chromosome styles are incompatible

### Quality Control System
The script includes **upstream control positions** (one base before each APOE SNP) in the BED intervals as an internal quality control mechanism:
- **rs429358** (pos 44,908,684): Includes upstream control at pos 44,908,683 
- **rs7412** (pos 44,908,822): Includes upstream control at pos 44,908,821

⚠️ **Critical QC Design**: These upstream positions should be **completely invariant** (no variants called) across all samples. This is a clever quality control feature that provides immediate verification of reference genome integrity.

**If variants appear at these control positions**, it indicates:
- Reference genome mismatch with the original mapping reference
- Coordinate system errors  
- Wrong genome build
- Corrupted reference file

**Expected behavior**: Only the actual SNP positions (44,908,684 and 44,908,822) should show variants. The upstream control positions must remain invariant to confirm proper reference alignment.

### Quality Filtering
- **Mapping Quality**: ≥20 (excludes poorly mapped reads)
- **Base Quality**: ≥20 (excludes low-quality base calls)
- **No depth downsampling**: Processes up to 100,000 reads per position

### Comprehensive Variant Information
- Reports allele depth (AD) for reference and alternate alleles
- Includes strand-specific allele depths (ADF, ADR)
- Provides total depth (DP) at each position
- Generates standard VCF format for downstream analysis

## Output Files

### 1. VCF File (`<prefix>.apoe.ALL.vcf.gz`)
- Compressed VCF containing variant calls at both APOE positions
- Includes all samples and comprehensive annotation
- Indexed with CSI format for rapid access

### 2. Index File (`<prefix>.apoe.ALL.vcf.gz.csi`)
- CSI index for the VCF file
- Enables efficient random access to variant data

### 3. Summary Table (`<prefix>.apoe.sites.withheader.tsv`)
- Tab-delimited file with clean, parsed genotype information
- **Columns**:
  - `chrom`: Chromosome name
  - `pos`: Genomic position
  - `ref`: Reference allele
  - `alt`: Alternate allele(s)
  - `gt`: Genotype call (0/0, 0/1, 1/1, etc.)
  - `dp`: Total read depth
  - `ad`: Allele depths (ref,alt1,alt2,...)

## Technical Details

### Genomic Coordinates (GRCh37/hg19)
```
rs429358: chr19:44,908,684 (0-based: 44,908,683-44,908,684)
rs7412:   chr19:44,908,822 (0-based: 44,908,821-44,908,822)

Quality Control Design:
Upstream control rs429358: chr19:44,908,683 (must show no variants)
Upstream control rs7412:   chr19:44,908,821 (must show no variants)
```

### BED File Structure
```bash
# Build both BED styles
bed_nochr="$(mktemp)"; cat > "$bed_nochr" <<EOF
19	44908683	44908684	rs429358
19	44908821	44908822	rs7412
EOF
bed_chr="$(mktemp)"; cat > "$bed_chr" <<EOF
chr19	44908683	44908684	rs429358
chr19	44908821	44908822	rs7412
EOF
```

**Important**: The BED coordinates include **upstream control positions** (one base before each SNP) as a quality control mechanism. These positions should show no variants called and serve as a reference integrity check.

### Variant Calling Parameters
- **bcftools mpileup**:
  - `-f`: Reference FASTA file
  - `-q 20`: Minimum mapping quality
  - `-Q 20`: Minimum base quality
  - `-a AD,ADF,ADR,DP`: Include allele depth annotations
  - `-d 100000`: Maximum depth (no downsampling)
  - `-R`: Restrict to specified regions (BED file)

- **bcftools call**:
  - `-m`: Multiallelic caller
  - `-A`: Keep all possible alternate alleles
  - `-Oz`: Compressed VCF output

### Error Handling
- **Exit Code 1**: Insufficient command-line arguments
- **Exit Code 2**: Missing required software (samtools/bcftools)
- **Exit Code 3**: Incompatible chromosome naming between BAM and FASTA

## Interpreting Results

### Genotype Calls
- **0/0**: Homozygous reference
- **0/1**: Heterozygous
- **1/1**: Homozygous alternate
- **./.**:  No call (insufficient coverage/quality)

### APOE Allele Determination
Use the genotype calls at both positions to determine APOE status:

| rs429358 | rs7412 | APOE Allele |
|----------|--------|-------------|
| T/T      | T/T    | ε2/ε2       |
| T/T      | C/T    | ε2/ε3       |
| T/C      | C/C    | ε3/ε4       |
| T/T      | C/C    | ε3/ε3       |
| C/C      | C/C    | ε4/ε4       |
| T/C      | T/C    | ε2/ε4       |

### Quality Assessment
- **High confidence**: DP ≥ 10, balanced allele depths
- **Medium confidence**: DP 5-9, some allelic imbalance acceptable
- **Low confidence**: DP < 5, consider manual review

## Limitations

1. **Reference genome dependency**: Coordinates are specific to GRCh37/hg19
2. **RNA-seq specific**: Designed for transcriptomic data (may need adjustment for WGS/WES)
3. **Expression dependency**: Requires APOE expression in the sample
4. **Allelic imbalance**: RNA-seq may show allele-specific expression effects

## Troubleshooting

### Common Issues
1. **"Missing samtools/bcftools in PATH"**
   - Install required software or add to PATH

2. **"Chromosome styles differ"**
   - Ensure BAM and FASTA use same naming (chr19 vs 19)
   - Use `samtools view -H` to check BAM chromosome names
   - Use `grep "^>" reference.fa | head` to check FASTA names

3. **No variants called**
   - Check if APOE is expressed in your sample
   - Verify alignment quality in the APOE region
   - Consider lowering quality thresholds for low-coverage samples
   - **Ensure reference genome matches**: Verify you're using the same reference FASTA that was used for initial read alignment

4. **Permission denied**
   - Make script executable: `chmod +x apoe_from_rnaseq.sh`

5. **Variants at control positions**
   - **Reference mismatch**: Most common cause is using a different reference genome than the one used for mapping
   - Check reference genome version and source (Ensembl vs UCSC)
   - Verify chromosome naming conventions match between all files
   - **Check quality control positions**: The upstream control positions should show NO variants called - if variants appear at these positions, you have a reference genome mismatch


⚠️ **Platform Note**: This script has been tested exclusively on macOS. While it should work on Linux systems, functionality on other platforms has not been verified.

## Citation
If using this script, please :

Arbones-Mainar, Jose M.(2025). APOE Genotyping from RNA-seq Data. [Computer software]. GitHub.
https://github.com/arbones/functions/tree/main/apoe_from_rnaseq