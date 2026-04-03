#!/bin/bash

# =========== Do not change this code ===========
# Script name must match the IBD caller format (e.g. hap-ibd.sh, phasedibd.sh)
# PONDEROSA uses the filename (minus .sh) to determine the IBD output format.

dirpath=${1}
prefix=${dirpath}/simulation

# Find the simulated VCF (ped-sim may output .vcf or .vcf.gz)
if [[ -f "${prefix}.vcf.gz" ]]; then
    vcffile="${prefix}.vcf.gz"
elif [[ -f "${prefix}.vcf" ]]; then
    vcffile="${prefix}.vcf"
else
    echo "ERROR: No simulated VCF found at ${prefix}.vcf or ${prefix}.vcf.gz"
    exit 1
fi
# ================================================

# ========== Edit these to your config ==========

basedir="/Users/cwilli50/Desktop/projects/"
hapibd_jar="$basedir/PONDEROSA/ponderosa/simulation/hap-ibd/hap-ibd.jar"

chr_map_dir="$basedir/PONDEROSA/ponderosa/simulation/hap-ibd/plink.GRCh37.map"

# Memory in gigs
mem=8

# No. of cpus to use
nthreads=1

# ================================================

# ======== Checking that files exist =============
if [[ ! -f $hapibd_jar ]]; then
    echo "ERROR: hap-ibd jar not found: $hapibd_jar"
    exit 1
fi

if [[ ! -d $chr_map_dir ]]; then
    echo "ERROR: map directory not found: $chr_map_dir"
    exit 1
fi
# ================================================

# Create the concatenated map file for hap-ibd
# PONDEROSA expects this file as "genetic.map" in plink format (chr, rsid, cM, bp)
cat "${chr_map_dir}/plink.chr1.GRCh37.map" > "${dirpath}/genetic.map"
for i in {2..22}; do
    cat "${chr_map_dir}/plink.chr${i}.GRCh37.map" >> "${dirpath}/genetic.map"
done

java -Xmx${mem}g -jar $hapibd_jar \
    gt=$vcffile \
    map=${dirpath}/genetic.map \
    nthreads=$nthreads \
    out=$prefix