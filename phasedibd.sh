module load plink/2.00
pop=${1}
rel=${2}
iter=${3}

for chrom in {1..22}
do

plink2 --vcf ${pop}/${rel}${iter}.vcf --chr ${chrom} --recode vcf --out ${pop}/${rel}${iter}_chr${chrom}

python pedigree_tools.py ${pop}/sim_chr${chrom}.map $rel $iter $pop $chrom pop_ibd

done