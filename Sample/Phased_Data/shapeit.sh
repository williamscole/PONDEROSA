for i in {1..22} 
do
shapeit --input-bed [.bed file prefix] -M [Genetic map file] -O Sample_shapeit.chr$i --window 5 --main 30 --prune 16 --burn 10 --thread 32  --noped
done
