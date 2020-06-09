import time
import random
import haplotype_scores as haps
king = open("/share/hennlab/projects/himba_pedigree/ALGORITHM/King_Relatedness_no2278.seg").readlines()[1:]

def pair_id(line):
	return min([line.split()[1],line.split()[3]]) + "_" + max([line.split()[1],line.split()[3]])

relatives = [pair_id(line) for line in king if float(line.split()[6]) > 0.30]

def run(num):
	start = time.time()
	l = random.sample(relatives,num)
	a = haps.get_hap_score(l,{"match_file":"/share/hennlab/projects/himba_pedigree/Germline/Himba_germline_chr1.output.match","map_file":"/share/hennlab/projects/himba_pedigree/plink_data/newHimba_shapeit.chr1.map","ped_file":"None","out":"test","num_chr":1,"cm_gap":1,"disc_homoz":1},"None")
	return time.time()-start
out = open("run_time.txt","w")	
for i in range(100,3200,200):
	times = []
	for iter in range(50):
		times.append(run(i))
	out.write("%s %s\n" % (i,sum(times)/len(times)))
	print(i)
