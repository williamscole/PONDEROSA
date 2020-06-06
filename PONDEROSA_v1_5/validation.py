import statistics as stat
po_data = open("BAGS_PO.txt").readlines()[1:]
fam = open("/share/hennlab/projects/himba_pedigree/BAGS_PONDEROSA/BAGS_convertedID.fam").readlines()

mom_dict = dict([[lines.split()[1],lines.split()[2]] for lines in fam])
dad_dict = dict([[lines.split()[1],lines.split()[3]] for lines in fam])
strength = {lines.split()[0]:float(lines.split()[-1]) for lines in po_data}
def ID(iid1,iid2):
	return min(iid1,iid2) + "_" + max(iid1,iid2)

fail,success = 0,0
f_strengths,s_strengths= [],[]
for lines in po_data:
	parent,child = lines.split()[9:11]
	if child in mom_dict and mom_dict[child] == parent:
		success += 1
		s_strengths.append(strength[ID(parent,child)])
	if parent in mom_dict and mom_dict[parent] == child:
		fail += 1
		f_strengths.append(strength[ID(parent,child)])
	if child in dad_dict and dad_dict[child] == parent:
		success += 1
		s_strengths.append(strength[ID(parent,child)])
	if parent in dad_dict and dad_dict[parent] == child:
		f_strengths.append(strength[ID(parent,child)])
		fail += 1
print(success,fail,stat.mean(s_strengths),stat.mean(f_strengths))
