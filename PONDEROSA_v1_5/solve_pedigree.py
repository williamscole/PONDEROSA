import pandas as pd
import pedigree_graphs as ped
import numpy as np
pedigree = ped.Pedigree()
pedigree.run_PONDEROSA("King_Relatedness_no2278.seg","Himba_allPO.fam","Himba",False)

second = pd.read_csv("Himba_second.txt",delim_whitespace=True)

print(second)

#If 2 individuals are HS but have different dummy/ungtd parents
def merge_parent(iid1,iid2,sex):
	print(iid1,iid2,sex)
	parent1,parent2 = pedigree.get_parent(iid1,sex),pedigree.get_parent(iid2,sex)
	plist = [parents[0] for parents in [parent1,parent2] if parents != []]
	print(len(plist))
	if len(plist) == 0:
		parent = pedigree.get_dummy_id()
	else:
		parent = plist[0]
		if len(plist) == 1:
			if parent1 == []:
				pedigree.add_po(parent,iid1,sex,True)
			else:
				pedigree.add_po(parent,iid2,sex,True)
		elif len(plist) == 2:
			rm_parent = plist[1]
			print(pedigree.get_pedigree_structure()[rm_parent][2])
			for children in pedigree.get_pedigree_structure()[rm_parent][2]:
				pedigree.add_po(parent,children[0],sex,True)

def add_halfsibs(df):
	for index,rel in enumerate(["PHS","MHS"]):
		for pairs in df[df["REL"] == rel][["OLDER","YOUNGER"]].values.tolist():
			merge_parent(pairs[0],pairs[1],index+1)
			print(pairs)

add_halfsibs(second)
