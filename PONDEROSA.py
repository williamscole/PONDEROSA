import pandas as pd
import numpy as np
import statistics as stat
import sys
import os.path
import time
import datetime
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def get_hap_score(relative_list,par_file,hap_file):
	pars = [par_file[par] for par in ["match_file","map_file","ped_file","out","num_chr","cm_gap","disc_homoz"]]
	match_file,map_file,ped_file,out,num_chr,cm_gap,disc_homoz = pars
	num_chr,cm_gap,disc_homoz = int(num_chr),float(cm_gap),int(disc_homoz)
	class GenotypeData:
		def __init__(self,map_file,ped_file,chrm):
			self.gts = dict()
			if ped_file != "None":
				for ind_gts in open(ped_file.replace("chr1","chr%s" % chrm)).readlines():
					ind_gts = ind_gts.split()
					iid,gts = ind_gts[0],ind_gts[6:]
					gt_out = list()
					for alleles in range(0,len(gts),2):
						pair=(gts[alleles],gts[alleles+1])
						gt_out.append(pair)
					self.gts[iid] = gt_out

			self.snp_pos = list()
			self.mb_cm = dict()
			for snps in open(map_file.replace("chr1","chr%s" % chrm)).readlines():
				snps = snps.split()
				cm, mb = float(snps[2]),int(snps[3])
				self.snp_pos.append(cm)
				self.mb_cm[mb] = cm

		def mb_to_cm(self,mb):
			return self.mb_cm[mb]

		def gap_discordance(self,iid1,iid2,gap_start,gap_end,threshold):
			status = False
			if self.gts != {}:
				start=self.snp_pos.index(gap_start)+1
				end=self.snp_pos.index(gap_end)
				num_disc_homoz = 0
				for index in range(start,end):
					gt1,gt2=self.gts[iid1][index],self.gts[iid2][index]
					if (gt1[0]==gt1[1]) and (gt2[0]==gt2[1]) and (gt1!=gt2):
						num_disc_homoz+=1
					if num_disc_homoz > threshold:
						status=True
						break
			return status

	class PairData:
		def __init__(self,relative_list,out):
			self.relative_list,self.out = relative_list,out
			self.pair_data = dict()
			for pairs in relative_list:
				iid1,iid2 = pairs.split("_")
				#pair maps to {hap tots}, total ibd (for hap score)
				self.pair_data[pairs] = [{iid1:0,iid2:0},0,0]

		def finish_chrm(self,pair,iid1,iid2,seg_list,num_segs):
			hap_data,tot_ibd = {iid1:{"0":0,"1":0},iid2:{"0":0,"1":0}},0
			for segs in seg_list:
				seg_len = segs[4]-segs[3]
				tot_ibd += seg_len
				hap_data[segs[1]][segs[5]] += seg_len
				hap_data[segs[2]][segs[6]] += seg_len
			self.pair_data[pair][1] += tot_ibd
			self.pair_data[pair][2] += num_segs
			self.pair_data[pair][0][iid1] += max(hap_data[iid1]["0"],hap_data[iid1]["1"])
			self.pair_data[pair][0][iid2] += max(hap_data[iid2]["0"],hap_data[iid2]["1"])

		def get_scores(self,pair):
			info = self.pair_data[pair]
			iid1,iid2 = pair.split("_")
			h1,h2,total,num = info[0][iid1],info[0][iid2],info[1],info[2]
			if total != 0:
				h1,h2 = h1/total,h2/total
				ratio = min(h1/h2,h2/h1)
			else:
				h1,h2,ratio = 0,0,0
			return [pair,iid1,h1,iid2,h2,ratio,num]

		def get_relatives(self):
			return self.relative_list

		def write_out(self,relative_list,out):
			out_df = list()
			for pairs in self.relative_list:
				out_df.append(self.get_scores(pairs))
			out_df = pd.DataFrame(out_df,columns=["PAIR_ID","IID1","H1","IID2","H2","HSR","N"])
			with open("%s.haps" % self.out,"w") as outfile:
				outfile.write(out_df.to_string(index=False))
			outfile.close()
			return out_df

	if hap_file != "None":
		sys.stdout.write("\rSkipping hap score computation.\n")
		hap_df = open(hap_file).readlines()
		hap_df = [i.split() for i in hap_df]
		hap_df = pd.DataFrame(hap_df[1:],columns=hap_df[0])
		for col in ["H1","H2","HSR","N"]:
			hap_df[col] = hap_df[col].astype(float)
		return hap_df

	hap_data = PairData(relative_list,out)

	for chrm in range(1,num_chr+1):
		sys.stdout.write("\rCalculating hap score...chr %s" % chrm)
		sys.stdout.flush()
		genotype_data = GenotypeData(map_file,ped_file,chrm)

		def create_match_df(relative_list):
			match = pd.read_csv(match_file.replace("chr1","chr%s" % chrm),delim_whitespace=True,header=None)
			match["IID1"] = match[1].apply(lambda x: x[:-2])
			match["IID2"] = match[3].apply(lambda x: x[:-2])
			match["PAIR_ID"] = match[["IID1","IID2"]].min(axis=1) + "_" + match[["IID1","IID2"]].max(axis=1)
			match = match[match["PAIR_ID"].isin(relative_list)]
			match["HAP1"] = match[1].apply(lambda x: x[-1])
			match["HAP2"] = match[3].apply(lambda x: x[-1])
			match["CM_START"] = match[5].apply(genotype_data.mb_to_cm)
			match["CM_END"] = match[6].apply(genotype_data.mb_to_cm)
			match = match[["PAIR_ID","IID1","IID2","CM_START","CM_END","HAP1","HAP2"]].values.tolist()
			match.sort()
			return match

		match = create_match_df(hap_data.get_relatives())

		while match != []:
			pair_segs = list()
			pair,iid1,iid2 = match[0][:3]
			while pair in match[0]:
				pair_segs.append(match[0])
				match = match[1:]
				if match == []:
					break
			start1,end1,index_list,num_segs = pair_segs[0][3],pair_segs[0][4],[0],0
			for i in range(1,len(pair_segs)):
				start2,end2 = pair_segs[i][3:5]
				if start1 < start2 < end1 and start1 < end2 < end1: #Completely overlap
					continue
				index_list.append(i)
				if start1 == start2: #Segs have same start but seg2 is as long or longer
					del index_list[-2]
					continue
				if start2 > end1 + cm_gap or (start2 <= end1 + cm_gap and genotype_data.gap_discordance(iid1,iid2,end1,start2,disc_homoz)): #New seg
					num_segs += 1
					start1,end1 = start2,end2
					continue
				end1 = end2
			num_segs += 1
			hap_data.finish_chrm(pair,iid1,iid2,[pair_segs[i] for i in index_list],num_segs)            

	sys.stdout.write("\rCalculating hap score...done   \n")
	return hap_data.write_out(relative_list,out)

def start_up(parameter_dict,file_dict,run_type):
	def logo():
		tree = ["\n           P O N D E R O S A",
		"                v.1.5          ",
		"                  &",                   
		"                 ##&",                   
		"                #####",                  
		"              &%#######",              
		"            &%&%#####%###",              
		"          &   &%######&   &",            
		"             &%########&",               
		"            %%%%#########",              
		"          %%%%%###########&",           
		"        %%&  &%#########  &&&",          
		"      *     %%###########.",             
		"          &%%##############",            
		"        &&&%%%############&",         
		"     #  &%%%%##/(###########&/",         
		"      %%%%%%##%((%%############%&",      
		"   &%%%%&&  &%%%&&&%#######&&&&&&&&&",   
		"&&.&%     &%%&   &((&&  ###&&&&&",        
		"                 %((&       ",      
		"                 %*(&",                  
		"                 #((&",                 
		"                .(((%",                  
		"                &%(((%",                 
		"               %%#(#((((",               
		"              &&&&(  &(&)",
		"                          ",
		"               (c) 2020",
		"C.M. Williams, C.R. Gignoux, B.M. Henn\n"]

		def print_color(color,text):
			if color == "bold":
				color = "\033[1m"
			if color == "brown":
				color = "\033[0;33m"
			elif color == "green":
				color = "\033[0;32m"
			colored_text = f"\033[{color}{text}\033[00m"
			return colored_text

		for index,i in enumerate(tree):
				if index in [k for k in range(19,26)]:
					print(print_color("brown",i))
				elif index in [k for k in range(2,19)]:
					print(print_color("green",i))
				elif index == 0:
					print(print_color("bold",i))
				elif index in [1,27,28]:
					print(i)

	def write_args():
		sys.stdout.write("Run type: %s\n\n" % run_type)
		sys.stdout.write("Arguments used:\n")
		for files in list(dict.fromkeys(file_dict)):
			file_name = file_dict[files]
			if file_name == "None":
				file_name = "None provided"
			sys.stdout.write("    %s: %s\n" % (files,file_name))
		for parameters in list(dict.fromkeys(parameter_dict)):
			parameter_val = parameter_dict[parameters]
			sys.stdout.write("    %s: %s\n" % (parameters,parameter_val))

	def locate_files():
		#Check for files
		file_type = list(dict.fromkeys(file_dict))
		not_found = [file_dict[file] for file in file_type if (not os.path.isfile(file_dict[file]) and file_dict[file] != "None")]
		if not_found != []:
			sys.stdout.write("\nWARNING!\nThe following files were not found:\n")
			for msg in not_found:
				sys.stdout.write("\t%s\n" % msg)
			sys.stdout.write("\nExiting PONDEROSA...\n\n")
			sys.exit()

	logo()
	write_args()
	sys.stdout.write("\nChecking files...done\n")
	locate_files()
	return {**parameter_dict,**file_dict}

class LogFile:
	def __init__(self,parameters,run_type):
		self.start_time = time.time()
		self.logfile = open("%s.log" % parameters["out"],"w")
		self.logfile.write("PONDEROSA v1.5\n\n%s\n" % datetime.datetime.now())
		self.logfile.write("\nRun type: %s\n\nParameters/files provided:\n" % run_type)        
		for pars in list(dict.fromkeys(parameters)):
			self.logfile.write("\t%s: %s\n" % (pars,parameters[pars]))

	def write(self,line):
		self.logfile.write(line)

	def validation(self,second_df):
		pass

	def mz_twins(self,twins):
		if twins != {}:
			self.write("\nMZ twins present. PONDEROSA collapses MZ twins into one individual.\n")
			self.write("The ID on the left has been replaced by the ID on the right.\n")
			for iid in twins:
				self.write("%s %s\n" % (iid,twins[iid]))

	def write_log(self):
		time_elapsed = round(time.time() - self.start_time,1)
		self.write("\nTime elapsed: %s seconds\n\n" % time_elapsed)
		self.logfile.close()

	def write_errors(self,error_dict):
		if 2 in error_dict and error_dict[2] != []:
			self.write("\nNon-critical errors detected. Please double check.\n")
			for errors in error_dict[2]:
				self.write(errors[0])
				for pairs in errors[1]:
					self.write("\t"+" ".join(pairs) + "\n")
		for errors in error_dict[1]:
			self.write("Critical errors detected.\n")
			self.write(errors[0])
			self.write_log()
			sys.stdout.write("\nCritical error detected. See log file.\nExiting PONDEROSA...\n")
			sys.exit()

class Pedigree:
	def __init__(self,pedigree=False):
		if not pedigree:
			self.pedigree_structure = dict()
		#can initilize with existing pedigree structure
		else:
			self.pedigree_structure = pedigree

		#list of all ppl
		self.ind_list = list()

		#keep track of dummy ids
		self.dummy_id = 1

		'''This pedigree structure is a dict with each ind mapping to their parents and their children.
		Moving from a child to a parent is said to be +1 and parent down to child -1. Thus to get from you to your
		full sibling you go up one and down one, so the coord is considered +1,-1. To go to your cousin, you go up
		two to your shared grandparent (+1,+1) and then down to your aunt/uncle and down to the cousin (-1,-1)'''
		self.ped_coord = {"PO":[1],
						  "FS":[1,-1],
						  "AV":[1,1,-1],
						  "CO":[1,1,-1,-1]}

		self.rel_to_deg = { "PO":"PO","FS":"FS",
							"PHS":"2nd","MHS":"2nd","GP":"2nd","AV":"2nd",
							"CO":"3rd","GGP":"3rd","HAV":"3rd",
							"HCO":"4th","GGGP":"4th"}

		#Type 1 error is critical and exits PONDEROSA; type 2 is added to log file but does not stop running
		self.errors = {1:[],2:[]}

	def pair_ID(self,df,iids=["IID1","IID2"]):
		df["MAXID"] = df[iids].max(axis=1)
		df["MINID"] = df[iids].min(axis=1)
		df["PAIR_ID"] = df["MINID"] + "_" + df["MAXID"]

	def get_dummy_id(self):
		self.dummy_id += 1
		return "Missing" + "%03d" % (self.dummy_id - 1)

	def add_person(self,iid,sex):
		#only add if not already in
		if iid not in self.pedigree_structure:
			self.ind_list.append(iid)
			#each maps to [[dad],[mom],children,sex]
			self.pedigree_structure[iid] = [[],[],[],sex]

	def add_po(self,parent,child,sex,override=False):
		if self.pedigree_structure[child][sex-1] == [] or override: #parent not reported yet
			self.pedigree_structure[child][sex-1] += [(parent,1)]
			self.add_person(parent,sex)
			self.pedigree_structure[parent][2] += [(child,-1)]

	def get_coord(self,rel_type):
		return self.ped_coord[rel_type]

	def get_parent(self,iid,sex):
		parent = self.pedigree_structure[iid][sex-1]
		try:
			parent = [parent[0][0]]
		except:
			parent = []
		return parent

	def add_error(self,error_msg,pairs,error_type):
		if error_type == 1 or pairs != []:
			self.errors[error_type].append([error_msg,pairs])

	def add_mz_twins(self,twin_list):
		twin_list = [twins.split("_") for twins in twin_list]
		self.mz_twins = dict(twin_list)

	def get_mz_twins(self):
		return self.mz_twins

	def check_mz(self,iid_list):
		out_list = []
		for iid in iid_list:
			if iid in self.mz_twins:
				iid = self.mz_twins[iid]
			out_list.append(iid)
		return out_list

	def add_king(self,king_file):
		self.king = pd.read_csv(king_file,sep="\t")
		self.pair_ID(self.king,["ID1","ID2"])
		self.add_mz_twins(self.king[self.king["InfType"] == "Dup/MZ"]["PAIR_ID"].values.tolist())
		self.king["DUP"] = self.king["ID1"].isin(self.mz_twins) | self.king["ID2"].isin(self.mz_twins)
		self.king = self.king[~self.king["DUP"]]
		self.gtd = list(dict.fromkeys(self.king["ID1"].values.tolist() + self.king["ID1"].values.tolist()))
		self.king = self.king.iloc[:,[12,6,7,8,9]]
		self.king.columns = ["PAIR_ID","IBD1","IBD2","PIHAT","KINGINF"]

	def make_from_fam(self,fam):
		self.fam_file = open(fam).readlines()
		for ppl in self.fam_file:
			iid,dad,mom,sex = ppl.split()[1:5]
			iid,dad,mom = self.check_mz([iid,dad,mom])
			if iid in self.mz_twins:
				iid = self.mz_twins[iid]
			self.add_person(iid,sex)
			def add_fam_po(parent,child,sex):
				if parent != "0":
					self.add_po(parent,child,sex)
			add_fam_po(dad,iid,1)
			add_fam_po(mom,iid,2)
		self.ind_list.sort()

	def get_first(self,iid,direction):
		if direction == -1:
			return_rels = [rels for rels in self.pedigree_structure[iid][2]]
		if direction == 1:
			return_rels = [rels for rels in self.pedigree_structure[iid][0]] + [rels for rels in self.pedigree_structure[iid][1]]
		return return_rels

	def get_relationships(self,numeric_path,iid,stop_list):
		return_rels = [iid]
		for direction in numeric_path:
			rel_lists = [[rel2[0] for rel2 in self.get_first(rel1,direction)] for rel1 in return_rels]
			return_rels = []
			for rel_l in rel_lists:
				return_rels += [rel for rel in rel_l if (direction == 1 or rel not in stop_list)]
		return return_rels

	def get_lineal(self,stop=4):
		out_list = []
		for iid in self.ind_list:
			if stop == 4:
				out_list += [[iid,self.get_parent(iid,sex)[0],"PO"] for sex in [1,2] if self.get_parent(iid,sex) != []]
			for gen in range(1,stop):
				coord = [1 for i in range(1,gen+2)]
				new_gp = self.get_relationships(coord,iid,[])
				if new_gp == []:
					break
				rel = "G" * gen + "P"
				out_list += [[iid,gparent,rel] for gparent in new_gp]
		return out_list

	def parent_comp(self,iid1,iid2,sex,phase): #returns True if 
		#We know that iid1 and iid2 are connected by one or two parents
		#Let's say sex = 1 (male); if phase 1 we want to know if they have different mothers (bc if yes --> PHS)
		#If phase 2: we want to know if they have the same father (bc if yes --> PHS)
		if phase == 1:
			sex = {1:2,2:1}[sex]
		p1,p2 = self.get_parent(iid1,sex),self.get_parent(iid2,sex)
		if phase == 1: #want to know if diff parents
			status = (p1 != p2) and (p1 != [] or p2 != [])
		elif phase == 2: #want to know if same parents
			status = (p1 == p2) and (p1 != [] or p2 != [])
		return status

	#returns 2 list: one of pairs connected by 2 lineages (full) and one connected by one (half)
	def one_and_two(self,iid,rel_list,rel_type):
		dups_rm = list(dict.fromkeys(rel_list))
		one = [rel for rel in dups_rm if (rel_list.count(rel) == 1 and (rel_type == "AV" or iid < rel))]
		two = [rel for rel in dups_rm if (rel_list.count(rel) == 2 and (rel_type == "AV" or iid < rel))]
		return [one,two]

	#Phase 1: there are FS with missing parents --> sibs connected by only 1 lineage CAN BE full sibs
	#Phase 2: no FS with missing parents --> sibs connected by only 1 lineage ARE NOT full sibs
	def get_sibling_rels(self,rel_type,phase=2): #sibling rel = sib, av, co
		out_list = list()
		def analyze_rels(relative_list):
			single_lin,double_lin = self.one_and_two(iid,relative_list,rel_type)
			for relative in double_lin: #either FS,CO,or AV
				out_list.append([iid,relative,rel_type])
			for relative in single_lin:
				if rel_type == "FS":
					if self.parent_comp(iid,relative,2,phase): #MHS
							out_list.append([iid,relative,"MHS"])
							continue
					elif self.parent_comp(iid,relative,1,phase): #PHS
							out_list.append([iid,relative,"PHS"])
							continue
					if phase == 1: #If phase 1 and they haven't been added, they can still be FS
						out_list.append([iid,relative,"UNK"])
				else:
					out_list.append([iid,relative,"H"+rel_type])
		'''Have to treat diff relationships separately. For siblings, we want to return to the children of IIDs parent.
		For CO and AV, we dont want this, otherwise IIDs sibs will be cousins and IIDs parents will be uncles/aunts.
		For CO and AV we find the rels twice one on the maternal and once on the paternal lineage. This is to distinguish
		lineages in cases where an individual is related through more than one parent. For instance, if an ind is HC from mom
		and HC from dad, if we dont treat them as different we will think they are full cousins'''
		for iid in self.ind_list:
			if rel_type == "FS":
				analyze_rels(self.get_relationships(self.get_coord(rel_type),iid,[]))
			else:
				for sex in [1,2]:
					parent = self.get_parent(iid,sex)
					if parent != []:
						analyze_rels(self.get_relationships(self.get_coord(rel_type)[1:],parent[0],parent))
		return out_list

	#Takes in a list of parents (e.g. all the mothers of a set of FS) and looks for 1) multiple mothers and 2) no mother
	def resolve_parents(self,parent_list):
		parent_list = [i[0] for i in parent_list if i != []]
		parent_list = list(dict.fromkeys(parent_list))
		if len(parent_list) == 0: #No mom/dad --> create dummy parent
			parent = [self.get_dummy_id(),[]]
		elif len(parent_list) == 1: #Exactly one mom/dad --> a-OK
			parent = [parent_list[0],[]]
		elif len(parent_list) > 1: #More than one mom/dad --> will raise error
			parent = [parent_list[0],parent_list]
		#returns list with parent 0th index and duplicate parents in 1st index (if any)
		return parent

	def resolve_siblings(self,trust_fs):
		'''In phase 1, we assume that there are FS with 1 (or more) missing parents. Thus,
		absence of a lineage does not necessarily mean sibs are HS. This function finds FS with BOTH lineages
		and sibs with one lineage and where the other lineage does NOT connect so are HS. It then
		resolves FS by making sure FS pairs 1) have both parents and 2) have the same set of parents. This function
		runs during phase 1 to solve single lineage siblings (e.g. place them in FS or HS category) and prepares
		the software for phase 2, where it is assumed that all single lineage siblings are HS.'''

		'''king_fs TRUE if there are no FS in dataset with both parents; used when pedigree structure is sparse.
		king_fs list contains all pairs inferred by KING as FS'''
		if trust_fs == "True":
			king_fs = self.king[(self.king["KINGINF"] == "FS") & (self.king["IBD2"] > 0.15)]["PAIR_ID"].values.tolist()
			fs_pairs = [[sibs.split("_")[0],sibs.split("_")[1],"FS"] for sibs in king_fs]
			low_IBD2 = self.king[(self.king["KINGINF"] == "FS") & (self.king["IBD2"] <= 0.15)]["PAIR_ID"].values.tolist()
			self.add_error("ERROR05: KING has inferred the following as FS, but have low IBD2 values (< 0.15):\n",[pairs.split("_") for pairs in low_IBD2],2)
		else:
			king_fs = self.king[self.king["KINGINF"] == "FS"]["PAIR_ID"].values.tolist()
			king_fs = [sibs.split("_")+["UNK"] for sibs in king_fs]

			#all siblings from the pedigree structure
			av_pairs = [rels for rels in self.get_sibling_rels("AV") if rels[2] == "AV"]
			all_siblings = self.get_sibling_rels("FS",1) + av_pairs + self.get_lineal(2)
			all_siblings = pd.DataFrame(all_siblings,columns=["IID1","IID2","TYPE"])
			self.pair_ID(all_siblings)

			#all full siblings as inferred by king
			king_fs = pd.DataFrame(king_fs,columns=["IID1","IID2","TYPE"])
			self.pair_ID(king_fs)

			#merge the two dataframes and drop king pairs that are already known
			all_siblings = pd.concat([all_siblings,king_fs]).drop_duplicates("PAIR_ID").reset_index(drop=True)

			#add IBD values to the dataframe
			all_siblings["NUMERIC_TYPE"] = all_siblings["TYPE"].replace({"FS":0,"MHS":1,"PHS":1,"GP":1,"AV":1,"UNK":2})
			all_siblings = pd.merge(all_siblings,self.king,how="left",on=["PAIR_ID","PAIR_ID"])
			all_siblings = all_siblings.dropna()

			#run unresolved pairs thru the LDA
			unresolved_pairs = all_siblings[all_siblings["NUMERIC_TYPE"] == 2][["IID1","IID2"]].values.tolist()
			unresolved_ibd = all_siblings[all_siblings["NUMERIC_TYPE"] == 2][["IBD1","IBD2"]].values.tolist()

			fs_pairs = all_siblings[all_siblings["NUMERIC_TYPE"] == 0][["IID1","IID2"]].values.tolist()
			if len(unresolved_pairs) > 0:
				#create training arrays and LDA
				training = all_siblings[all_siblings["NUMERIC_TYPE"] != 2]
				train_labs,train_vals = training["NUMERIC_TYPE"].values.tolist(),training[["IBD1","IBD2"]].values.tolist()
				if train_labs.count(0) == 0 or train_labs.count(1) == 0:
					self.add_error("ERROR04: Sparse pedigree warning: not enough training pairs. Try rerunning with king_fs as True\n",[],1)
					return

				classif = LinearDiscriminantAnalysis().fit(train_vals,train_labs)

				#predict whether each sibship is FS or HS
				predicted_rels = classif.predict(unresolved_ibd)
				predicted_rels = [[unresolved_pairs[i][0],unresolved_pairs[i][1],predicted_rels[i]] for i in range(len(unresolved_pairs))]
				fs_pairs += [[sibs[0],sibs[1]] for sibs in predicted_rels if sibs[2] == 0]
			
		'''Here creating a list of all full sib sets (i.e. where every set is FS) to make sure they have parents/
		have same parents'''
		fs_sets = list()
		for pairs in fs_pairs:
			iid1,iid2 = pairs[:2]
			new_set = True
			for sets in fs_sets:
				if iid1 in sets and iid2 not in sets:
					new_set = False
					sets.append(iid2)
					break
				if iid2 in sets and iid1 not in sets:
					new_set = False
					sets.append(iid1)
					break
			if new_set:
				fs_sets.append([iid1,iid2])

		error_pairs = []
		for sets in fs_sets:
			dad = self.resolve_parents([self.get_parent(sibs,1) for sibs in sets])
			mom = self.resolve_parents([self.get_parent(sibs,2) for sibs in sets])
			if dad[1] != [] or mom[1] != []:
				error_pairs.append(sets)
				continue
			dad,mom = dad[0],mom[0]
			for sibs in sets:
				self.add_po(dad,sibs,1)
				self.add_po(mom,sibs,2)
		self.add_error("ERROR03: The following FS sets have different parents. The problem has been ignored but please double check.\n",error_pairs,2)

	#Creates a pandas table of all pedigree relationships present in the structure at the time of running it
	def get_all_pairs(self):
		rel_list = []
		for rels in ["FS","AV","CO"]:
			rel_list += self.get_sibling_rels(rels)
		rel_list += self.get_lineal()
		self.relatives = pd.DataFrame(rel_list,columns=["IID1","IID2","REL"])
		self.pair_ID(self.relatives)
		self.relatives = self.relatives[["PAIR_ID","IID1","IID2","REL"]]
		self.relatives = pd.merge(self.relatives,self.king,how="left",on=["PAIR_ID","PAIR_ID"])
		self.relatives["GTD"] = self.relatives["PIHAT"].notnull()
		self.relatives["DEGREE"] = self.relatives["REL"].replace(self.rel_to_deg)
		self.king["FOUND"] = self.king["PAIR_ID"].isin(self.relatives["PAIR_ID"])

	def get_relatives(self,rel,gt_only=True):
		if gt_only:
			out_list = len(self.relatives[(self.relatives["REL"]==rel) & (self.relatives["GTD"] == True)][["IID1","IID2"]].values.tolist())
		else:
			out_list = len(self.relatives[self.relatives["REL"]==rel][["IID1","IID2"]].values.tolist())
		return out_list

	def check_po(self):
		king_po = self.king[self.king["KINGINF"] == "PO"]["PAIR_ID"].values.tolist()
		fam_po = self.relatives[(self.relatives["REL"] == "PO") & (self.relatives["GTD"] == True)]["PAIR_ID"].values.tolist()
		missing = [[pairs.split("_")[0],pairs.split("_")[1]] for pairs in king_po if pairs not in fam_po]
		self.add_error("ERROR02: The following PO pairs have been inferred by KING but are not reported in the .fam file.\n",missing,2)

	def get_king(self):
		return self.king

	def get_rels(self):
		return self.relatives

	def get_pedigree_structure(self):
		return self.pedigree_structure

	def ruleout_hs(self,iid1,iid2,sex):
		parent1,parent2 = self.get_parent(iid1,sex),self.get_parent(iid2,sex)
		for parents in [parent1,parent2]:
			if parents != [] and parents[0] in self.gtd:
				return True
		return False

	def print_out(self,out):
		with open("%s_pairs.txt" % out,"w") as outfile:
			outfile.write(self.relatives.to_string(index=False,na_rep="NA",columns=["PAIR_ID","IID1","IID2","GTD","IBD1","IBD2","PIHAT","KINGINF","REL","DEGREE"]))

	def run_PONDEROSA(self,king_file,fam_file,out,trust_fs):
		self.add_king(king_file)
		self.make_from_fam(fam_file)
		self.resolve_siblings(trust_fs)
		self.get_all_pairs()
		self.print_out(out)
		self.check_po()
		return self.errors

def main():
	def init(par_file):
		par_file = [lines.strip() for lines in open(par_file).readlines()]
		run_type = [lines.split()[0] for lines in par_file[1:4] if "True" in lines][0]
		parameters = {lines.split()[0]:lines.split()[1] for lines in par_file[15:]}
		files = {lines.split()[0]:lines.split()[1] for lines in par_file[6:13]}
		pars = start_up(parameters,files,run_type)
		return pars,run_type

	def run_hapscores(kingf,hap_file):
		if [lines for lines in open(kingf).readlines()[1:] if float(lines.split()[6]) > 0.90] == []:
			log.write_errors({1:[["No KING-defined PO pairs."]]})
		relative_list = [[lines.split()[1],lines.split()[3]] for lines in open(kingf).readlines()[1:] if float(lines.split()[6]) > 0.30]
		relative_list = [min(pair) + "_" + max(pair) for pair in relative_list]
		return get_hap_score(relative_list,pars,hap_file)

	#Input if a df with IID1,IID2 and hap scores calculated; output is the same df with 4 new fields:
	#AGE1,AGE2,younger ind,older ind. If age data is available for both, younger/older is determined by age.
	#Else: determined by individual haplotype scores
	def resolve_generations(df,agef,young_lab,old_lab,out):
		def add_age(df):
			age_data = {}
			if agef != "None":
				age_data = {lines.split()[0]:float(lines.split()[1]) for lines in open(agef).readlines()}
			df["AGE1"],df["AGE2"] = df["IID1"].map(age_data),df["IID2"].map(age_data)
			df["USE_H"] = df["AGE1"].isna() | df["AGE2"].isna()
			df["USE_H"] = np.where(df["AGE1"]==df["AGE2"],True,df["USE_H"])

		def resolve_h(df):
			hap_data = df[["IID1","H1","IID2","H2"]].values.tolist()
			df[old_lab] = [{i[1]:i[0],i[3]:i[2]}[min(i[1],i[3])] for i in hap_data]
			df[young_lab] = [{i[1]:i[0],i[3]:i[2]}[max(i[1],i[3])] for i in hap_data]
			df["METHOD"] = "H"
			df["STRENGTH"] = abs(df["H1"] - df["H2"])
			return df
		def resolve_age(df):
			age_data = df[["IID1","AGE1","IID2","AGE2"]].values.tolist()
			df[old_lab] = [{i[1]:i[0],i[3]:i[2]}[max(i[1],i[3])] for i in age_data]
			df[young_lab] = [{i[1]:i[0],i[3]:i[2]}[min(i[1],i[3])] for i in age_data]
			df["METHOD"] = "AGE"
			df["STRENGTH"] = abs(df["AGE1"] - df["AGE2"])
			return df

		add_age(df)
		w_age = resolve_age(df.copy()[~df["USE_H"]])
		wo_age = resolve_h(df.copy()[df["USE_H"]])
		df = pd.concat([w_age,wo_age]).drop("USE_H",axis=1)
		return df

	def po_analysis(hap_df,agef,out):
		po = [[lines.split()[1],lines.split()[3]] for lines in open(pars["king_file"]).readlines()[1:] if lines.split()[-1] == "PO"]
		po_data = pd.DataFrame([min(pairs) + "_" + max(pairs) for pairs in po],columns=["PAIR_ID"])
		po_data = pd.merge(po_data,hap_df,on="PAIR_ID",how="left")
		po_data = resolve_generations(po_data,agef,"CHILD","PARENT",out)
		with open("%s_PO.txt" % out,"w") as outfile:
			  outfile.write(po_data.to_string(index=False,columns=["PAIR_ID","IID1","IID2","H1","H2","AGE1","AGE2","CHILD","PARENT","METHOD","STRENGTH"]))

	def phase3_checkpoint(out):
		pairs_df = pd.read_csv("%s_pairs.txt" % out, delim_whitespace=True)
		pairs_df = pairs_df[pairs_df["GTD"]]
		num_rels = [pairs_df[pairs_df["REL"]==rel].shape[0]>0 for rel in ["MHS","PHS","AV","GP"]] + [pairs_df[pairs_df["DEGREE"]=="3rd"].shape[0] > 0]
		if sum(num_rels) == 5:
			return
		problem_rels = [rel for rel,status in zip(["MHS","PHS","AV","GP","3rd"],num_rels) if not status]
		error_msg = "ERROR01: Not enough of the following pairs found: " + " ".join(problem_rels) + "\n"
		log.write_errors({1:[[error_msg]]})

	def infer_second(king_df,relative_df,hap_df,threshold,mhs_gap,gp_gap):
		class Data:
			def __init__(self,king,rel,hap_scores,threshold,mhs_gap,gp_gap):
				#remove outliers; main goal is to remove duplicate ind (who have >1 rel)
				degrees = ["FS","2nd","3rd"]
				self.second = ["AV","GP","MHS","PHS"]
				def remove_outliers(df):
					df = df[df["GTD"] & df["DEGREE"].isin(degrees)]
					third_df = df[df["DEGREE"]=="3rd"].copy()
					df = df[df["DEGREE"] != "3rd"]
					third_mean = stat.mean(third_df["IBD1"])
					third_sd = stat.stdev(third_df["IBD1"])
					third_df = third_df[third_df.apply(lambda x: abs((x.IBD1-third_mean)/third_sd) < 2,axis=1)]
					third_df = third_df[~third_df["PAIR_ID"].isin(df["PAIR_ID"])]
					df = pd.concat([df,third_df])
					df = df.drop(["IID1","IID2","IBD1","IBD2","PIHAT","KINGINF","GTD"],axis=1)
					return df

				self.putative = king[(king["IBD1"] < 0.75) & (king["IBD1"] > 0.30)]
				self.putative = pd.merge(self.putative,hap_scores,on="PAIR_ID",how="left").dropna()
				rel = remove_outliers(rel)
				self.putative = pd.merge(self.putative,rel,on="PAIR_ID",how="left")

				self.training = self.putative[self.putative["DEGREE"].isin(degrees)].copy()
				self.putative = self.putative[self.putative["DEGREE"].isna()]
				self.gp_gap = gp_gap
				self.mhs_gap = mhs_gap

			def check_error(self,train_lab,lab_types,error_msg):
				for labs in lab_types:
					if train_lab.count(labs) == 0:
						log.write_errors({1:[["ERROR01: Not enough %s pairs to train %s classifier" % (labs,error_msg),[]]]})

			def get_training(self,df,lab_type,val_list):
				vals,labs = df[val_list].values.tolist(),df[lab_type].values.tolist()
				return vals,labs

			def find_putative(self):
				train_val,train_lab = self.get_training(self.training,"DEGREE",["IBD1","IBD2"])
				self.check_error(train_lab,["2nd","3rd"],"degree")
				classif = LinearDiscriminantAnalysis().fit(train_val,train_lab)
				self.putative["SECOND_PROB"] = self.putative.apply(lambda x: classif.predict_proba([[x.IBD1,x.IBD2]])[0][0],axis=1)
				self.putative = self.putative[self.putative["SECOND_PROB"] > threshold]

			def classify_second(self,train_df,put_df):
				train_val,train_lab = self.get_training(train_df[train_df["DEGREE"]=="2nd"],"REL",["HSR","N"])
				self.check_error(train_lab,["AV","MHS","PHS","GP"],"2nd degree")
				classif = LinearDiscriminantAnalysis().fit(train_val,train_lab)
				probs = classif.predict_proba(put_df[["HSR","N"]].values.tolist())
				for index,rel in enumerate(self.second):
					put_df[rel] = [p[index] for p in probs]
				return put_df

			def set_zero(self,col,df): #where condition is true, sets the probability to 0
				df[col] = np.where(df["CONDITION"],0,df[col])

			def set_conditional(self,df): #recomputes conditional prob
				df["SUM"] = df.loc[:,self.second].sum(axis=1)
				for rels in self.second:
					df[rels] = df[rels]/df["SUM"]

			def ages(self,df): #constrain rels with age info
				df["CONDITION"] = abs(df["AGE1"]-df["AGE2"]) > self.mhs_gap
				self.set_zero("MHS",df)
				df["CONDITION"] = abs(df["AGE1"]-df["AGE2"]) < self.gp_gap
				self.set_zero("GP",df)

			def parent(self,df): #constrain HS rels with parent data
				for i,rel in enumerate(["PHS","MHS"]):
					df["CONDITION"] = df.apply(lambda x: pedigree.ruleout_hs(x.IID1,x.IID2,i+1),axis=1)
					self.set_zero(rel,df)

			def av_error(self,h1,age1,h2,age2):
				return {age1:h1,age2:h2}[max(age1,age2)] > {age1:h1,age2:h2}[min(age1,age2)]

			def write_out(self,out,df):
				self.set_conditional(df)
				df["REL"] = df[self.second].idxmax(axis=1)
				df["PROB"] = df[self.second].max(axis=1)
				df["AV_ERROR"] = np.where(df["REL"] == "AV",df.apply(lambda x: self.av_error(x.H1,x.AGE1,x.H2,x.AGE2),axis=1),False)
				with open("%s_second.txt" % out,"w") as outfile:
					outfile.write(df.to_string(index=False,na_rep="NA",columns=["PAIR_ID","YOUNGER","OLDER","METHOD","REL","SECOND_PROB","PROB","HSR","N","AV","GP","MHS","PHS","AV_ERROR"]))

			def run(self):
				self.find_putative()
				self.putative = self.classify_second(self.training,self.putative)
				self.putative = resolve_generations(self.putative,pars["age_file"],"YOUNGER","OLDER",pars["out"])
				self.parent(self.putative)
				self.ages(self.putative)
				self.set_conditional(self.putative)
				self.write_out(pars["out"],self.putative)

			def validation(self,out):
				second_test = self.training[self.training["DEGREE"]=="2nd"].copy()
				second_test = self.classify_second(self.training,second_test)
				second_test["PREDICTED"] = second_test[self.second].idxmax(axis=1)
				with open("%s_training_data.txt" % out,"w") as outfile:
					outfile.write(self.training.to_string(index=False,na_rep="NA"))
				log.validation(second_test)

		data = Data(king_df,relative_df,hap_df,threshold,gp_gap,mhs_gap)
		data.run()

	#Step 1: check files
	pars,run_type = init(sys.argv[-1])

	#Init log file
	log = LogFile(pars,run_type)

	#Skip hap score computation if ped only
	if run_type in ["po_only","run_all"]:
		#Step 2: compute hap scores
		hap_df = run_hapscores(pars["king_file"],pars["hap_file"])

		#Step 3: analyze PO pairs
		sys.stdout.write("Analyzing PO pairs...")
		po_analysis(hap_df,pars["age_file"],pars["out"])
		sys.stdout.write("\rAnalyzing PO pairs...done\n")

	#Quit program if po only
	if run_type == "po_only":
		log.write_log()
		sys.exit()

	#Step 4: create ped structure; resolve amb sibships, add missing parents, creates 2 df: king_df, rel_df
	sys.stdout.write("Building pedigree graphs...")
	pedigree = Pedigree()
	errors = pedigree.run_PONDEROSA(pars["king_file"],pars["fam_file"],pars["out"],pars["trust_fs"])
	log.write_errors(errors)
	relative_df = pedigree.get_rels()
	king_df = pedigree.get_king()
	log.mz_twins(pedigree.get_mz_twins())
	sys.stdout.write("\rBuilding pedigree graphs...done\n")

	#Quit program if ped only
	if run_type == "ped_only":
		log.write_log()
		sys.exit()
	else:
		phase3_checkpoint(pars["out"])

	#Step 5: infer second deg pairs
	sys.stdout.write("Finding and inferring 2nd degree pairs...")
	infer_second(king_df,relative_df,hap_df,float(pars["likelihood"]),int(pars["gp_gap"]),int(pars["mhs_gap"]))
	sys.stdout.write("\rFinding and inferring 2nd degree pairs...done\n")

	#Finish
	sys.stdout.write("Writing log...")
	log.write_log()
	sys.stdout.write("\rWriting log...done\n")

main()
