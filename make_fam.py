import sys
import pandas as pd
import numpy as np

def make_fam(fam,po_output,king,out,po_gap):
	king = pd.read_csv(king,delim_whitespace=True)
	king["PAIR_ID"] = king.apply(lambda x: min([x.ID1,x.ID2]) + "_" + max([x.ID1,x.ID2]),axis=1)
	king.set_index("PAIR_ID",inplace=True)

	try:
		po = pd.read_csv(po_output,delim_whitespace=True)
	except:
		print("PO output from PONDEROSA not found")
		sys.exit()
	if po.shape[0] == 0:
		print("PO output provided but no pairs present")
		sys.exit()
	po["INDEX"] = po["PAIR_ID"]
	po.set_index("INDEX",inplace=True)
	po["MIN_H"] = po.apply(lambda x: min([x.H1,x.H2]),axis=1)
	po["LOCKED"] = po["STRENGTH"] > 15.0
	po["PARENT_AGE"] = po.apply(lambda x: {x.IID1:x.AGE1,x.IID2:x.AGE2}[x.PARENT],axis=1)

	fam = pd.read_csv(fam,delim_whitespace=True,names=["FID","IID","PID","MID","SEX","PHE"])
	fam["INDEX"] = fam["IID"]
	fam.set_index("INDEX",inplace=True)
	fam["PID"] = fam["PID"].astype(str)
	fam["MID"] = fam["MID"].astype(str)
	fam["SEX"] = fam["SEX"].astype(int)

	problemParents = {"sex":{},"dup":{}}

	logFile = {}
	def writeLog(subj,msg):
		if subj == "write":
			out = open("%s_fam.log" % msg,"w")
			errors = {"multipleSex":"Each group (1 & 2) should be different genders, but provided sex assignments do not match.\n",
					  "duplicate":"The following individuals each have 2+ parents of the same gender.\n",
					  "ageGap":"The following PO pairs have small age gaps.\n",
					  "relatedParents":"The following individual has parents who are at least 3rd degree related according to KING.\n"}
			for subj in logFile:
				out.write(errors[subj])
				for msg in logFile[subj]:
					out.write(msg)
		if subj not in logFile:
			logFile[subj] = []
		logFile[subj].append(msg)

	def getSex(iid):
		return fam.at[iid,"SEX"]

	def writeParent(child,parent,df,problemParents):
		def addProblem(child,parent,problem):
			if child not in problemParents[problem]:
				problemParents[problem][child] = []
			problemParents[problem][child].append(parent)
		sex = getSex(parent)
		if sex == 0:
			addProblem(child,parent,"sex")
		else:
			if df.at[child,{1:"PID",2:"MID"}[sex]] == "0":
				df.at[child,{1:"PID",2:"MID"}[sex]] = parent
			else:
				cur_parent = df.at[child,{1:"PID",2:"MID"}[sex]] 
				if cur_parent != parent:
					addProblem(child,parent,"dup")

	def getPairID(iid1,iid2):
		return min([iid1,iid2]) + "_" + max([iid1,iid2])

	def switchPO(child,parent,lock=False):
		pair_id = getPairID(child,parent)
		if not po.at[pair_id,"LOCKED"]:
			po.at[pair_id,"CHILD"],po.at[pair_id,"PARENT"] = child,parent
			if lock:
				po.at[pair_id,"LOCKED"] = True

	def writeFam(df):
		problemParents = {"sex":{},"dup":{}}
		po.apply(lambda x: writeParent(x.CHILD,x.PARENT,df.copy(),problemParents),axis=1)
		return problemParents

	def getIBD2():
		def two_parents(iid):
			return po[(po["CHILD"]==iid)&(po["LOCKED"])].shape[0] == 2
		fs_pairs = king[king["InfType"]=="FS"][["ID1","ID2","IBD2Seg"]].values.tolist()
		ibd2_val = []
		for pairs in fs_pairs:
			iid1,iid2,ibd2 = pairs
			if two_parents(iid1) and two_parents(iid2):
				ibd2_val.append(ibd2)
		if ibd2_val == []:
			ibd2_val = [0.18]
		return np.mean(ibd2_val)-(1.5*np.std(ibd2_val))

	def getSets():
		ibd2 = getIBD2()
		sets,map_sets,new_index = [[None]],{},1
		def getIndex(iid):
			if iid not in map_sets:
				map_sets[iid] = 0
			return map_sets[iid]
		for pairs in king[king["IBD2Seg"] > ibd2][["ID1","ID2"]].values.tolist():
			iid1,iid2 = pairs
			index = max([getIndex(iid1),getIndex(iid2)])
			if index == 0:
				sets.append([iid1,iid2])
				map_sets[iid1],map_sets[iid2] = new_index,new_index
				new_index += 1
			else:
				map_sets[iid1],map_sets[iid2] = index,index
				sets[index] = list(set(sets[index] + [iid1,iid2]))
		return sets[1:]

	def getPO(iid,rtype):
		lines = po[(po["IID1"]==iid) | (po["IID2"]==iid)][["CHILD","PARENT"]].values.tolist()
		children,parents = [i[0] for i in lines if i[0] != iid],[i[1] for i in lines if i[1] != iid]
		return {"all":children+parents,"children":children,"parents":parents}[rtype]

	def addFS():
		ibd2 = getIBD2()
		for sets in getSets():
			all_plist = [getPO(iid,"all") for iid in sets]
			parents = list(set(all_plist[0]) & set(all_plist[1]))
			for iid,iid_plist in zip(sets,all_plist):
				children = [i for i in iid_plist if i not in parents]
				for parent in parents:
					switchPO(iid,parent,True)
				for child in children:
					switchPO(child,iid,True)

	def KINGdegree(iid1,iid2):
		pair_id = getPairID(iid1,iid2)
		try:
			return king.at[pair_id,"InfType"]
		except:
			return "UN"

	def youngParent():
		for pair in po[po["PARENT_AGE"]<po_gap].values.tolist():
			parent,child = pair[7:9]
			switchPO(child,parent,True)

	def addUnrelated():
		def findUnrelated(iid):
			plist = getPO(iid,"all")
			return_list = [plist[:]]
			for iid1 in plist:
				plist = plist[1:]
				for iid2 in plist:
					if KINGdegree(iid1,iid2) in ["4th","UN"]:
						return return_list + [iid1,iid2]
			return []
		for iid in fam["IID"]:
			unrelatedPair = findUnrelated(iid)
			if unrelatedPair != []:
				parent1,parent2 = unrelatedPair[1:]
				plist = [i for i in unrelatedPair[0] if i not in [parent1,parent2]]
				switchPO(iid,parent1,True)
				switchPO(iid,parent2,True)
				switch = [switchPO(i,iid,True) for i in plist]

	def estSex(l1,l2,sex_dict):
		def meanMin(l):
			return np.mean(po[po["PARENT"].isin(l)]["MIN_H"])
		sex1,sex2 = np.unique([getSex(i) for i in l1] + [0])[1:], np.unique([getSex(i) for i in l2] + [0])[1:]
		if len(sex1) == 1 and len(sex2) == 1 and sex1 != sex2:
			sex_dict[sex1[0]] += l1
			sex_dict[sex2[0]] += l2
			return sex_dict
		if len(sex1) > 1 or len(sex2) > 1:
			line1 = "\t".join(["Group 1:"] + ["%s (%s)" % (i,getSex(i)) for i in l1 if getSex(i) != 0] + ["\n"])
			line2 = "\t".join(["Group 2:"] + ["%s (%s)" % (i,getSex(i)) for i in l2 if getSex(i) != 0] + ["\n"])
			writeLog("multipleSex",line1+line2)
		l = [[meanMin(l1)]+l1,[meanMin(l2)]+l2]
		l.sort()
		sex_dict[1] += l[1][1:]
		sex_dict[2] += l[0][1:]
		return sex_dict

	def orientSex(parentSets):
		sex = {0:[],1:[]}
		def makeDict():
			return_dict = {j:[] for j in [i[0] for i in parentSets]+[i[1] for i in parentSets]}
			for sets in parentSets:
				iid1,iid2 = sets
				return_dict[iid1].append(iid2)
				return_dict[iid2].append(iid1)
			return {i:list(np.unique(return_dict[i])) for i in return_dict}
		def added(iid,sex_dict):
			return iid in sex_dict[1] or iid in sex_dict[2]
		sex_dict = {1:[],2:[]}
		parentSets = makeDict()
		def iterSet(iid1,num):
			counter = parentSets[iid1]
			if iid1 in sex[num]:
				return
			sex[num].append(iid1)
			if counter == []:
				return
			else:
				num += 1
				for iid2 in counter:
					if iid2 not in sex[num%2]:
						iterSet(iid2,num%2)
				return
		for iid in parentSets:
			if added(iid,sex_dict):
				continue
			iterSet(iid,0)
			sex_dict = estSex(sex[0],sex[1],sex_dict)
			sex = {0:[],1:[]}
		return sex_dict

	def noSex(problemParents):
		sets = [problemParents[i] for i in problemParents]
		singleSets = [i[0] for i in sets if len(i) == 1]
		doubleSets = [i for i in sets if len(i) == 2]
		sex_dict = orientSex(doubleSets)
		male,female = sex_dict[1],sex_dict[2]
		def addSingle(iid,fMean,mMean):
			if fMean == 0 and mMean == 0:
				male.append(iid)
			else:
				iidMean = np.mean(po[po["PARENT"]==iid]["MIN_H"])
				if abs(iidMean-fMean) < abs(iidMean-mMean):
					female.append(iid)
				else:
					male.append(iid)
		fMean,mMean = 0,0
		if male != [] and female != []:
			fMean = np.mean(po[po["PARENT"].isin(female)]["MIN_H"])
			mMean = np.mean(po[po["PARENT"].isin(male)]["MIN_H"])
		for iid in singleSets:
			if iid in male or iid in female:
				continue
			addSingle(iid,fMean,mMean)
		return [male,female]

	def writeDuplicates(iid,parent1):
		sex = getSex(parent1[0])
		parent2 = fam.at[iid,{1:"PID",2:"MID"}[sex]][0]
		writeLog("duplicate","\t".join(["%s: " % iid] + ["%s" % i for i in parent2] + ["\n"]))

	def addSex(l,sex):
		for iid in l:
			fam.at[iid,"SEX"] = sex

	def ageIssues():
		for pairs in po[(po["METHOD"]=="AGE") & (po["STRENGTH"]<=15.0)][["IID1","IID2","AGE1","AGE2"]].values.tolist():
			writeLog("ageGap","%s %s %s %s\n" % tuple(pairs))

	def relatedParents(df):
		for line in df[(df["PID"] != "0")&(df["MID"] != "0")][["IID","PID","MID"]].values.tolist():
			iid,pid,mid = line
			deg = KINGdegree(pid,mid)
			if deg in ["2nd","FS","PO"]:
				writeLog("relatedParents","%s: %s %s (%s)\n" % (iid,pid,mid,deg))

	def run():
		youngParent()
		addUnrelated()
		addFS()
		problemParents = writeFam(fam)
		if len(problemParents["sex"]) > 0: #if there are any PO pairs with no sex
			fam_noSex = fam.copy()
			fam_noSex["SEX"] = [0 for _ in range(fam.shape[0])]
			problemParents = writeFam(fam_noSex)
			pid,mid = noSex(problemParents["sex"])
			addSex(pid,1)
			addSex(mid,2)
		problemParents = {"sex":{},"dup":{}}
		po.apply(lambda x: writeParent(x.CHILD,x.PARENT,fam,problemParents),axis=1)
		duplicateParents = problemParents["dup"]
		for iid in duplicateParents:
			writeDuplicates(iid,duplicateParents[iid])
		ageIssues()
		relatedParents(fam)
		writeLog("write",out)
		outFam = open("%s_PONDEROSA.fam" % out,"w")
		lines = fam.values.tolist()
		for i in lines:
			i = [str(j) for j in i]
			outFam.write(" ".join(i) + "\n")

	run()

par_file = sys.argv[-1]
par_file = pd.read_csv(par_file,delim_whitespace=True)
params = {"fam_file":None,"out":None,"king_file":None,"po_gap":None}
print("Running with the following input:")
for i in params:
	params[i] = par_file.at[i,"Run_type"]
	print(params[i])
params["po_output"] = "%s_PO.txt" % params["out"]
print(params["po_output"])

make_fam(params["fam_file"],params["po_output"],params["king_file"],params["out"],float(params["po_gap"]))
