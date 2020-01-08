from optparse import OptionParser
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import random
import statistics as st

parser = OptionParser()

parser.add_option("--pedigree",dest="pedigree")
parser.add_option("--king",dest="king")
parser.add_option("--ersa",dest="ersa")
parser.add_option("--training",dest="training")

(option,args) = parser.parse_args()

def pair_id(iid1,iid2):
	temp = [iid1,iid2]
	temp.sort()
	iid1,iid2 = temp[0],temp[1]
	iid1 = iid1.split("HMB")[1]
	iid2 = iid2.split("HMB")[1]
	if iid1 == "678b":
		iid1 = "679"
	if iid2 == "678b":
		iid2 = "679"
	return "%04d" % int(iid1) + "%04d" % int(iid2)

#Convert the ERSA file
out = open("new" + option.ersa, "w")
out.write("ID ERSAdeg ERSArel\n")
for pairs in open(option.ersa).readlines()[3:]:
	pairs = pairs.split()
	iid1,iid2,ancestors,degree = pairs[0],pairs[1],pairs[2],pairs[3]
	unique_id = pair_id(iid1,iid2)
	if degree == "1" and ancestors == "0":
		degree,rel = "PO","PO"
	elif degree == "1" and ancestors == "2":
		degree,rel = "FS","FS"
	elif degree == "2":
		degree = "2nd"
		if ancestors == "0":
			rel = "GP"
		if ancestors == "1":
			rel = "HS"
		if ancestors == "2":
			rel = "AV"
	elif degree == "3":
		degree = "3rd"
		if ancestors == "0":
			rel = "GGP"
		if ancestors == "1":
			rel = "HV"
		if ancestors == "2":
			rel = "CO"
	elif degree == "4":
		degree,rel = "4th","NA"
	elif degree not in ["PO","FS","2nd","3rd","4th"]:
		degree,rel = "UN","NA"
	out.write(" ".join([unique_id,degree,rel,"\n"]))
out.close()

#Convert the KING file
out = open("new" + option.king, "w")
out.write("ID IBD1 IBD2 KINGdeg\n")
for pairs in open(option.king).readlines()[1:]:
	pairs = pairs.split()
	iid1,iid2,ibd1,ibd2,rel = pairs[1],pairs[3],pairs[6],pairs[7],pairs[9]
	unique_id = pair_id(iid1,iid2)
	out.write(" ".join([unique_id,ibd1,ibd2,rel,"\n"]))
out.close()

#Convert the PONDEROSA outputs
rel_to_deg = {"PO":"PO","FS":"FS",
			  "PHS":"2nd","MHS":"2nd","GP":"2nd","AV":"2nd",
			  "GGP":"3rd","CO":"3rd","HV":"3rd",
			  "GGGP":"4th","HC":"4th"}
pair_list = {"2nd":[],"3rd":[],"4th":[]}
for pairs in open(option.pedigree).readlines()[1:]:
	pairs = pairs.split()
	iid1,iid2,rel,confidence = pairs[0],pairs[1],pairs[2],pairs[3]
	unique_id = pair_id(iid1,iid2)
	if confidence != "1.0" or rel == "FS" or rel == "PO":
		continue
	pair_list[rel_to_deg[rel]].append(unique_id)

ibd_list = {"2nd":[],"3rd":[],"4th":[]}
for pairs in open("new" + option.king).readlines()[1:]:
	pairs = pairs.split()
	unique_id,ibd1 = pairs[0],float(pairs[1])
	for degrees in ["2nd","3rd","4th"]:
		if unique_id in pair_list[degrees]:
			ibd_list[degrees].append(ibd1)

keep_list = {}
for degrees in ["2nd","3rd","4th"]:
	std = st.stdev(ibd_list[degrees])
	avg = st.mean(ibd_list[degrees])
	lower,upper = avg - 2.32635 * std, avg + 2.32635 * std
	keep = []
	for pairs in open("new" + option.king).readlines()[1:]:
		pairs = pairs.split()
		unique_id,ibd1,ibd2 = pairs[0],float(pairs[1]),float(pairs[2])
		if unique_id in pair_list[degrees] and lower <= ibd1 <= upper:
			keep.append([unique_id,ibd1,ibd2])
		#elif unique_id in pair_list[degrees] and not lower <= ibd1 <= upper:
			#print(unique_id,ibd1,lower,upper)
	keep_list[degrees] = keep

out = open("new" + option.pedigree,"w")
out.write("ID ActualDeg PONDdeg\n")
for degrees1 in ["2nd","3rd","4th"]:
	for iid1 in keep_list[degrees1]:
		test_val = [iid1[1],iid1[2]]
		unique_id = iid1[0]
		train_val,train_lab = [],[]
		for degrees2 in ["2nd","3rd","4th"]:
			for iid2 in keep_list[degrees2]:
				if iid2 != iid1:
					train_val.append([iid2[1],iid2[2]])
					train_lab.append(degrees2)
		lda_clf = LinearDiscriminantAnalysis().fit(train_val,train_lab)
		pred_deg = lda_clf.predict([test_val])[0]
		out.write(" ".join([unique_id,degrees1,pred_deg,"\n"]))
out.close()

second_rels = {"PHS":[],"MHS":[],"AV":[],"GP":[]}
for pairs in open(option.training).readlines()[1:]:
	pairs = pairs.split()
	iid1,iid2,ratio,num,rel = pairs[0],pairs[2],float(pairs[4]),int(pairs[5]),pairs[6]
	unique_id = pair_id(iid1,iid2)
	second_rels[rel].append([unique_id,ratio,num])

out = open("new" + option.training,"w")
out.write("ID ActualRel SimpleRel PONDrel\n")
for rels1 in ["PHS","MHS","GP","AV"]:
	for iid1 in second_rels[rels1]:
		test_val = [iid1[1],iid1[2]]
		unique_id = iid1[0]
		train_val,train_lab = [],[]
		for rels2 in ["PHS","MHS","GP","AV"]:
			for iid2 in second_rels[rels2]:
				if iid2 != iid1:
					train_val.append([iid2[1],iid2[2]])
					train_lab.append(rels2)
		lda_clf = LinearDiscriminantAnalysis().fit(train_val,train_lab)
		pred_rel = lda_clf.predict([test_val])[0]
		if "HS" in rels1:
			ersa_rel = "HS"
		if "HS" in pred_rel:
			ersa_pred_rel = "HS"
		out.write(" ".join([unique_id,rels1,ersa_rel,pred_rel,ersa_pred_rel,"\n"]))
out.close()

