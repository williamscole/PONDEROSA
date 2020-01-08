from optparse import OptionParser
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import random

parser = OptionParser()

parser.add_option("--pedigree",dest="pedigree")
parser.add_option("--king",dest="king")
parser.add_option("--ersa",dest="ersa")

(option,args) = parser.parse_args()

class Relatives:
	def __init__(self):
		self.pairs = list()
		self.pair_data = dict()
		self.rel_to_deg = {"PO":"PO","FS":"FS",
							"PHS":"2","MHS":"2","GP":"2","AV":"2",
							"HV":"3","CO":"3","GGP":"3",
							"HC":"4","GGGP":"4"}

	def order_pair(self,iid1,iid2):
		temp = [iid1,iid2]
		temp.sort()
		return temp[0],temp[1]

	def add_pedigree(self,line):
		line = lines.split()
		iid1,iid2,rel,confidence = line[0],line[1],line[2],line[3]
		degree = self.rel_to_deg[rel]
		if confidence == "1.0":
			pair = self.order_pair(iid1,iid2)
			if pair in self.pairs:
				degree2 = self.pair_data[pair][1]
				if degree < degree2:
					self.pair_data[pair][1] = degree
					self.pair_data[pair][0] = rel
			else:
				self.pairs.append(pair)
				self.pair_data[pair] = [rel,degree,"KINGdeg","IBD1","IBD2","ERSAdeg","ERSArel"]

	def add_king(self,line):
		line = line.split()
		iid1,iid2,ibd1,ibd2,inf_deg = line[1],line[3],float(line[6]),float(line[7]),line[9]
		pair = self.order_pair(iid1,iid2)
		if pair in self.pairs:
			if inf_deg in ["2nd","3rd","4th"]:
				inf_deg = inf_deg[0]
			self.pair_data[pair][2] = inf_deg
			self.pair_data[pair][3] = ibd1
			self.pair_data[pair][4] = ibd2

	def get_ersa_rel(self,ancestors,degree):
		if degree not in ["1","2","3","4"]:
			degree = "UN"
			rel = "NA"
		if degree == "1":
			if ancestors == "0":
				degree = "PO"
				rel = "PO"
			if ancestors == "2":
				degree = "FS"
				rel = "FS"	
		if degree == "2":
			if ancestors == "0":
				rel = "GP"
			if ancestors == "1":
				rel = "HS"
			if ancestors == "2":
				rel = "AV"
		if degree == "3":
			if ancestors == "0":
				rel = "GGP"
			if ancestors == "1":
				rel = "HV"
			if ancestors == "2":
				rel = "CO"
		if degree =="4":
			rel = "NA"
		return degree,rel		
		

	def add_ersa(self,line):
		line = line.split()
		iid1,iid2,ancestors,degree = line[0],line[1],line[2],line[3]
		ersa_rel = self.get_ersa_rel(ancestors,degree)
		degree,rel = ersa_rel[0],ersa_rel[1]
		pair = self.order_pair(iid1,iid2)
		if pair in self.pairs:
			self.pair_data[pair][5] = degree
			self.pair_data[pair][6] = rel

	def create_classifier(self):
		labs,vals = [],[]
		for pairs in self.pairs:
			data = self.pair_data[pairs]
			deg,ibd1,ibd2 = data[1],data[3],data[4]
			if deg == "PO":
				deg = 0
			if deg == "FS":
				deg = 1
			labs.append(deg)
			vals.append([ibd1,ibd2])
		return LinearDiscriminantAnalysis().fit(vals,labs)

	def write_out(self):
		out = open("Comparisons.txt","w")
		out.write("ActualRel ActualDeg KINGDeg IBD1 IBD2 ERSADeg ERSARel PONDDeg\n")
		lda_clf = self.create_classifier()
		for pairs in self.pairs:
			data = self.pair_data[pairs]
			ponderosa = lda_clf.predict([[data[3],data[4]]])[0]
			if ponderosa == "0":
				ponderosa = "PO"
			if ponderosa == "1":
				ponderosa = "FS"
			data.append(ponderosa)
			data[3] = str(data[3])
			data[4] = str(data[4])

			data.append("\n")
			out.write(" ".join(data))

relatives = Relatives()

for lines in open(option.pedigree).readlines()[1:]:
	relatives.add_pedigree(lines)

for lines in open(option.king).readlines()[1:]:
	relatives.add_king(lines)

for lines in open(option.ersa).readlines()[3:]:
	relatives.add_ersa(lines)

relatives.write_out()

