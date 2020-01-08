from optparse import OptionParser
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import random

parser = OptionParser()

#parser.add_option("--pedigree",dest="pedigree")
#parser.add_option("--king",dest="king")
parser.add_option("--training",dest="training")
#parser.add_option("--ersa",dest="ersa")

(option,args) = parser.parse_args()

class SecondRelatives:
	def __init__(self):
		self.second = {0:[],1:[],2:[],3:[]}
		self.second_rels = {"PHS":0,"MHS":1,"GP":2,"AV":3}
		self.min_length = 0

	def add_training(self,line):
		line = line.split()
		ratio,num,rel = float(line[4]),int(line[5]),line[6]
		self.second[self.second_rels[rel]].append([ratio,num])

	def stdize(self):
		self.min_length = min([len(self.second[i]) for i in range(4)])

	def get_indexes(self,length,num):
		train = random.sample(range(0,length),num)
		test = []
		for i in range(0,length):
			if i not in train:
				test.append(i)
		return train,test

	def num_range(self):
		start = 5 * round((0.25 * self.min_length)/5)
		end = 5 * round((0.27 * self.min_length)/5)
		return start,end

	def get_sets(self,rel,num):
		rel_list = self.second[rel]
		keep_index = random.sample(range(0,len(rel_list)),self.min_length)
		keep_list = [rel_list[i] for i in keep_index]
		indexes = self.get_indexes(self.min_length,num)
		train_index,test_index = indexes[0],indexes[1]
		return [keep_list[i] for i in train_index],[keep_list[i] for i in test_index]

second = SecondRelatives()

for line in open(option.training).readlines()[1:]:
	second.add_training(line)

second.stdize()

for pairs in range(second.num_range()[0],second.num_range()[1],5):
	data = {0:[0,0,0,0,0],1:[0,0,0,0,0],2:[0,0,0,0,0],3:[0,0,0,0,0]}
	for iterations in range(100):
		phs,mhs,gp,av = second.get_sets(0,pairs),second.get_sets(1,pairs),second.get_sets(2,pairs),second.get_sets(3,pairs)
		train_labs,train_vals = [],[]
		for labs,rels in enumerate([phs,mhs,gp,av]):
			for vals in rels[0]:
				train_labs.append(labs)
				train_vals.append(vals)
		lda_clf = LinearDiscriminantAnalysis().fit(train_vals,train_labs)
		for labs,rels in enumerate([phs,mhs,gp,av]):
			for vals in rels[1]:
				predicted_rel = lda_clf.predict([vals])[0]
				data[labs][predicted_rel] += 1
				data[labs][4] += 1
	deg_rel = {0:"PHS",1:"MHS",2:"GP",3:"AV"}
	for i in range(4):
		for k in range(4):
			print(deg_rel[i],deg_rel[k],float(data[i][k])/data[i][4])
		'''sensitivity = float(data[i][i])/data[i][4]
		specificity = float(data[i][i])/sum([data[k][i] for k in range(4)])
		print(i,sensitivity,specificity)'''

