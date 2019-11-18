from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from optparse import OptionParser
import sys
import time
import statistics as stat
import os

parser = OptionParser()
#Required args
parser.add_option("--map",dest="map_file")
parser.add_option("--king",dest="king_file")
parser.add_option("--match",dest="match_file")
parser.add_option("--fam",dest="fam_file")
#Optional args
parser.add_option("--ped",dest="ped_file",default="")
parser.add_option("--ages",dest="age_file",default="")
parser.add_option("--out",dest="out",default="PONDEROSA")
parser.add_option("--mhs_gap",dest="mhs_gap",default=30)
parser.add_option("--gp_gap",dest="gp_gap",default=30)
parser.add_option("--po_gap",dest="po_gap",default=15)
parser.add_option("--haps",dest="haps",default="")
parser.add_option("--ilash",action="store_true",dest="ilash")

(option,args) = parser.parse_args()


degrees = ["PO","FS","2nd","3rd","4th"]
ped_rels = {"PO":["PO"],
            "FS":["FS"],
            "2nd":["PHS","MHS","GP","AV"],
            "3rd":["GGP","HV","CO"],
            "4th":["GGGP","HC"]}

def ilash2germline(match_file):
    for i in range(1,23):
        ilash = open(match_file  % str(i)).readlines()
        germline = open("GERMLINE" + match_file % str(i),"w")
        for lines in ilash:
            lines = lines.split()
            lines[1] = lines[1].split("_")[0] + "." + lines[1].split("_")[1]
            lines[3] = lines[3].split("_")[0] + "." + lines[3].split("_")[1]
            lines[7] = "\n"
            cols = [0,1,2,3,4,5,6,7]
            germline.write(" ".join([lines[i] for i in cols]))
    return "GERMLINE" + match_file

def find_hap_score1(rel_list,match_file,map_file,ped_file,out,ilash):
    class GenotypeData:
        def __init__(self,map_file,ped_file,chrm):
            self.gts = dict()
            if ped_file != "":
                for ind_gts in open(ped_file % str(chrm)).readlines():
                    ind_gts = ind_gts.split()
                    iid,gts = ind_gts[0],ind_gts[6:]
                    gt_out = list()
                    for alleles in range(0,len(gts),2):
                        pair=(gts[alleles],gts[alleles+1])
                        gt_out.append(pair)
                    self.gts[iid] = gt_out

            self.snp_pos = list()
            self.mb_cm = dict()
            for snps in open(map_file % str(chrm)).readlines():
                snps = snps.split()
                cm, mb = float(snps[2]),int(snps[3])
                self.snp_pos.append(cm)
                self.mb_cm[mb] = cm

        def mb_to_cm(self,mb):
            return self.mb_cm[mb]

        def gap_discordance(self,iid1,iid2,gap_start,gap_end,threshold=1):
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
            return not status


    class PairData:
        def __init__(self,rel_list,out):
            self.pair_data = dict()
            self.keep_rels = list()
            for pairs in rel_list:
                iid1,iid2 = pairs[0],pairs[1]
                temp = [iid1,iid2]
                temp.sort()
                self.keep_rels.append((temp[0],temp[1]))
                #pair maps to {hap tots}, total ibd, num segs
                self.pair_data[(temp[0],temp[1])] = [{iid1:0,iid2:0},0,0]
            self.out = open("%s_haplotype_scores.txt" % out,"w")

        def order_pair(self,iid1,iid2):
            temp = [iid1,iid2]
            temp.sort()
            return (temp[0],temp[1])

        def move_on(self,iid1,iid2):
            pair = self.order_pair(iid1,iid2)
            return pair not in self.keep_rels

        def finish_chrm(self,iid1,iid2,seg_list,num_segs):
            hap1_0,hap1_1,hap2_0,hap2_1,tot = 0,0,0,0,0
            for segs in seg_list:
                l = segs[3] - segs[2]
                tot += l
                if segs[4] == 0:
                    hap1_0 += l
                elif segs[4] == 1:
                    hap1_1 += l
                if segs[5] == 0:
                    hap2_0 += l
                elif segs[5] == 1:
                    hap2_1 += l
            pair = self.order_pair(iid1,iid2)
            self.pair_data[pair][1] += tot
            self.pair_data[pair][2] += num_segs
            hap1 = max(hap1_0,hap1_1)
            hap2 = max(hap2_0,hap2_1)
            self.pair_data[pair][0][iid1] += hap1
            self.pair_data[pair][0][iid2] += hap2

        def get_scores(self,iid1,iid2):
            pair = self.order_pair(iid1,iid2)
            info = self.pair_data[pair]
            h1,h2,total,num = info[0][iid1],info[0][iid2],info[1],info[2]
            h1,h2 = h1/total,h2/total
            ratio = min(h1/h2,h2/h1)
            return [iid1,h1,iid2,h2,ratio,total,num]

        def write_out(self,rel_list):
            self.out.write(" ".join(["IID1","Hap1","IID2","Hap2","Ratio","NumSegs","\n"]))
            out_list = list()
            for pairs in rel_list:
                iid1,iid2 = pairs[0],pairs[1]
                scores = self.get_scores(iid1,iid2)
                self.out.write(" ".join([scores[0],str(round(scores[1],3)),scores[2],str(round(scores[3],3)),str(round(scores[4],3)),str(scores[6]),"\n"]))
                out_list.append(scores)
            return out_list

    pairs = PairData(rel_list,out)

    if ilash:
        match_file = ilash2germline(match_file)

    for chrm in range(1,23):
        ticker_char = {0: "|", 1: "/", 2: "-", 3: "\\"}
        ticker = 1
        sys.stdout.write("\rCalculating hap score...chr %s %s" % (chrm,ticker_char[ticker % 4]))
        sys.stdout.flush()
        genotype_data = GenotypeData(map_file,ped_file,chrm)
        match = list()
        for segs in open(match_file % str(chrm)).readlines():
            iid1,iid2 = segs.split()[1].split(".")[0],segs.split()[3].split(".")[0]
            if pairs.move_on(iid1,iid2):
                continue
            hap1,hap2 = int(segs.split()[1].split(".")[1]),int(segs.split()[3].split(".")[1])
            mb_start,mb_end = int(segs.split()[5]),int(segs.split()[6])
            cm_start,cm_end = genotype_data.mb_to_cm(mb_start),genotype_data.mb_to_cm(mb_end)
            match.append([iid1,iid2,cm_start,cm_end,hap1,hap2])
            ticker += 1
            sys.stdout.write("\rCalculating hap score...chr %s %s" % (chrm,ticker_char[ticker % 4]))
            sys.stdout.flush()
        if ilash:
            os.remove(match_file % str(chrm))
        match.sort()

        while match != []:
            pair_segs = list()
            line = match[0]
            iid1,iid2 = line[0],line[1]
            while iid1 == match[0][0] and iid2 == match[0][1]:
                pair_segs.append(match[0])
                match = match[1:]
                if match == []:
                    break
            ticker += 1
            sys.stdout.write("\rCalculating hap score...chr %s %s" % (chrm,ticker_char[ticker % 4]))
            start1,end1,index_list,num_segs = pair_segs[0][2],pair_segs[0][3],[0],0
            for i in range(1,len(pair_segs)):
                start2,end2 = pair_segs[i][2],pair_segs[i][3]
                if start1 == start2:
                    del index_list[-1]
                elif end2 <= end1:
                    continue
                elif (start2 > end1 + 1) or (0 < start2 - end1 <= 1 and genotype_data.gap_discordance(iid1,iid2,end1,start2)):
                    num_segs += 1
                end1 = end2
                index_list.append(i)
            num_segs += 1
            pairs.finish_chrm(iid1,iid2,[pair_segs[i] for i in index_list],num_segs)

    sys.stdout.write("\rCalculating hap score...done   ")
    return pairs.write_out(rel_list)

def find_hap_score2(haps):
    sys.stdout.write("\rHap score file provided...")
    sys.stdout.flush()
    haps = open(haps).readlines()[1:]
    out_file = list()
    for pairs in haps:
        pairs = pairs.split()
        scores = [pairs[0],float(pairs[1]),pairs[2],float(pairs[3]),float(pairs[4]),0,int(pairs[5])]
        out_file.append(scores)
    sys.stdout.write("\rHap score file provided...done\n")
    return out_file

class Data:
    def __init__(self,out):
        #Output files
        self.pedigree = open("%s.pedigree" % out,"w")
        self.score = open("%s.score" % out,"w")
        self.training = open("%s.training" % out,"w")
        self.log = open("%s.log" % out,"w")
        self.fam = open("%s_updated.fam" % out,"w")

        #Init output files
        self.pedigree.write("\t".join(["IID1","IID2","Rel","Confidence"]) + "\n")
        self.score.write("\t".join(["IID1","HapScore","IID2","HapScore","Ratio","NumSegs","P(PHS)","P(MHS)","P(GP)","P(AV)"]) + "\n")
        self.training.write("\t".join(["IID1","HapScore","IID2","HapScore","Ratio","NumSegs","Rel"]) + "\n")

        #Lists of various individuals/pairs of ind
        #List of all individuals (either gt'd or not)
        self.ind_list = list()
        #List of rels to run through hap score finder
        self.run_haps = list()
        #List of sibs who may be FS or HS
        self.unresolved_sibs = list()
        #Need to make sure that sibs w/o BOTH parents are counted
        self.king_sibs = list()
        #Dict of phenotypic data of individuals
        #iid maps to [gt status,age,sex,father,mother,[children]]
        self.individuals = dict()

        #Dict of pairs and their mutual info
        #pair maps to [ibd1,ibd2,{iid1:hap1,iid2:hap2},ratio,num,total]
        self.pairs = dict()

        #Dict keeps track of all inferred relatives
        #Degree --> ped rel --> list of pairs
        self.relationships = dict()
        for degs in degrees:
            for rels in ped_rels[degs]:
                self.relationships[rels] = []
        self.dummy_id = 0

        #Gap thresholds
        self.po_gap = option.po_gap
        self.gp_gap = option.gp_gap
        self.mhs_gap = option.mhs_gap

    #To make sure that the pair is in the right order for dict lookup
    def order_pair(self,iid1,iid2):
        temp = [iid1,iid2]
        temp.sort()
        return (temp[0],temp[1])

    def add_king_line(self,line):
        line = line.split()
        fid1,iid1,fid2,iid2,ibd1,ibd2 = line[0],line[1],line[2],line[3],float(line[6]),float(line[7])
        if ibd2 > 0.95:
            ibd1,ibd2 = 0.50,0.25
        if "FS" in line or "Dup/MZ" in line:
            self.king_sibs.append((iid1,iid2))
        self.pairs[self.order_pair(iid1,iid2)] = [ibd1,ibd2]
        if iid1 not in self.individuals:
            self.ind_list.append(iid1)
            self.individuals[iid1] = [1,"NA",0,"0","0",[],fid1]
        if iid2 not in self.individuals:
            self.ind_list.append(iid2)
            self.individuals[iid2] = [1,"NA",0,"0","0",[],fid2]
        if ibd1 > 0.1875:
            self.run_haps.append((iid1,iid2))

    def compute_hap_scores(self,haps):
        if haps == "":
            hap_scores = find_hap_score1(self.run_haps,option.match_file,option.map_file,option.ped_file,option.out,option.ilash)
        else:
            hap_scores = find_hap_score2(haps)
        for pairs in hap_scores:
            iid1,hap1,iid2,hap2,ratio,total,num = pairs[0],pairs[1],pairs[2],pairs[3],pairs[4],pairs[5],pairs[6]
            pair = self.order_pair(iid1,iid2)
            self.pairs[pair].append({iid1:hap1,iid2:hap2})
            self.pairs[pair].append([ratio,num])
            self.pairs[pair].append(total)

    #Add fam line (add IID, mom, dad to the data structure)
    def add_fam_line(self,line):
        line = line.split()
        fid,iid,dad,mom,sex = line[0],line[1],line[2],line[3],int(line[4])
        if iid in self.individuals:
            self.individuals[iid][2],self.individuals[iid][3],self.individuals[iid][4],self.individuals[iid][6] = sex,dad,mom,fid
        else:
            self.ind_list.append(iid)
            self.individuals[iid] = [0,"NA",sex,dad,mom,[],fid]

        def add_parent(parent_sex):
            parent = line[parent_sex + 1]
            if parent != "0":
                if parent in self.individuals:
                    self.individuals[parent][5].append(iid)
                else:
                    self.ind_list.append(parent)
                    self.individuals[parent] = [0,"NA",parent_sex,"0","0",[iid],fid]
        add_parent(1)
        add_parent(2)

    def add_age(self,line):
        line = line.split()
        iid,age = line[0],int(line[1])
        if iid in self.individuals:
            self.individuals[iid][1] = age

    #If true, individual is genotyped
    def gt_status(self,iid):
        return self.individuals[iid][0] == 1

    #Can return IBD1 or IBD2, depending on what is specified 
    def get_ibd(self,iid1,iid2,type):
        if self.gt_status(iid1) and self.gt_status(iid2):
            pair = self.order_pair(iid1,iid2)
            if pair in self.pairs:
                ibd = self.pairs[pair][type-1]
            else:
                ibd = 0.0
        else:
            ibd = "NA"
        return ibd

    #Get tot num of IBD segs and hap ratio for a pair
    def get_hap_data(self,iid1,iid2):
        pair = self.order_pair(iid1,iid2)
        try:
            data = self.pairs[pair][3]
        except:
            data = ["NA","NA"]
        return data

    #Get hap score for IID1 (in relation to IID2)
    def get_hap_score(self,iid1,iid2):
        pair = self.order_pair(iid1,iid2)
        try:
            score = self.pairs[pair][2][iid1]
        except:
            score = "NA"
        return score

    def get_dad(self,iid):
        return self.individuals[iid][3]

    def get_mom(self,iid):
        return self.individuals[iid][4]

    def get_children(self,iid):
        return self.individuals[iid][5]

    def get_age(self,iid):
        return self.individuals[iid][1]

    def get_sex(self,iid):
        return self.individuals[iid][2]

    def get_fid(self,iid):
        return self.individuals[iid][6]

    def get_individuals(self):
        return self.ind_list

    def add_relationship(self,iid1,iid2,rel,prob):
        prob = round(prob,3)
        if not ((iid1,iid2) in self.relationships[rel] or (iid2,iid1) in self.relationships[rel]):
            self.relationships[rel].append((iid1,iid2))
            self.pedigree.write("\t".join([iid1,iid2,rel,str(prob)]) + "\n")

    #Recursive function that finds all lineal relatives
    def lineal_relationship(self,grandparent,iid,gen,prob):
        clist = self.get_children(iid)
        for child in clist:
            if gen == 0:
                self.add_relationship(child,grandparent,"PO",prob)
            elif gen > 0:
                self.add_relationship(child,grandparent,gen*"G"+"P",prob)
            self.lineal_relationship(grandparent,child,gen+1,prob)

    #Det whether IID1 and IID2 have same mother (0 = ambiguous, 1 = same mom, 2 = diff mom)
    def mom_status(self,iid1,iid2):
        if self.get_mom(iid1) == self.get_mom(iid2):
            if self.get_mom(iid1) == "0":
                status = 0
            else:
                status = 1
        else:
            status = 2
        return status

    def dad_status(self,iid1,iid2):
        if self.get_dad(iid1) == self.get_dad(iid2):
            if self.get_dad(iid1) == "0":
                status = 0
            else:
                status = 1
        else:
            status = 2
        return status

    def add_siblings(self,iid1,iid2,prob,rel="NA"):
        if rel == "NA":
            dad,mom = self.dad_status(iid1,iid2),self.mom_status(iid1,iid2)
            if dad == 1 and mom == 1:
                rel = "FS"
            elif dad == 1 and mom == 2:
                rel = "PHS"
            elif dad == 2 and mom == 1:
                rel = "MHS"
            elif dad == 0:
                self.unresolved_sibs.append((iid1,iid2,"MHS"))
            elif mom == 0:
                self.unresolved_sibs.append((iid1,iid2,"PHS"))
        if rel != "NA":
            self.add_relationship(iid1,iid2,rel,prob)
            if rel == "FS":
                co,av = "CO","AV"
            else:
                co,av = "HC","HV"
            clist1,clist2 = self.get_children(iid1),self.get_children(iid2)
            for child1 in clist1:
                self.add_relationship(child1,iid2,av,prob)
                for child2 in clist2:
                    self.add_relationship(child1,child2,co,prob)
            for child2 in clist2:
                self.add_relationship(child2,iid1,av,prob)

    def get_relatives(self,rel):
        return self.relationships[rel]

    def get_unresolved_sibs(self):
        return self.unresolved_sibs

    def get_king_unresolved_sibs(self):
        out_sibs = list()
        for sibs in self.king_sibs:
            iid1,iid2 = sibs[0],sibs[1]
            dad1,dad2,mom1,mom2 = self.get_dad(iid1),self.get_dad(iid2),self.get_mom(iid1),self.get_mom(iid2)
            if [dad1,dad2,mom1,mom2] == ["0","0","0","0"]:
                out_sibs.append((iid1,iid2))
        return out_sibs

    def get_dummyID(self):
        self.dummy_id += 1
        return "Missing" + "%03d" % self.dummy_id

    def resolve_siblings(self):
        fs_sets = list()
        full_sibs = self.get_relatives("FS")
        for sibs in full_sibs:
            sib1,sib2 = sibs[0],sibs[1]
            new_set = True
            for sets in fs_sets:
                if sib1 in sets and sib2 not in sets:
                    new_set =  False
                    sets.append(sib2)
                    break
                elif sib2 in sets and sib1 not in sets:
                    new_set = False
                    sets.append(sib1)
                    break
            if new_set:
                fs_sets.append([sib1,sib2])

        def resolve(parent_list):
            dummy = True
            for parents in parent_list:
                if parents != "0":
                    parent = parents
                    dummy = False
                    break
            if dummy:
                parent = self.get_dummyID()
            return parent
        for sets in fs_sets:
            dad = resolve([self.get_dad(iid) for iid in sets])
            mom = resolve([self.get_mom(iid) for iid in sets])
            fid = self.get_fid(sets[0])
            for iid in sets:
                self.add_fam_line(" ".join([fid,iid,dad,mom,str(self.get_sex(iid)),"-9"]))
                if "Missing" in dad:
                    self.add_fam_line(" ".join([fid,dad,"0","0","1","-9"]))
                if "Missing" in mom:
                    self.add_fam_line(" ".join([fid,mom,"0","0","2","-9"]))
        for iid in self.ind_list:
            fid,dad,mom,sex = self.get_fid(iid),self.get_dad(iid),self.get_mom(iid),self.get_sex(iid)
            self.fam.write(" ".join([fid,iid,dad,mom,str(sex),"-9","\n"]))


    def generate_training_data(self,type):
        lab_list,val_list = [],[]
        if type == "2nd":
            for labs,rels in enumerate(ped_rels["2nd"]):
                for pairs in self.get_relatives(rels):
                    iid1,iid2 = pairs[0],pairs[1]
                    hap_data = self.get_hap_data(iid1,iid2)
                    h1,h2 = self.get_hap_score(iid1,iid2),self.get_hap_score(iid2,iid1)
                    lab_list.append(labs)
                    val_list.append(hap_data)
                    self.training.write("\t".join([iid1,str(h1),iid2,str(h2),str(hap_data[0]),str(hap_data[1]),rels,"\n"]))
            lower_list,upper_list = [],[]

        elif type == "IBD":
            ibd1_vals = {0:[],1:[],2:[],3:[],4:[]}
            for labs,degs in enumerate(degrees):
                for rels in ped_rels[degs]:
                    for pairs in self.get_relatives(rels):
                        iid1,iid2 = pairs[0],pairs[1]
                        ibd1,ibd2 = self.get_ibd(iid1,iid2,1),self.get_ibd(iid1,iid2,2)
                        ibd1_vals[labs].append(ibd1)
                        lab_list.append(labs)
                        val_list.append([ibd1,ibd2])
            #print(ibd1_vals)
            new_labs,new_vals = [],[]
            lower_list,upper_list = [],[]
            for i in range(5):
                avg = stat.mean(ibd1_vals[i])
                sd = stat.stdev(ibd1_vals[i])
                upper,lower = avg+(2.32635*sd),avg-(2.32635*sd)
                lower_list.append(lower)
                upper_list.append(upper)
                for index,labs in enumerate(lab_list):
                    if labs == 0 and i == 0:
                        new_labs.append(i)
                        new_vals.append(val_list[index])
                        continue
                    if labs == i and lower < val_list[index][0] < upper:
                        new_labs.append(i)
                        new_vals.append(val_list[index])
            val_list,lab_list = new_vals,new_labs
        return val_list,lab_list,lower_list,upper_list

    def get_putative_second(self):
        out_list = list()
        for pairs in self.run_haps:
            iid1,iid2 = pairs[0],pairs[1]
            inferred = False
            for rels in ["PO","FS","PHS","MHS","GP","AV"]:
                if ((iid1,iid2) in self.relationships[rels]
                    or (iid2,iid1) in self.relationships[rels]):
                    inferred = True
                    break
            if not inferred:
                out_list.append((iid1,iid2))
        return out_list

    def infer_relationship(self,iid1,iid2,ped_prob,prob_2nd):
        rels = {0:"PHS",1:"MHS",2:"GP",3:"AV"}
        h1,h2 = self.get_hap_score(iid1,iid2),self.get_hap_score(iid2,iid1)
        def have_parent(iid,sex):
            if sex == 1:
                parent = self.get_dad(iid)
            elif sex == 2:
                parent = self.get_mom(iid)
            return parent != "0" and "Missing" not in parent
        if have_parent(iid1,1) or have_parent(iid2,1):
            ped_prob[0] = 0
        if have_parent(iid1,2) or have_parent(iid2,2):
            ped_prob[1] = 0
        age1,age2 = self.get_age(iid1),self.get_age(iid2)
        if age1 != "NA" and age2 != "NA":
            if abs(age1-age2) < self.gp_gap:
                ped_prob[2] = 0
            if abs(age1-age2) < self.mhs_gap:
                ped_prob[1] = 0
        tot = sum(ped_prob)
        ped_prob = [ped_prob[i]/tot for i in range(4)]
        ped_rel = ped_prob.index(max(ped_prob))
        prob = prob_2nd * ped_prob[ped_rel]
        ped_rel = rels[ped_rel]
        age1,age2 = self.get_age(iid1),self.get_age(iid2)
        if age1 != "NA" and age2 != "NA":
                if age1 > age2:
                    iid1,iid2 = iid2,iid1
        else:
            if h1 < h2:
                    iid1,iid2 = iid2,iid1
        if ped_rel == "GP" or ped_rel == "AV":
            self.add_relationship(iid1,iid2,ped_rel,prob)
        else:
            self.add_siblings(iid1,iid2,prob,ped_rel)

    def write_score(self,iid1,iid2,h1,h2,hap_data,probs):
        out = [iid1,str(round(h1,3)),iid2,str(round(h2,3)),str(round(hap_data[0],3)),str(hap_data[1])]
        for i in probs:
            out.append(str(round(i,3)))
        out.append("\n")
        self.score.write("\t".join(out))

    def write_log(self,line):
        self.log.write(line)

    def close_files(self):
        self.pedigree.close()
        self.score.close()
        self.training.close()
        self.log.close()
        self.fam.close()

def main(fam_file,match_file,king_file,map_file,
        ped_file,age_file,out,mhs_gap,gp_gap,po_gap,haps,ilash):
    start_time = time.time()
    data = Data(out)

    for lines in open(king_file).readlines()[1:]:
        data.add_king_line(lines)

    data.compute_hap_scores(haps)
    sys.stdout.write("\rFinding high-confidence relationships...")
    sys.stdout.flush()
    for lines in open(fam_file).readlines():
        data.add_fam_line(lines)

    if age_file != "":
        for lines in open(age_file).readlines():
            data.add_age(lines)
    for iid in data.get_individuals():
        data.lineal_relationship(iid,iid,0,1.0)
        clist = data.get_children(iid)
        for child1 in clist:
            clist = clist[1:]
            for child2 in clist:
                data.add_siblings(child1,child2,1.0)
    sys.stdout.write("\rFinding high-confidence relationships...done\n")

    #Create IBD1 vs IBD2 classifier for degrees of rel
    sys.stdout.write("\rResolving ambiguous sibships...")
    sys.stdout.flush()
    training_data = data.generate_training_data("IBD")
    ibd_lda = LinearDiscriminantAnalysis().fit(training_data[0],training_data[1])
    for pairs in data.get_unresolved_sibs():
        iid1,iid2,rel = pairs[0],pairs[1],pairs[2]
        ibd1,ibd2 = data.get_ibd(iid1,iid2,1),data.get_ibd(iid1,iid2,2)
        prob = ibd_lda.predict_proba([[ibd1,ibd2]])[0]
        if prob[1] > prob[2]:
            data.add_siblings(iid1,iid2,prob[1]/(prob[1]+prob[2]),"FS")
        else:
            data.add_siblings(iid1,iid2,prob[2]/(prob[1]+prob[2]),rel)
    for pairs in data.get_king_unresolved_sibs():
        iid1,iid2 = pairs[0],pairs[1]
        ibd1,ibd2 = data.get_ibd(iid1,iid2,1),data.get_ibd(iid1,iid2,2)
        prob = ibd_lda.predict_proba([[ibd1,ibd2]])[0]
        if prob[1] > prob[2]:
            data.add_siblings(iid1,iid2,prob[1]/(prob[1]+prob[2]),"FS")
    data.resolve_siblings()
    sys.stdout.write("\rResolving ambiguous sibships...done\n")

    #Recreate the classifier for the added pairs
    sys.stdout.write("\rClassifying putative second degree relatives...")
    sys.stdout.flush()
    training_data = data.generate_training_data("IBD")
    ibd_lda = LinearDiscriminantAnalysis().fit(training_data[0],training_data[1])
    data.write_log("Number of training pairs for 1st LDA:\n")
    for degs in degrees:
        num = 0 
        for rels in ped_rels[degs]:
            num += len(data.get_relatives(rels))
        data.write_log("%s degree: %s\n" % (degs,num))
    data.write_log("\nIBD1 98% confidence intervals:\n")
    lower_list,upper_list = training_data[2],training_data[3]
    for index,degs in enumerate(degrees):
        lower,upper = lower_list[index],upper_list[index]
        data.write_log("%s degree: (%s, %s)\n" % (degs,lower,upper))


    #Look for IBD discrepancies
    problem_pairs = list()
    for lab,deg in enumerate(degrees[:3]):
        for rels in ped_rels[deg]:
            for pairs in data.get_relatives(rels):
                iid1,iid2 = pairs[0],pairs[1]
                ibd1,ibd2 = data.get_ibd(iid1,iid2,1),data.get_ibd(iid1,iid2,2)
                if (lab == 0 and ibd1 < 0.8) or (lab == 1 and ibd1 < 0.3) or (lab == 2 and ibd1 < 0.3):
                        problem_pairs.append([iid1,iid2,rels,str(ibd1),str(ibd2),"\n"])
    if problem_pairs != []:
        data.write_log("\nThe following relationship pairs were inferred by the given .fam file. ")
        data.write_log("Their IBD values are not consistent with the given pedigree structure.\n")
        data.write_log("Please note that (1) these are 2nd/FS pairs with IBD1 < 0.30 or PO pairs with IBD1 < 0.80 ")
        data.write_log("and (2) these relationships are inconsistent because of incorrect PO assignment(s) in the .fam file.\n")
        data.write_log("Please double-check and rerun if necessary.\n")
        for pairs in problem_pairs:
            data.write_log(" ".join(pairs))

    #Look for age
    problem_pairs = list()
    for rels in ["PO","GP"]:
        for pairs in data.get_relatives(rels):
            younger,older = pairs[0],pairs[1]
            young_age,old_age = data.get_age(younger),data.get_age(older)
            if young_age != "NA" and old_age != "NA": 
                if young_age > old_age:
                    problem_pairs.append([rels,younger,str(young_age),older,str(old_age),"\n"])
                if rels == "PO" and abs(young_age-old_age) < po_gap:
                    problem_pairs.append([rels,younger,str(young_age),older,str(old_age),"\n"])
                if rels == "GP" and abs(young_age-old_age) < gp_gap:
                    problem_pairs.append([rels,younger,str(young_age),older,str(old_age),"\n"])
            else:
                h1,h2 = data.get_hap_score(younger,older),data.get_hap_score(older,younger)
                if h1 != "NA" and h2 != "NA" and h2 > h1 and h2 > 0.80:
                    problem_pairs.append([rels,younger,str(h1),older,str(h2),"\n"])
    if problem_pairs != []:
        data.write_log("The following pairs are PO or GP pairs in which the reported age is not consistent with the relationship.\n")
        data.write_log("Either the child/granchild is older than the parent/grandparent or the age gap is too small.\n")
        data.write_log("If age data is unavailable, hap scores are used. If a parent/grandparent has a hap score greater than the child/grandchild ")
        data.write_log("and the hap score is >0.80, then a warning will be printed.")

    #Create the classifier for 2nd degree relatives
    training_data = data.generate_training_data("2nd")
    second_lda = LinearDiscriminantAnalysis().fit(training_data[0],training_data[1])
    data.write_log("\nNumber of training  pairs for 2nd LDA:\n")
    for rels in ped_rels["2nd"]:
        num = len(data.get_relatives(rels))
        data.write_log("%s: %s\n" % (rels,num))

    #Determine the correct relationship for each putative second degree relative
    for pairs in data.get_putative_second():
        iid1,iid2 = pairs[0],pairs[1]
        ibd1,ibd2 = data.get_ibd(iid1,iid2,1),data.get_ibd(iid1,iid2,2)
        prob_2nd = ibd_lda.predict_proba([[ibd1,ibd2]])[0][2]
        if prob_2nd > 0.80:
            hap_data = data.get_hap_data(iid1,iid2)
            h1,h2 = data.get_hap_score(iid1,iid2),data.get_hap_score(iid2,iid1)
            ped_prob = list(second_lda.predict_proba([hap_data])[0])
            data.write_score(iid1,iid2,h1,h2,hap_data,ped_prob)
            data.infer_relationship(iid1,iid2,ped_prob,prob_2nd)
    sys.stdout.write("\rClassifying putative second degree relatives...done\n")

    data.write_log("\nRuntime: %s seconds\n" % round(time.time() - start_time,0))

main(option.fam_file,option.match_file,option.king_file,option.map_file,
        option.ped_file,option.age_file,option.out,
        option.mhs_gap,option.gp_gap,option.po_gap,option.haps,option.ilash)
