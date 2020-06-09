import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sys

'''Errors that can be risen here:
1) If there are not enough FS/HS training pairs for the classifier
        Action: quit PONDEROSA
2) If there is a full sib set with parent issues
        Action: ignore but append to log
3) If there are KING PO not reported in .fam
        Action: ignore but append to log
4) MZ twins
'''

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

    def add_po(self,parent,child,sex):
        if self.pedigree_structure[child][sex-1] == []: #parent not reported yet
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
            self.add_error("KING has inferred the following as FS, but have low IBD2 values (< 0.15):\n",[pairs.split("_") for pairs in low_IBD2],2)
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
                    self.add_error("Sparse pedigree warning: not enough training pairs. Try rerunning with king_fs as True\n",[],1)
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

        #Takes in a list of parents (e.g. all the mothers of a set of FS) and looks for 1) multiple mothers and 2) no mother
        def resolve_parents(parent_list):
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

        error_pairs = []
        for sets in fs_sets:
            dad = resolve_parents([self.get_parent(sibs,1) for sibs in sets])
            mom = resolve_parents([self.get_parent(sibs,2) for sibs in sets])
            if dad[1] != [] or mom[1] != []:
                error_pairs.append(sets)
                continue
            dad,mom = dad[0],mom[0]
            for sibs in sets:
                self.add_po(dad,sibs,1)
                self.add_po(mom,sibs,2)
        self.add_error("The following FS sets have different parents. The problem has been ignored but please double check.\n",error_pairs,2)

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
        self.add_error("The following PO pairs have been inferred by KING but are not reported in the .fam file.\n",missing,2)

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
            outfile.write(self.relatives.to_string(index=False,na_rep="NA"))

    def run_PONDEROSA(self,king_file,fam_file,out,trust_fs):
        self.add_king(king_file)
        self.make_from_fam(fam_file)
        self.resolve_siblings(trust_fs)
        self.get_all_pairs()
        self.print_out(out)
        self.check_po()
        return self.errors


