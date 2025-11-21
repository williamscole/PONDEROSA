import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import argparse
import os
from collections import namedtuple
import itertools as it
import yaml
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from copy import deepcopy

from pedigree_tools import ProcessSegments, Pedigree, PedigreeHierarchy, introduce_phase_error

''' Ponderosa can take as input a file (e.g., containing IBD segments) that has all chromosomes
represented, in which case it does not expect chr1 to be in the file name. Otherwise, if the files
are split by chromosome, the user should supply the path and file name for chromosome 1, in which
case the chromosome number should be denoted as chr1'''
def get_file_names(file_name):

    if "chr1" in file_name:
        return [file_name.replace("chr1", f"chr{chrom}") for chrom in range(1, 23)]
    else:
        return [file_name]
    


''' SampleData is a class that holds a networkx graph where nodes are individuals
and edges describe the relationship between those individuals. It can be queried easily
using the get_edges and get_nodes functions, which allow the user to supply functions
that act as conditionals for returning node pairs and nodes, resp.'''
class SampleData:

    def __init__(self, fam_file, **kwargs):

        ### load the fam file
        fam = pd.read_csv(fam_file, delim_whitespace=True, dtype=str, header=None)

        ### the default information for adding a person
        self.default = {"sex": 0, "mother": -1, "father": -1, "children": [], "age": np.nan, "popl": "pop1"}

        ### init the graph and set defaults
        g = nx.Graph()
        g.add_nodes_from(fam.apply(lambda x: (x[1], {"sex": x[4],
                                                    "mother": x[3] if x[3] != "0" else -1,
                                                    "father": x[2] if x[2] != "0" else -1,
                                                    "children": [],
                                                    "age": np.nan,
                                                    "popl": "pop1"}), axis=1).values)

        ### get the population
        popln = kwargs.get("population", "pop1")

        # add the children
        for sex_col in [2, 3]:
            # assume that missing parent is coded as 0
            for parent, children_df in fam[fam[sex_col] != "0"].groupby(sex_col):
                # get a list of the children
                children = children_df[1].values.tolist()
                # add to attribute
                nx.set_node_attributes(g, {parent: children}, "children")

        # add the population data
        pop_file = kwargs.get("pop_file", "")
        if os.path.exists(pop_file):
            pops = pd.read_csv(pop_file, delim_whitespace=True, dtype=str, header=None)
            pop_attrs = {iid: pop for iid, pop in pops[[0,1]].values}
            nx.set_node_attributes(g, pop_attrs, "popl")
        else:
            print("No population file has been provided.")

        # add age data
        age_file = kwargs.get("age_file", "")
        if os.path.exists(age_file):
            ages = pd.read_csv(age_file, delim_whitespace=True, dtype={0: str, 1: float}, header=None)
            age_attrs = {iid: age for iid, age in ages[[0,1]].values}
            nx.set_node_attributes(g, age_attrs, "age")
        else:
            print("No age data has been provided.")

        ### get the genome len
        map_file = kwargs.get("map_file", "")
        if os.path.exists(map_file):
            map_files = get_file_names(map_file)

            map_df = pd.concat([pd.read_csv(filen, delim_whitespace=True, header=None) for filen in map_files])

            genome_len = 0
            for chrom, chrom_df in map_df.groupby([0]):
                genome_len += (chrom_df.iloc[-1][2] - chrom_df.iloc[0][2])
        # no map files provided, use default genome length
        else:
            print("No map files have been provided. Assuming genome length of 3545 cM.")
            genome_len = 3545
                # store the genome_len

        self.genome_len = genome_len

        ### load the king file
        king_file = kwargs.get("king_file", "")
        mz_twins = set()
        if os.path.exists(king_file):
            master_hier = PedigreeHierarchy("tree_codes.yaml")

            king_df = pd.read_csv(king_file, delim_whitespace=True, dtype={"ID1": str, "ID2": str})

            # get set of MZ twins; remove the 2nd twin in the pair
            mz_twins = set(king_df[king_df.InfType=="Dup/MZ"]["ID2"].values.tolist())

            # create an IBD graph and add var attrs to it
            g.add_edges_from(king_df.apply(lambda x: [x.ID1, x.ID2, 
                                                {"k_prop": x.PropIBD,
                                                "k_ibd1": x.IBD1Seg,
                                                "k_ibd2": x.IBD2Seg,
                                                "k_degree": x.InfType,
                                                "probs": deepcopy(master_hier)}],
                                                axis=1).values)

        ### subset to only samples from the population and who are not part of the twin set
        popl_samples = {nodes for nodes, attrs in g.nodes(data=True) if attrs.get("popl", "")==popln}
        self.g = g.subgraph(popl_samples - mz_twins)


        ### load ibd, process ibd segs
        ibd_file = kwargs.get("ibd_file", "")
        if os.path.exists(ibd_file):
            print("Processing IBD segments...")
            ibd_files = get_file_names(ibd_file)
            # load ibd files
            ibd_df = pd.concat([pd.read_csv(filen, delim_whitespace=True, dtype={"id1": str, "id2": str})
                                for filen in ibd_files])

            ibd_df["l"] = ibd_df["end_cm"] - ibd_df["start_cm"]

            # iterate through the pairs of individuals, compute ibd stats
            for (id1, id2), pair_df in ibd_df.groupby(["id1", "id2"]):

                # only compute hsr, other stats for 3rd+ degree relatives
                k = self.g.get_edge_data(id1, id2, {}).get("k_prop", 0)

                # 2^-3.5 is the lower limit for 3rd degree relatives
                if k < 2**-3.5:
                    continue

                # compute the IBD data
                pair_ibd = ProcessSegments(pair_df)
                ibd_data = pair_ibd.ponderosa_data(genome_len, inter_phase=False)

                # add phase errors
                pair_ibd_pe = ProcessSegments(introduce_phase_error(pair_df, 50))
                ibd_data_pe = pair_ibd_pe.ponderosa_data(genome_len, inter_phase=False)

                # set up the pedigree hierarchy, which will store the probs and info on how the probs were computed
                hier = deepcopy(master_hier)
                # TODO implement this somewhere else
                # hier.update_attr_from([[node, "ordering", sorted([id1, id2])] for node in hier.init_nodes])
                # hier.update_attr_from([[node, "ordering_method", "sorted"] for node in hier.init_nodes])

                cur_edge_data = self.g.get_edge_data(id1, id2)


                cur_edge_data["ibd1"] = ibd_data.ibd1
                cur_edge_data["ibd2"] = ibd_data.ibd2
                cur_edge_data["h"] = {id1: ibd_data.h1, id2: ibd_data.h2}
                cur_edge_data["h_error"] = {id1: ibd_data_pe.h1, id2: ibd_data_pe.h2}
                cur_edge_data["n"] = ibd_data.n
                cur_edge_data["probs"] = hier


                # # add ibd1 data and initialze probs
                # self.g.add_edge(id1, id2, ibd1=ibd_data.ibd1,
                #                      ibd2=ibd_data.ibd2,
                #                      h={id1: ibd_data.h1, id2: ibd_data.h2},
                #                      h_pe={id1: ibd_data_pe.h1, id2: ibd_data_pe.h2},
                #                      n=ibd_data.n,
                #                      k=(ibd_data.ibd1/2 + ibd_data.ibd2),
                #                      probs=hier,
                #                      segments=pair_df)

        else:
            print("No IBD files have been provided.")


    # returns a list of node pair edges
    # optional func arg is a function that takes as input the data dict and returns true if wanted
    def get_edges(self, f=lambda x: True):
        return [(id1, id2) for id1, id2, data in self.g.edges(data=True) if f(data)]

    # returns a nested np array of [id1, id2, attr1, attr2, ...] for all edges that pass the function
    def get_edge_data(self, attr_list, f=lambda x: True):
        edges, out_attrs = [], []
        for id1, id2, data in self.g.edges(data=True):
            if f(data):
                edges.append([id1, id2])
                out_attrs.append([data[attr] for attr in attr_list])
        return np.array(edges, dtype=str), np.array(out_attrs)

    # updates a bunch of edges at once
    # edge_list looks like [(A, B), (B, C), (D, E)]
    # attr_list looks like [1, 2, 3]
    def update_edges(self, attr_list, attr_name):

        attr_dict = {(id1, id2): attr for (id1, id2), attr in attr_list}
        nx.set_edge_attributes(self.g, attr_dict, attr_name)

    # same as above but returns a set of nodes
    # e.g., get_nodes(lambda x: x["age"] > 90) returns all nodes that are >90 yo
    def get_nodes(self, f=lambda x: True):
        return [id1 for id1, data in self.g.nodes(data=True) if f(data)]
    
    def get_subset(self, f):
        # get the nodes based on the conditional f
        node_subset = self.get_nodes(f)
        # return the subset for which the conditional is true
        return self.g.subgraph(set(node_subset))
    
    # subsets or returns a subset of the g based on a iterable of edges; inplace=True means it will modify self.g
    def edge_subset(self, edges, inplace):
        if inplace:
            self.g = self.g.edge_subgraph(edges).copy()

        else:
            return self.g.edge_subgraph(edges).copy()

    
    def to_dataframe(self, edges, include_edges):
        if len(edges) == 0:
            return nx.to_pandas_edgelist(self.g, source="id1", target="id2")
        # make copy of the graph
        tmp = self.g.copy()
        # subgraph of the edges
        if include_edges:
            tmp = tmp.edge_subgraph(edges)
        # subgraph of the complement of the edges
        else:
            tmp.remove_edges_from(edges)
        return nx.to_pandas_edgelist(tmp, source="id1", target="id2")

    # update the hier probabilities
    def update_probs(self, probs_list, prob_labs, method, prob_objs):
        # iterate through the probabilities
        for probs, prob_obj in zip(probs_list, prob_objs):
            prob_obj.set_attrs({i: p for i,p in zip(prob_labs, probs)}, "p_con")
            prob_obj.set_attrs({i: method for i in prob_labs}, "method")

    def node_attr(self, node, attr, nan=np.nan):
        try:
            return self.g.nodes[node][attr]
        except:
            return nan


class Classifiers:
    def __init__(self, pairs):

        # store the training data for each classifier; purpose is to allow leave one out training
        self.training = {}

        ### Train the degree classifier
        degree_train = pairs.get_pair_df("relatives").dropna(subset=["k_ibd1"])

        self.training["degree"] = [np.array(degree_train[["k_ibd1", "k_ibd2"]].values.tolist()),
                                   np.array(degree_train["degree"].values.tolist())]
        
        ### Train the hap classifier
        hap_train = pairs.get_pair_df_from(["HS", "GPAV"]).dropna(subset=["h"])

        # get the phase error classifier
        X_train = hap_train.apply(lambda x: [x.h_error[x.pair[0]], x.h_error[x.pair[1]]], axis=1).values.tolist()
        y_train = ["Phase error" for _ in X_train]

        # now add the actual haplotype scores
        X_train += hap_train.apply(lambda x: [x.h[x.pair[0]], x.h[x.pair[1]]], axis=1).values.tolist()
        y_train += hap_train["requested"].values.tolist()

        self.training["hap"] = [np.array(X_train), np.array(y_train)]

        ### Train the degree classifier
        n_train = pairs.get_pair_df_from(["MGP", "MHS", "PHS", "PGP", "AV"]).dropna(subset=["n"])

        # also train on the kinship coefficient
        n_train["ibd_cov"] = n_train.apply(lambda x: x.k_ibd2 + x.k_ibd1, axis=1)

        self.training["n"] = [np.array(n_train[["ibd_cov", "n"]].values.tolist()), np.array(n_train["requested"].values.tolist())]

        self.loo = LeaveOneOut()

    # if n samples in X_train, will train the classifier n times, leaving a different sample out each time
    def leaveOneOut(self, X_train, y_train, func, labels):
            return_probs = []; return_labels = []
            lda = LinearDiscriminantAnalysis()

            for train, test in self.loo.split(X_train):

                lda.fit(X_train[train], y_train[train])

                # predict
                return_probs.append(func(lda, X_train[test])[0])
                return_labels.append(lda.classes_)

            # just predicting the label; don't need the probabilities
            if labels:
                return return_probs, iter(return_labels)

            return return_probs

    def return_classifier(self, classif):

        X_train, y_train = self.training[classif]

        lda = LinearDiscriminantAnalysis()

        lda.fit(X_train, y_train)

        return lda
    
    # returns an array of the probabilities and an iterator of the classes that each probability describes
    def predict_proba(self, classif, X=[]):

        # get the training data for the classifier
        X_train, y_train = self.training[classif]

        # no training data supplied --> perform leave one out
        if len(X)==0:
            return self.leaveOneOut(X_train, y_train, lambda lda, X: lda.predict_proba(X), labels=True)

        lda = self.return_classifier(classif)
 
        return lda.predict_proba(X), it.cycle([lda.classes_])

    # predicts just the label
    def predict(self, classif, X=[]):
        lda = LinearDiscriminantAnalysis()
        X_train, y_train = self.training[classif]

        # perform leave-one-out
        if len(X)==0:
            return self.leaveOneOut(X_train, y_train, lda, lambda lda, X: lda.predict(X), labels=False)

        lda = self.return_classifier(classif)

        return lda.predict(X)
    
    def write_pkl(self, classif, output):

        lda = self.return_classifier(classif)

        i = open(output, "wb")
        pkl.dump(lda, i)
        i.close()


class ResultsData:
    def __init__(self, samples, pairs, df=pd.DataFrame()):

        self.samples = samples
        self.pairs = pairs
        self.df = self.samples.to_dataframe([], include_edges=False) if df.shape[0]==0 else df

    # writes out a human readable output; can specificy the columns
    def write_readable(self, output, **kwargs):

        cols = kwargs.get("columns",
                          ["id1", "id2", "k_ibd1", "k_ibd2", "most_probable", "probability", "degree"])
        
        if "h1" in cols:
            self.df["h1"], self.df["h2"] = zip(*self.df.apply(lambda x: [x.h[x.id1], x.h[x.id2]], axis=1))

        if "degree" in cols:
            self.df["degree"] = self.df["probs"].apply(lambda x: x.degree_nodes[np.argmax([x.hier.nodes[node]["p"] for node in x.degree_nodes])])

        self.df[cols].to_csv(output, index=False, sep="\t")

    # creates the dataframe
    def to_dataframe(self, edges, include_edges, inplace=False):

        tmp = self.samples.to_dataframe(edges, include_edges)

        if inplace:
            self.df = tmp
        else:
            return tmp

    # takes as input a function that works on the dataframe and subsets it
    def subset_dataframe(self, func, inplace=False):
        tmp = self.df[self.df.apply(lambda x: func(x), axis=1)]

        if inplace:
            self.df = tmp
        else:
            return tmp
        
    def subset_samples(self, func):

        self.subset_dataframe(func, inplace=True)

        self.samples.edge_subset(self.df[["id1", "id2"]].apply(tuple, axis=1).values, inplace=True)

    # sets the min prob for readable format; can be rerun with new probs
    def set_min_probability(self, min_p, update_attrs=False):
        self.df["most_probable"], self.df["probability"] = zip(*self.df["probs"].apply(lambda x: x.most_probable(min_p)))

        if update_attrs:
            self.samples.update_edges(self.df[["id1","id2"]].values, self.samples["most_probable"].values, "most_probable")
            self.samples.update_edges(self.df[["id1","id2"]].values, self.samples["most_probable"].values, "probability")

    # recomputes probs across all rows of df
    def compute_probs(self):
        self.df["probs"].apply(lambda x: x.compute_probs())

    # pickles the object
    def pickle_it(self, output):
        self.df = None
        # self.samples = None
        i = open(output, "wb")
        pkl.dump(self, i)
        i.close()

    def most_likely_among(self, nodes, update_attrs=False):
        likely_among = self.df["probs"].apply(lambda x: nodes[np.argmax([x.hier.nodes[node]["p"] for node in nodes])]).values

        if update_attrs:
            self.samples.update_edges(zip(self.df[["id1","id2"]].values, likely_among), "likely_among")

        else:
            return likely_among



def PONDEROSA(samples, **kwargs):

    pedigree = Pedigree(samples=samples, pedigree_file="pedigree_codes.yaml", tree_file="tree_codes.yaml")
    pedigree.find_all_relationships()

    pairs = pedigree.hier

    training = Classifiers(pairs)

    if kwargs.get("assess", False):

        unknown_df = samples.to_dataframe(pairs.get_pairs("relatives"), include_edges=True).dropna(subset=["k_prop"])

    else:

        # get all close relatives to exclude
        found_close = list(pairs.get_nodes_from(["2nd", "PO", "FS"]))

        # make a dataframe of all pairs that are not close relatives
        unknown_df = samples.to_dataframe(found_close, include_edges=False)

    unknown_df = unknown_df.reset_index()


    unknown_df["ibd_cov"] = unknown_df.apply(lambda x: x.ibd1 + x.ibd2, axis=1)

    # predict the degree of relatedness
    # unknown_df["predicted_degree"] = training.predict("degree", unknown_df[["k_ibd1", "k_ibd2"]].values)


    # get the degree probabilities
    probs, labels = training.predict_proba("degree", unknown_df[["k_ibd1", "k_ibd2"]].values)

    # add the probabilities to the tree
    for (_, row), prob in zip(unknown_df.iterrows(), probs):
        row["probs"].add_probs(list(zip(next(labels), prob)), "ibd")

    if kwargs.get("assess", False):

        second_pairs = samples.to_dataframe(pairs.get_pairs("2nd"), include_edges=True)[["id1", "id2"]]
        second = unknown_df.merge(second_pairs, on=["id1", "id2"], how="inner")

    else:
        second = unknown_df[unknown_df["probs"].apply(lambda x: x.hier.nodes["2nd"]["p_con"] > 0.2)]


    # get the n_ibd segs classifier probabilities
    probs, labels = training.predict_proba("n", second[["ibd_cov", "n"]].values)
    # add the probabilities to the tree
    for (_, row), prob in zip(second.iterrows(), probs):
        row["probs"].add_probs(list(zip(next(labels), prob)), "nsegs")

        # for each of these nodes, take the sum of the two children
        for node in ["GP", "GPAV", "HS"]:
            # array of the children probabilities
            child_probs = [row["probs"].hier.nodes[i]["p_con"] for i in row["probs"].hier.successors(node)]
            # add the probability
            row["probs"].add_probs(node, p_con=sum(child_probs), method="nsegs")


    # get the probabilities from the hap score classifier
    probs, labels = training.predict_proba("hap", second.apply(lambda x: sorted([h for _,h in x.h.items()], reverse=True), axis=1).values.tolist())
    
    for (_, row), prob in zip(second.iterrows(), probs):
        # get the index of the Phase error class
        classes = list(next(labels))
        pe_index = classes.index("Phase error"); del classes[pe_index]
        # the chance of high Phase error is high; do not update the HS or GPAV probabilities
        if prob[pe_index] < 0.2:
            row["probs"].add_probs(list(zip(classes, np.delete(prob, pe_index))), "hap")

    samples.edge_subset(unknown_df[["id1", "id2"]].apply(lambda x: tuple(x), axis=1).values, inplace=True)

    relatives_obj = ResultsData(samples=samples, pairs=pairs, df=unknown_df)
    relatives_obj.compute_probs()
    relatives_obj.set_min_probability(kwargs.get("min_p", 0.5))
    relatives_obj.write_readable(f"{kwargs.get('output', 'output')}.txt")
    if not kwargs.get("assess", False):
        relatives_obj.subset_samples(lambda x: x.probs.hier.nodes["2nd"]["p"] > 0.5)
    relatives_obj.pickle_it(f"{kwargs.get('output', 'output')}_results.pkl")

def parse_args():
    parser = argparse.ArgumentParser()
    # Required file arguments
    parser.add_argument("--ibd", help = "IBD segment file. If multiple files for each chromosome, this is the path for chromosome 1.")
    parser.add_argument("--fam", help="PLINK-formated .fam file")
    parser.add_argument("--king", help="KING .seg file.")

    # Optional file arguments
    parser.add_argument("--ages", help="Age file. First column is the IID, second column is the age", default="")
    parser.add_argument("--map", help = "PLINK-formatted .map file.", default="")
    parser.add_argument("--populations", help="Path and file name of .txt file where col1 is the IID and col2 is their population.", default="")
    parser.add_argument("--yaml", help="YAML file containing all arguments (optional). Can be combined with CLI arguments.", default="")
    parser.add_argument("--pedigree_codes", default="pedigree_codes.yaml")
    
    # Other arguments
    parser.add_argument("--output", help = "Output prefix.", default="Ponderosa")
    parser.add_argument("--min_p", help="Minimum posterior probability to output the relationship.", default=0.50, type=float)
    parser.add_argument("--population", help="Population name to run Ponderosa on.", default="pop1")
    parser.add_argument("--assess", help="For assessing the performance of Ponderosa.", action="store_true")

    parser.add_argument("-training", help = "Path and name of the 'degree' pkl file for training.", default="")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    return args

'''
Takes as input a yaml_file and returns an object whose attributes are various Ponderosa arguments
'''
def load_yaml_args(yaml_file, args):

    i = open(yaml_file)
    yaml_dict = yaml.safe_load(i)

    # iterate through the arguments and their values (possibly default)
    for argname, arg in vars(args).items():
        # if the arg is already in the yaml_dict, it will take it, otherwise take the default
        yaml_dict[argname] = yaml_dict.get(argname, arg)

    Args = namedtuple("Args", list(yaml_dict))
    args = Args(*[j for _,j in yaml_dict.items()])
    return args


if __name__ == "__main__":

    args = parse_args()

    # yaml file provided as the argument file
    if os.path.exists(args.yaml):
        args = load_yaml_args(args.yaml, args)

    if args.debug and os.path.exists(f"{args.output}_samples.pkl"):
        print("Samples pkl file provided.")
        i = open(f"{args.output}_samples.pkl", "rb")
        samples = pkl.load(i)
        i.close()

    # Get the samples for Ponderosa input
    else:
        samples = SampleData(fam_file=args.fam,
                    king_file=args.king,
                    ibd_file=args.ibd,
                    map_file=args.map)
        if args.debug:
            i = open(f"{args.output}_samples.pkl", "wb")
            pkl.dump(samples, i)
            i.close()

    PONDEROSA(samples=samples,
              min_p=args.min_p,
              assess=args.assess,
              output=args.output)
