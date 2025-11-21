import networkx as nx
import pandas as pd
import itertools as it
import numpy as np
from datetime import datetime
import time
import argparse
import logging
import yaml
from math import floor, ceil
# import phasedibd as ibd
import os
from collections import Counter
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.mixture import GaussianMixture
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
import pickle as pkl
import sys
import concurrent.futures
import subprocess
from matplotlib import colors
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib as mpl
import heapq


def split_regions(region_dict, new_region):
    # returns the overlap of 2 regions (<= 0 if no overlap)
    def overlap(region1, region2):
        start1, end1 = region1
        start2, end2 = region2
        return min(end1,end2) - max(start1,start2)
    
    # out region will be returned; is a dict of regions mapping to members of region    
    out_region = dict()
    # overlapped keeps track of all the regions that overlap with the new region
    overlapped = {tuple(new_region[:2]):[new_region[2]]}
    # iterate through the existing regions
    for region in sorted(region_dict):
        # if overlap
        if overlap(region, new_region[:2]) > 0:
            # the regions completely overlap, just add the member and return region dict
            if tuple(region) == tuple(new_region[:2]):
                region_dict[region] += [new_region[2]]
                return region_dict
            # bc the region overlaps, add it to overlapped
            overlapped[region] = region_dict[region]
        # no overlap, but add the region to the out_region dict
        else:
            out_region[region] = region_dict[region]
    # all the segments in overlapped overlap, so each consecutive pairs of coordinates in sites should/could have different members
    sites = sorted(set(it.chain(*overlapped)))
    # iterate thru consecutive sites
    for start, stop in zip(sites, sites[1:]):
        # get the members of the regions that overlap the consecutive sites
        info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]
        # unpack the membership
        out_region[(start,stop)] = sorted(it.chain(*info))
    return out_region


def visualize_classifiers(classif_name, classif, ax):

    def plot_degree_classifier(classif, ax):

        labs = list(classif.classes_)

        XX, YY = np.meshgrid(np.linspace(0, 1, 500), np.linspace(0, 0.5, 250))

        Z = classif.predict(np.c_[XX.ravel(), YY.ravel()])

        # plot colors
        Z = np.array([labs.index(i) for i in Z])
        Z = Z.reshape(XX.shape)
        ax.pcolormesh(XX, YY, Z, cmap="rainbow")

        # plot countour lines
        Z = classif.predict_proba(np.c_[XX.ravel(), YY.ravel()])
        for i in range(len(labs)):
            Zi = Z[:,i].reshape(XX.shape)
            ax.annotate(labs[i], (XX[np.where(Zi>0.5)].mean(), YY[np.where(Zi>0.5)].mean()))
            ax.contour(XX, YY, Zi, [0.5], linewidths=4, colors="white")

        ax.set_xlabel("IBD1")
        ax.set_ylabel("IBD2")
        ax.set_title("degree classifier")

    def plot_hap_classifier(classif, ax):

        labs = list(classif.classes_)

        XX, YY = np.meshgrid(np.linspace(0.5, 1, 500), np.linspace(0.5, 1, 500))

        Z = classif.predict(np.c_[XX.ravel(), YY.ravel()])

        colorsList = [(198, 30, 49),(150, 130, 88),(38, 1, 90)]
        cmap = colors.ListedColormap(colorsList)

        Z = np.array([labs.index(i) for i in Z])
        Z = Z.reshape(XX.shape)
        ax.pcolormesh(XX, YY, Z, cmap=cmap)

        Z = classif.predict_proba(np.c_[XX.ravel(), YY.ravel()])
        for i in range(len(labs)):
            Zi = Z[:,i].reshape(XX.shape)
            ax.annotate(labs[i], (XX[np.where(Zi>0.5)].mean(), YY[np.where(Zi>0.5)].mean()))
            ax.contour(XX, YY, Zi, [0.5], linewidths=4, colors="white")
            
        ax.set_xlabel(r"$h_1$")
        ax.set_ylabel(r"$h_2$")
        ax.set_title("hap classifier")

    def plot_nsegs_classifier(classif, ax):

        X = np.arange(10, 120)

        labs = classif.classes_

        for index, lab in enumerate(labs):

            Y = classif.predict_proba([[0.25, i] for i in X])[:, index]

            ax.plot(X, Y, label=lab)

        ax.legend()
        ax.set_xlabel("Number of segments")
        ax.set_ylabel("Probability")
        ax.set_title("nsegs classifier")

    if classif_name == "degree":
        plot_degree_classifier(classif, ax)

    elif classif_name == "hap":
        plot_hap_classifier(classif, ax)

    elif classif_name == "nsegs":
        plot_nsegs_classifier(classif, ax)

def plot_prediction(df, classif, x, y):
    df["predicted"] = classif.predict(df[[x, y]].values)
    df["probability"] = [max(i) for i in classif.predict_proba(df[[x, y]].values)]

    fig, axs = plt.subplots(1, 2, figsize=(18, 8))

    sns.scatterplot(data=df, x=x, y=y, hue="probability", ax=axs[1])
    sns.scatterplot(data=df, x=x, y=y, hue="predicted", legend=False, alpha=0.4, ax=axs[0])

    for lab, tmp in df.groupby("predicted"):
        axs[0].text(x=tmp["ibd1"].mean(), y=tmp["ibd2"].mean(), s=lab, fontsize="medium")

    return fig, axs


# perform various computations on ibd segments for a pair of individuals
# takes as input a phasedibd segment data frame
class ProcessSegments:
    def __init__(self, pair_df):
        self.segs = pair_df

    '''Function takes as input a region_dict, which has the following format:
        {(start, stop): [obj1, obj2, obj3, ...]}
        Also takes as input a new_region which has format: [start, stop, obj]
        This function sees if there is overlap between the new region and existing regions
        If there is overlap, it splits the current region into new regions and concats the objs'''
    def split_regions(self, region_dict, new_region):
        # returns the overlap of 2 regions (<= 0 if no overlap)
        def overlap(region1, region2):
            start1, end1 = region1
            start2, end2 = region2
            return min(end1,end2) - max(start1,start2)

        # out region will be returned; is a dict of regions mapping to members of region    
        out_region = dict()
        # overlapped keeps track of all the regions that overlap with the new region
        overlapped = {tuple(new_region[:2]):[new_region[2]]}

        # iterate through the existing regions
        for region in sorted(region_dict):
            # if overlap
            if overlap(region, new_region[:2]) > 0:
                # the regions completely overlap, just add the member and return region dict
                if tuple(region) == tuple(new_region[:2]):
                    region_dict[region] += [new_region[2]]
                    return region_dict
                # bc the region overlaps, add it to overlapped
                overlapped[region] = region_dict[region]
            # no overlap, but add the region to the out_region dict
            else:
                out_region[region] = region_dict[region]
        
        # all the segments in overlapped overlap, so each consecutive pairs of coordinates in sites should/could have different members
        sites = sorted(set(it.chain(*overlapped)))
        # iterate thru consecutive sites
        for start, stop in zip(sites, sites[1:]):
            # get the members of the regions that overlap the consecutive sites
            info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]
            # unpack the membership
            out_region[(start,stop)] = sorted(it.chain(*info))
        
        return out_region

    # stitches together segments that are at most max_gap apart
    def segment_stitcher(self, segment_list, max_gap = 1):
        regions = {}

        # iterate through the segments
        for start, stop in segment_list:

            '''Alg works by adding start/stop positions to overlapped
                We init overlapped with the start/stop of the segment
                Next we iterate through the current regions
                If there is overlap, we add the start/stop to overlapped
                At the end, we create a new region taking the min of overlapped and the max of overlapped'''
            overlapped = {start, stop}

            # we rewrite the regions
            updated_regions = set()

            # iterate through the current regions
            for r1, r2 in regions:

                # if there is a overlap with the ibd segment
                if min(stop, r2) - max(start, r1) > -max_gap:
                    # add the start/stop of the region
                    overlapped |= {r1, r2}
                # no overlap, so add the region to the updated regions
                else:
                    updated_regions |= {(r1, r2)}

            # add the new segment/new region to updated regions
            updated_regions |= {(min(overlapped), max(overlapped))}

            # for the next iteration
            regions = updated_regions.copy()

        # return the regions
        return regions

    # returns ibd1, ibd2 values for the pair
    def get_ibd1_ibd2(self, n=False):

        # init with ibd1 of 0 cM and ibd2 of 0 cM
        ibd1, ibd2 = 0, 0
        n_ibd1, n_ibd2 = 0, 0

        # iterate through the chromosomes
        for chrom, chrom_df in self.segs.groupby("chromosome"):

            # rdata frame
            r = {}

            # iterate through the segments
            for _, row in chrom_df.iterrows():
                # add the segments from the perspective of id1; name hap index as 0 --> 1 and 1 --> 2
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id1_haplotype"]+1])
                # add the segments from the perspective of id2; rename the haplotype index as 0 --> 3 and 1 --> 4
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id2_haplotype"]+3])

            # iterate through the regions
            for (start, end), hap in r.items():
                # get the length of the region
                l = end - start
                # r is covered on all 4 haplotypes --> IBD2
                if sum(set(hap)) == 10:
                    ibd2 += l
                    n_ibd2 += 1
                # not ibd2
                else:
                    ibd1 += l
                    n_ibd1 += 1
        
        if n:
            return ibd1, n_ibd1, ibd2, n_ibd2
        
        return ibd1, ibd2

    # returns the number of IBD segments
    def get_n_segments(self):
        n = 0

        # run segment stitcher for each chromosome
        for _, chrom_df in self.segs.groupby("chromosome"):
            n += len(self.segment_stitcher(chrom_df[["start_cm", "end_cm"]].values))

        return n

    # returns the haplotype score of the pair
    def get_h_score(self, inter_phase=False):
        hap, tot = {0:0, 1:0}, 0
        chrom_tots = {0: [], 1: []}

        # iterate through the chromosome
        for _, chrom_df in self.segs.groupby("chromosome"):
            r= {}
            # iterate through segments
            for _, row in chrom_df.iterrows():
                # on id1
                r = self.split_regions(r, [row["start_cm"], row["end_cm"], row["id1_haplotype"]+1])
                # on id2
                r = self.split_regions(r, [row["start_cm"], row["end_cm"], row["id2_haplotype"]+3])

            # holds the hap information for the chromosome
            temp = {1:0, 2:0, 3:0, 4:0}
            for (start, end), hapl in r.items():
                # present on 1+ haplotype for at least one in the pair
                if len(hapl) > 2:
                    continue
                # get length of region and add to total
                l = end - start
                tot += l

                # add the hap information
                for h in hapl:
                    temp[h] += l
            # for id1
            hap[0] += max(temp[1], temp[2])
            # for id2
            hap[1] += max(temp[3], temp[4])

            chrom_tots[0].append([temp[1], temp[2]])
            chrom_tots[1].append([temp[3], temp[4]])

        t1 = sum([i[0] if inter_phase else max(i) for i in chrom_tots[0]])
        t2 = sum([i[0] if inter_phase else max(i) for i in chrom_tots[1]])

        h1 = tot and t1/tot or 0; h2 = tot and t2/tot or 0
        # return hap score
        # h1, h2 = tot and hap[0]/tot or 0, tot and hap[1]/tot or 0
        return h1 if h1 > 0.5 else 1-h1, h2 if h2 > 0.5 else 1-h2

    def ponderosa_data(self, genome_len, inter_phase, empty=False):
        # creates an empty class
        class ponderosa: pass
        ponderosa = ponderosa()

        if empty:
            return ponderosa

        # add ibd1, ibd2 data
        ibd1, n_ibd1, ibd2, n_ibd2 = self.get_ibd1_ibd2(n=True)
        ponderosa.ibd1 = ibd1 / genome_len
        ponderosa.ibd2 = ibd2 / genome_len
        ponderosa.n_ibd1 = n_ibd1
        ponderosa.n_ibd2 = n_ibd2

        # get the number of ibd segments
        ponderosa.n = self.get_n_segments()

        # get haplotype scores
        h1, h2 = self.get_h_score(inter_phase)
        ponderosa.h1 = h1
        ponderosa.h2 = h2

        # return the data
        return ponderosa

# pair_df is a dataframe of a pair of relatives
# mean_d is the mean distance between switch errors
def introduce_phase_error(pair_df, mean_d):
    
    # given a mean distance between switch error returns a list of randomly drawn sites
    def generate_switches(mean_d, index):
        #start the switch at 0
        switches, last_switch = [], 0
        
        # longest chrom is 287 cM 
        while last_switch < 300:
            
            # add the new site to the previous site
            switches.append(np.random.exponential(mean_d) + last_switch)
            
            # previous site is now the new site
            last_switch = switches[-1]
            
        # return
        return [(i, index) for i in switches]
    
    # store the newly create segments
    new_segments = []
    
    for chrom, chrom_df in pair_df.groupby("chromosome"):
    
        # generate the switch locations
        s1 = generate_switches(mean_d, 0)
        s2 = generate_switches(mean_d, 1)
        switches = np.array(sorted(s1 + s2))
        switch_index, switches = switches[:,1], switches[:,0]

        # old segments
        segments = chrom_df[["start_cm", "end_cm", "id1_haplotype", "id2_haplotype", "id1", "id2", "chromosome"]].values

        # iterate through the segments
        for start, stop, hap1, hap2, id1, id2, chrom in segments:

            # get number of switches before the segment
            n_dict = {0: len(np.where(np.logical_and(switches<start, switch_index==0))[0]),
                    1: len(np.where(np.logical_and(switches<start, switch_index==1))[0])}


            # get the index of switches within the segment
            b = np.where(np.logical_and(switches>=start, switches<=stop))[0]

            # iterate through the switches and the switch index in the segment
            for s, index in zip(switches[b], switch_index[b]):
                # add the broken segment as the current start --> s and add n, which is the number of preceding switches
                new_segments.append([chrom, id1, id2, hap1, hap2, start, s, n_dict[0], n_dict[1]])

                # new start
                start = s

                # increase the number of switches by 1 but only on the relevant switch
                n_dict[index] += 1

            # add the final segment
            new_segments.append([chrom, id1, id2, hap1, hap2, start, stop, n_dict[0], n_dict[1]])

    pair_df = pd.DataFrame(new_segments, columns = ["chromosome", "id1", "id2", "id1_haplotype", "id2_haplotype", "start_cm", "end_cm", "n1", "n2"])
    pair_df["l"] = pair_df["end_cm"] - pair_df["start_cm"]
    pair_df = pair_df[pair_df.l >= 2.5]
    pair_df["id1_haplotype"] = pair_df[["id1_haplotype", "n1"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)
    pair_df["id2_haplotype"] = pair_df[["id2_haplotype", "n2"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)

    return pair_df


'''
Takes as input a file name that has two columns which are the ordered edges of the hierarchical pedigree structure.
This has support both for holding pairwise data and probabilities for individual pairs
'''
class PedigreeHierarchy:
    
    def __init__(self, hier_file):

        def create_digraph(yaml_file):
            # open yaml file
            i = open(yaml_file, "r")
            y = yaml.safe_load(i)

            # get the edges
            edges = [[data.get("parent", data["degree"]), r] for r, data in y.items()]

            # init graph and add the nodes
            g = nx.DiGraph()
            g.add_nodes_from([(node, {"p": 0, "p_con": 0, "method": "None"}) for node in it.chain(*edges)])
            g.add_edges_from(edges)

            # get relationships that are roots and add a root to root all the subtrees
            roots = [nodes for nodes in g.nodes if g.in_degree(nodes)==0]
            g.add_node("relatives", p=1, p_con=1, method="None")
            g.add_edges_from([["relatives", nodes] for nodes in roots])

            return g
            
        self.hier = create_digraph(hier_file)

        self.init_nodes = set(self.hier.nodes)

        self.degree_nodes = list(set(self.hier.successors("relatives"))-{"relatives"})

        self.levels = {}

        # iterate through the degrees of relatedness
        for degree in self.hier.successors("relatives"):

            # get the subgraph
            degree_tree = self.hier.subgraph(set(nx.descendants(self.hier, degree)) | {degree})

            cur_nodes = [node for node in degree_tree if degree_tree.out_degree(node)==0]

            levels = []
            while len(cur_nodes) > 0:
                
                levels.append(cur_nodes)

                nxt_nodes = []
                parent_nodes = {next(self.hier.predecessors(node)) for node in cur_nodes} - {"relatives"}
                for nxt in parent_nodes:
                    add = True
                    for cur in cur_nodes:
                        try:
                            if nx.dijkstra_path_length(self.hier, nxt, cur) > 1:
                                add = False
                                break
                        except:
                            continue

                    if add:
                        nxt_nodes.append(nxt)

                cur_nodes = nxt_nodes[:]

            self.levels[degree] = levels
    

    ### This set of functions is for holding/managing/plotting the hierarchy for a pair of individuals
    ##################################################################################################

    def add_probs(self, node, method, **kwargs):
        # a list of probs has been provided, not a single node
        if type(node) == list:
            add_list = node[:]
        # information for a single node has been provided (along with other pos args); convert to a list    
        else:
            add_list = [[node, kwargs["p_con"]]]

        # add the new information
        for node, p in add_list:
            self.hier.nodes[node]["p_con"] = p
            self.hier.nodes[node]["method"] = method
        
    # after all probabilities have been added, cascade the probabilities down the tree
    def compute_probs(self):
        # iterate through each node in a breadthfirst order
        for parent in nx.bfs_tree(self.hier, "relatives"):

            # get the child nodes of the parent node
            children = nx.descendants_at_distance(self.hier, parent, 1)
            
            # sum the probability of the all the children for rescaling
            p = np.nansum([self.hier.nodes[node]["p_con"] for node in children])

            if p == 0:
                continue
            
            # for each child, rescale the conditional probability and compute the probability by multiply the conditional by the prob of the parent
            for node in children:
                self.hier.nodes[node]["p_con"] /= p
                self.hier.nodes[node]["p"] = self.hier.nodes[node]["p_con"]*self.hier.nodes[parent]["p"]

    def most_likely_among(self, nodes):
        return nodes[np.argmax([self.hier.nodes[node]["p"] for node in nodes])]
    
    def top2(self, nodes):
        ps = [self.hier.nodes[node]["p"] for node in nodes]
        return [j for _,j in heapq.nlargest(2, zip(ps, nodes))]
    
    # starting at the root, traverses the path of most probable relationships until it reaches a probability below min_p
    def most_probable(self, min_p):

        # get the most probable degree of relatedness
        degree_nodes = [(self.hier.nodes[node]["p"], node) for node in nx.descendants_at_distance(self.hier, "relatives", 1)]
        degree_nodes.sort(reverse=True, key=lambda x: -1 if np.isnan(x[0]) else x[0])

        degree_p, degree = degree_nodes[0]

        if degree_p < min_p:
            return degree, degree_p

        levels = self.levels[degree]
        
        # matrix holds the probabilities
        prob_matrix = np.zeros((len(levels), max([len(i) for i in levels])))

        # add the probabilities to the prob matrix
        for index, level in enumerate(levels):
            prob_matrix[index,:len(level)] += np.array([self.hier.nodes[node]["p"] for node in level])

        # get lowest level where there is a at least one relationship prob > min_p
        level_index = np.where(prob_matrix > min_p)[0][0]

        # at the lowest level, get the index of the rel with the highest prob
        node_index = np.argmax(prob_matrix, axis=1)[level_index]

        # return the relationship and the prob
        return levels[level_index][node_index], prob_matrix[level_index][node_index]

 
    # plots the hierarchy and the associated probabilities
    def plot_hierarchy(self, in_g, min_display_p=-1):

        # remove nodes where the prob is below a certain probability
        keep_nodes = [node for node, attr in in_g.nodes(data=True) if attr["p"] > min_display_p]

        g = in_g.subgraph(sorted(keep_nodes))

        # get the position of the nodes
        pos = graphviz_layout(g, prog='dot')

        # get max number of nodes per level
        level_size = Counter([y_coord for _,(_, y_coord) in pos.items()])
        max_width = sorted([num for _, num in level_size.items()], reverse=True)[0]

        # set the size of the plot according to the max width of the tree
        f, ax = plt.subplots(figsize=(1.8*max_width, 1.8*max_width))

        # draw the edges
        for i, j in g.edges:
            x1, y1 = pos[i]; x2, y2 = pos[j]
            ax.plot([x1, x2], [y1, y2], color="black", zorder=0)

        # colormap used for choosing the color of the nodes
        cmap = mpl.colormaps['autumn_r']._resample(8)

        # add the probabilities
        for node, coords in pos.items():
            pos[node] = list(coords) + [g.nodes[node]["p"], g.nodes[node]["p_con"], g.nodes[node]["method"]]

        # plot the nodes and their probabilities
        for lab, (x, y, p, p_con, method) in pos.items():
            ax.text(x, y, lab + f"\np={round(p, 3)}\n{method}\n{round(p_con, 3)}", color="black", va='center', ha='center', 
            bbox=dict(edgecolor=cmap(float(p)), facecolor='white', boxstyle='round'))

        # minor graph aspects
        ax.set_aspect(1)
        ax.axis('off')
        return ax

    # plots the probability tree for the degree
    def plot_degree(self, show_zero_p):

        g = self.hier.subgraph(nx.descendants_at_distance(self.hier, "relatives", 1) | {"relatives"})

        return self.plot_hierarchy(g, -1 if show_zero_p else 0.0)

    # plots the probability tree for second degree relatives
    def plot_second(self, show_zero_p):

        g = self.hier.subgraph(nx.descendants(self.hier, "2nd") | {"2nd"})

        return self.plot_hierarchy(g, -1 if show_zero_p else 0.0)

    # plots the entire hierarchy
    def plot_all(self, show_zero_p):

        return self.plot_hierarchy(self.hier, -1 if show_zero_p else 0.0)
    
    def save_plot(self, which_plot, outfile, show_zero_p=True):

        if which_plot == "degree":
            self.plot_degree(show_zero_p)
            plt.savefig(outfile, dpi=300)

        elif which_plot == "second":
            self.plot_second(show_zero_p)
            plt.savefig(outfile, dpi=300)

        elif which_plot == "all":
            self.plot_all(show_zero_p)
            plt.savefig(outfile, dpi=300)

    ### This set of functions is for holding pairs known relatives
    ##################################################################################################

    # adds a pair of relatives
    def add_pair(self, pair, rel, attrs):
        # add the edge
        self.hier.add_edge(rel, pair)

        # add the pair attributes
        for attr, attr_val in attrs.items():
            self.hier.nodes[pair][attr] = attr_val

    # add pairs from a list of individuals
    def add_pairs_from(self, pair_list):
        for pair, rel, attrs in pair_list:
            self.add_pair(pair, rel, attrs)

    # add/update an attribute
    def update_attr(self, pair, attr, val):
        self.hier.nodes[pair][attr] = val

    # add/update an attribute from a list
    def update_attr_from(self, pair_list):
        for pair, attr, val in pair_list:
            self.update_attr(pair, attr, val)

    # given a relationship, returns the the relative pairs under that relationship
    def get_pairs(self, node):
        return nx.descendants(self.hier, node) - self.init_nodes
    
    # given a parent node, returns a df with all the descendant nodes and all attributes they have as columns
    def get_pair_df(self, node):

        # get the pairs who are descendant pairs
        pair_list = list(self.get_pairs(node))

        # holds all columns and column values
        columns = {"pair": pair_list,
                   "degree": [nx.shortest_path(self.hier, "relatives", node)[1] for node in pair_list],
                   "rel": [next(self.hier.predecessors(node)) for node in pair_list]}

        # iterate through each person and add
        for index, pair in enumerate(pair_list):
            for attr, val in self.hier.nodes[pair].items():
                # attribute does not exist yet, so add it
                if attr not in columns:
                    columns[attr] = list(np.arange(len(pair_list)) * np.nan)

                # add the attribute
                columns[attr][index] = val
        
        # create the df
        out_df = pd.DataFrame(columns)

        # add a col with the requested node
        out_df["requested"] = node

        return out_df

    def get_relative_nodes(self, node, include=False):
        return (nx.descendants(self.hier, node) & self.init_nodes) | ({node} if include else set())

    # given a list of relationship nodes, returns all pairs under
    def get_nodes_from(self, node_list):
        return set(it.chain(*[list(self.get_pairs(node)) for node in node_list]))

    # given a list of nodes, returns a dataframe of all the attrs of the pairs
    def get_pair_df_from(self, node_list):
        return pd.concat([self.get_pair_df(node) for node in node_list]).reset_index(drop=True)
    

class Relationship:
    def __init__(self, data):
        p1 = sorted([[1] + path for path in data.get(1, [])], key=lambda x: x[1:])
        p2 = sorted([[2] + path for path in data.get(2, [])], key=lambda x: x[1:])

        # support for inter-generational relationships
        longest_path = max([len(path) for path in p1 + p2])
        matrix = np.array(sorted([path + [0 for _ in range(longest_path - len(path))] for path in p1 + p2]))

        # set attrs
        self.sex_specific = data["sex"]
        self.mat = matrix if self.sex_specific else matrix[:,1:]

    # returns boolean if it is the given relationship
    def is_relationship(self, mat):
        mat = mat if self.sex_specific else mat[:,1:]
        mat = np.array(sorted(mat.tolist()))

        # not the same shape == can't be the relationship
        if self.mat.shape != mat.shape:
            return False

        # returns True if the relationship matches
        return np.abs((self.mat - mat)).sum() == 0
        

class RelationshipCodes:
    def __init__(self, yaml_file):
        # open the code yaml
        i = open(yaml_file, "r")
        code_yaml = yaml.safe_load(i)

        # load each relationship, create the matrix
        self.codes = []
        for name, data in code_yaml.items():
            r = Relationship(data)
            self.codes.append([name, r])

        # sort such that the sex-specific relationships come first
        self.codes.sort(key=lambda x: int(x[1].sex_specific), reverse=True)

    # converts the path dictionary to a matrix
    def convert_to_matrix(self, path_dict):
        # get the paths
        p1 = sorted([[1] + [direction for _, direction in path] for path in path_dict[1]], key=lambda x: x[1:])
        p2 = sorted([[2] + [direction for _, direction in path] for path in path_dict[2]], key=lambda x: x[1:])

        # get the kinships
        k1 = sum([0.5**(len(path)-2) for path in p1])
        k2 = sum([0.5**(len(path)-2) for path in p2])

        # support for inter-generational relationships
        longest_path = max([len(path) for path in p1 + p2])
        mat = np.array([path + [0 for _ in range(longest_path - len(path))] for path in p1 + p2])

        # returns the matrix and expected ibd1, ibd2 values
        return mat, k1*(1-k2) + k2*(1-k1), k1*k2

    # returns the relationship
    def determine_relationship(self, path_dict):
        # get the matrix and expected ibd values
        mat, ibd1, ibd2 = self.convert_to_matrix(path_dict)

        # the pair are the same generation
        same_gen = mat[:,1:].sum() == 0

        # iterate through the relationships
        for name, robj in self.codes:
            # boolean if the relationship is true
            if robj.is_relationship(mat):
                return name, ibd1, ibd2, mat, same_gen

        ### haven't found a relationship, need to make sure that it's not a reversed code

        # get the first column of the matrix
        pcol = mat[:,:1]
        # reverse the direction of the rest of the matrix and flip along the horizontal
        tmp = np.flip(mat[:,1:]*-1)
        # add the parent column to the flipped matrix
        rv_mat = np.append(pcol, tmp, axis=1)
        # bound at least one possible relationship that it could be
        found = sum([robj.is_relationship(rv_mat) for _, robj in self.codes]) > 0

        # relationship is not found
        return "nan" if found else "unknown", ibd1, ibd2, mat, same_gen
    
class Classifier:
    def __init__(self, X, y, ids, name, lda=None):
        def equal_sample_size(X, y, count):
            out_X, out_y = [], []
            for lab in set(y):
                lab_index = np.where(y==lab)[0]
                X_to_keep = np.random.choice(lab_index, count, replace=False)
                out_X.append(X[X_to_keep,:])
                out_y += [lab]*count
            return np.concatenate(out_X), np.array(out_y)
        
        # supplied the lda directly
        if lda != None:
            self.lda = lda
            self.X = lda.means_
            self.y = lda.classes_
            self.train_ids = np.arange(len(self.y)).astype(str)
            self.name = name

        else:
            self.train_ids = np.array(["_".join(sorted([i,j])) for i,j in ids])

            # Check that there are enough training pairs
            min_count = np.inf
            for lab, count in Counter(y).items():
                min_count = count if count < min_count else min_count
                # if only 1, must remove from the training
                if count == 1:
                    print(f"Only 1 {lab} found. Removing from training set.")
                    self.train_ids = self.train_ids[np.where(y!=lab)]
                    X = X[np.where(y!=lab)[0],:]
                    y = y[np.where(y!=lab)]
                elif count < 5:
                    print(f"Only {count} {lab} found. Retaining in dataset.")

            X, y = equal_sample_size(X, y, min_count)

            self.lda = LinearDiscriminantAnalysis().fit(X, y)
            self.X = X; self.y = y
            self.name = name

    def loo(self, ids):
        proba_arr = []

        # iterate through the training pairs
        for id1 in ids:
            lda = LinearDiscriminantAnalysis().fit(self.X[np.where(self.train_ids!=id1)[0],:],
                                                    self.y[np.where(self.train_ids!=id1)[0]])
            # get proba
            proba = lda.predict_proba(self.X[np.where(self.train_ids==id1)[0]])
            proba_arr.append(list(proba[0]))
        
        return np.array(proba_arr)
    
    def train_label(self, id1, id2):
        str_id = "_".join(sorted([id1, id2]))
        if str_id in self.train_ids:
            return self.y[np.where(self.train_ids==str_id)[0][0]]
        return "NA"

    '''Takes as input an nxm matrix X (n is number of samples, m is number of features)
    and an array of length n, ids. Returns a nxp matrix (p is the number of classes)
    and an array of length p that contains the classes which correspond to the probabilities in each output row'''
    def predict_proba(self, X, ids):
        # make sure the ids are ordered
        ids = np.array(["_".join(sorted([i,j])) for i,j in ids])

        # index of the training data versus the test data
        train_idx = np.where(np.in1d(ids, self.train_ids))[0]
        test_idx = np.where(~np.in1d(ids, self.train_ids))[0]

        # the list of probs to concatenate
        to_concat = []

        # testing pairs are included
        if len(test_idx) > 0:
            test_proba = self.lda.predict_proba(X[test_idx, :])
            to_concat.append(test_proba)

        # training pairs to be assessed are included
        if len(train_idx) > 0:
            print(f"Trained {self.name} classifier ({len(train_idx)} samples).")
            t1 = time.time()
            # and the leave-one-out probabilities
            train_proba = self.loo(ids[train_idx])
            to_concat.append(train_proba)
            print(f"\tPerformed leave-one-out x-validation for {len(train_idx)} samples ({round(time.time()-t1, 2)} seconds).")

        # must concat the two
        if len(to_concat) > 1:
            proba = np.concatenate(to_concat)
        
        else:
            proba = to_concat[0]
        
        proba = proba[np.where(ids[np.concatenate((test_idx, train_idx))]==ids[:, np.newaxis])[1], :]
        
        return proba, self.lda.classes_
    
    '''See predict_proba for the input. Returns the most likely relationship'''
    def predict(self, X, ids):
        proba, classes = self.predict_proba(X, ids)

        return np.array([classes[np.argmax(prob)] for prob in proba])



class Pedigree:
    def __init__(self, **kwargs):

        po_list = kwargs.get("po_list", [])
        samples = kwargs.get("samples", None)

        # store the samples g
        self.samples = samples

        # a list of po-pairs has been supplied
        if len(po_list) > 0:
            po = nx.DiGraph()
            po.add_edges_from(po_list)

        # samples graph has been supplied
        elif len(self.samples.g.nodes) > 0:
            tmp = nx.DiGraph()
            tmp.add_edges_from(it.chain(*[[[data["mother"], node], [data["father"], node]] for node, data in self.samples.g.nodes(data=True)]))
            # add sex and optionally, other attrs
            nx.set_node_attributes(tmp, {node: {attr: data[attr] for attr in ["sex", "age"]} for node, data in self.samples.g.nodes(data=True)})
            po = tmp.subgraph(set(tmp.nodes) - {-1})

        self.po = po

        # must supply the yaml file
        self.R = RelationshipCodes(kwargs["pedigree_file"])

        # create the pedigree hierarchy
        self.hier = PedigreeHierarchy(kwargs["tree_file"])

        # keep track of the dummy id number
        self.n = 1

        self.logger = logging.getLogger("pedigree")

        log = logging.FileHandler("pedigree.log", "w")

        self.logger.addHandler(log)
        self.logger.setLevel(logging.INFO)

        self.logger.info(f"PONDEROSA pedigree mode\nStarted {time.strftime('%Y-%m-%d %H:%M')}\n")


    # for a given focal individual, finds all relationships
    def focal_relationships(self, focal):
        # recursive function that, given a focal individual returns all paths to relatives
        def get_paths(cur_relative, path_ids, path_dirs, paths):
            # init the next set of relatives
            next_set = []

            # past the first iteration, so we can get down nodes, but only down nodes that are not in the path
            if len(path_dirs) > 1:
                next_set += [(nxt_relative,-1) for nxt_relative in self.po.successors(cur_relative) if nxt_relative not in path_ids]

            # we're still moving up, so we can get the up nodes
            if path_dirs[-1] > 0:
                next_set += [(nxt_relative, 1) for nxt_relative in self.po.predecessors(cur_relative)]

            # we can't keep going; base case
            if len(next_set) == 0:
                paths.append((path_ids, path_dirs))
                return paths

            # iterate through the new set of relatives
            for nxt_relative, direction in next_set:
                paths = get_paths(nxt_relative, path_ids + [nxt_relative], path_dirs + [direction], paths)
            return paths

        # given the output of get_paths, creates all sub-paths
        def merge_paths(paths):
            # init the dict to store each relative pair and the paths along each parental lineage
            rel_pairs = {id2: {1: set(), 2: set()} for id2 in it.chain(*[path_ids for path_ids,_ in paths])}

            # iterate through the paths
            for path_ids, path_dirs in paths:
                # zip the path ids and the path directions
                path = [(id2, pdir) for id2, pdir in zip(path_ids, path_dirs)]

                # iterate through each person in the path
                for index in range(1, len(path)):
                    # get the id of the relative
                    id2 = path[index][0]
                    # get the subpath from the focal to the current id2
                    subpath = path[1:index+1]
                    # determine which parent they are related through
                    parent_sex = self.po.nodes[subpath[0][0]]["sex"]
                    # add to the rel pairs dictionary
                    rel_pairs[id2][int(parent_sex)] |= {tuple(path[1:index+1])}
            
            return rel_pairs

        # get all paths
        path_list = get_paths(focal, [focal], [1], [])

        # keep track of unknown relationships
        unknown_rels = []

        # iterate through each relative of the focal individual
        for id2, path_dict in merge_paths(path_list).items():
            if focal == id2:
                continue

            # get the relationship 
            rname, e_ibd1, e_ibd2, mat, same_gen = self.R.determine_relationship(path_dict)

            # don't want to add if rname is nan or they are same generation and id2 > id1
            if rname == "nan" or (same_gen and focal > id2):
                continue

            if rname == "unknown":
                unknown_rels.append(np.array2string(mat[:,1:]))
                continue

            edge_data = self.samples.g.get_edge_data(focal, id2, {})
            attrs = {attr: edge_data.get(attr, np.nan) for attr in ["k_ibd1", "k_ibd2", "h", "h_error", "n"]}
            
            # add attrs from the pedigree
            attrs["ibd1"] = e_ibd1; attrs["ibd2"] = e_ibd2; attrs["mat"] = mat

            self.hier.add_pair((focal, id2), rname, attrs)

        return unknown_rels
    
    '''
    Finds sets of full-siblings and makes sure that they have the same parents
    '''
    def resolve_sibships(self):
        '''
        Takes as input a list of indviduals who are full-siblings
        If they all have both parents listed --> do nothing
        If they all have only one parent --> add another parent
        If they all have no parents --> add two parents
        Returns a list of new edges to add
        TODO: add support for sex
        '''
        def add_parent(fam):
            # for each individual in the "family", make a set of the number of parents that can be added (0, 1, or 2)
            new_parents = {2-self.po.in_degree(node) for node in fam}
            if len(new_parents) > 1:
                return []
            new_nodes = []
            n_parents_to_add = new_parents.pop()
            if n_parents_to_add == 1:
                sex = self.po.nodes[next(self.po.predecessors(next(iter(fam))))]["sex"]
                new_nodes += [[f"Missing{self.n}", node, {"1": "2", "2": "1"}[sex]] for node in fam]
                self.n += 1
            elif n_parents_to_add == 2:
                new_nodes += [[f"Missing{self.n}", node, "1"] for node in fam]
                new_nodes += [[f"Missing{self.n+1}", node, "2"] for node in fam]
                self.n += 2
            return new_nodes
        
        # all pairs of individuals who share the same parent 
        sib_pairs = set(it.chain(*[[tuple(sorted([id1, id2])) for id1, id2 in it.combinations(self.po.successors(node), r=2)] for node in self.po.nodes]))
        sib_pairs |= {tuple(sorted([id1, id2])) for id1, id2 in self.samples.get_edges(lambda x: x["k_degree"]=="FS")}

        # get the data frame of the sib pairs
        sib_df = self.samples.to_dataframe(sib_pairs, include_edges=True)

        # train gaussian mixture model
        gmm = GaussianMixture(n_components=2,
                              means_init=[[0.5, 0], [0.5, 0.25]],
                              covariance_type="spherical").fit(sib_df[["k_ibd1", "k_ibd2"]].values.tolist())

        self.logger.info(f"\nTrained sibling GMM. Half-sibling means are ibd1={round(gmm.means_[0][0], 2)}, ibd2={round(gmm.means_[0][1], 2)}")
        self.logger.info(f"Full-sibling means are ibd1={round(gmm.means_[1][0], 2)}, ibd2={round(gmm.means_[1][1], 2)}\n")
        
        # predict whether FS or HS
        sib_df["predicted"] = gmm.predict(sib_df[["k_ibd1", "k_ibd2"]].values.tolist())

        # the label 1 corresponds to full-siblings
        fs_g = nx.Graph()
        fs_g.add_edges_from(sib_df[sib_df.predicted==1][["id1", "id2"]].values)

        # get a list of the new parent ids to add
        new_parents = list(it.chain(*[add_parent(fam) for fam in nx.connected_components(fs_g)]))

        # add new parents
        tmp = self.po.copy()
        tmp.add_edges_from([i[:2] for i in new_parents])
        nx.set_node_attributes(tmp, {parent: {"sex": sex, "age": np.nan} for parent, _, sex in new_parents})
        
        self.logger.info(f"Added {self.n} missing parents that are shared between full-siblings.\n")

        # reassign
        self.po = tmp

    # finds all relationships for nodes in the graph
    def find_all_relationships(self):

        self.resolve_sibships()

        unknown_rels = it.chain(*[self.focal_relationships(focal) for focal in self.po.nodes])

        self.logger.info("The following unknown relationships were found:")

        for unkr, n in Counter(unknown_rels).items():
            self.logger.info(f"{n} of the following were found:")
            self.logger.info(unkr + "\n")

        return self.hier

    
class TrainPonderosa:
    def __init__(self, real_data, **kwargs):
        
        # obj will store pair information
        self.pairs = PedigreeHierarchy()

        # stores args
        self.args = kwargs
        
        # stores if real data
        self.read_data = real_data

        # get the name of the population
        self.popln = kwargs.get("population", "pop1")
        
        # if using real data, there are additional kwargs that are required. These include
        # samples: Samples object
        # pedigree: Pedigree object
        if real_data:

            # get the objects
            samples = kwargs["samples"]
            pedigree = kwargs["pedigree"]
            genome_len = samples.genome_len

            # iterate through all the pairs
            for id1, id2 in pedigree.pairs.get_nodes("relatives"):

                # get the ibd data from the pair
                edge_data = samples.g.get_edge_data(id1, id2)

                # no IBD data; continue
                if edge_data == None:
                    continue

                # get the ibd segments
                pair_df = edge_data["segments"]

                # make a process segments obj
                s = ProcessSegments(pair_df)

                # get an empty ponderosa data obj and fill with the edge data 
                data = s.ponderosa_data(None, empty=True)
                data.ibd1 = edge_data["ibd1"]; data.ibd2 = edge_data["ibd2"]
                data.h1 = edge_data["h"][id1]; data.h2 = edge_data["h"][id2]
                data.n = edge_data["n"]

                # add the data
                rel = nx.predecessor(pedigree.pairs,(id1,id2))
                self.pairs.add_relative(rel, data)

                # if a 2nd degree relative
                if nx.predecessor(pedigree.pairs, rel) == "2nd":
                        pair_df_se = introduce_phase_error(pair_df, kwargs.get("mean_d", 30))
                        
                        # get h1 and h2
                        s = ProcessSegments(pair_df_se)
                        data = s.ponderosa_data(genome_len)
                        
                        # add to pedigree hierarchy
                        self.pairs.add_relative("PhaseError", data)
        
        else:
            # load the simulated segments
            segment_files = os.listdir("segments/")
            # only take the segment file for the population
            sim_df = pd.concat([pd.read_feather(f"segments/{i}") for i in segment_files if i.split("_")[1] == self.popln])

            # get the genome len
            pop = sim_df["pop"].values[0]
            
            genome_len = 0
            for chrom in range(1, 23):
                temp = pd.read_csv(f"{pop}/sim_chr{chrom}.map", delim_whitespace=True, header=None)
                genome_len += (temp[2].max() - temp[2].min())

            # iterate through the different relationships
            for rel, rel_df in sim_df.groupby("relative"):
                
                # iterate through the pairs
                for pair, pair_df in rel_df.groupby(["id1", "id2"]):
                    
                    # get the data needed for ponderosa
                    s = ProcessSegments(pair_df)
                    data = s.ponderosa_data(genome_len)
                    
                    # add the data to the pair obj
                    self.pairs.add_relative(rel, data)
                    
                    # get switch errors if second
                    if rel in ["MGP", "PGP", "MHS", "PHS", "AV"]:
                        
                        # get switch errors
                        pair_df_se = introduce_phase_error(pair_df, kwargs.get("mean_d", 30))
                        
                        # get h1 and h2
                        s = ProcessSegments(pair_df_se)
                        data = s.ponderosa_data(genome_len)
                        
                        # add to pedigree hierarchy
                        self.pairs.add_relative("PhaseError", data)
                
    # given a set of relative nodes (e.g., 2nd, HS, MHS, etc.), return dataframe with the attributes in attrs
    def get_dataframe(self, rel_nodes, attrs):
        pair_df = pd.DataFrame()
        
        for node in rel_nodes:
            
            df = pd.DataFrame(self.pairs.get_pairs(node), columns = ["pair"])
            df["node"] = node
            pair_df = pd.concat([pair_df, df])            
        
        for col in attrs:
            pair_df[col] = pair_df["pair"].apply(lambda x: getattr(x, col))
            
        return pair_df[attrs + ["node"]]
    
    # given a classifier, plots the probability space
    def plot_classifier(self, classif, classif_name, xlim, ylim, xlab, ylab, data, n_comp):
        
        cmaps = ["Blues", "Greens", "Purples", "Reds"]
        colors = ["blue", "green", "purple", "red"]
        
        # plot classifier
        x = np.linspace(xlim[0], xlim[1], 100)
        y = np.linspace(ylim[0], ylim[1], 100)
        X, Y = np.meshgrid(x, y)
        
        for cat in range(n_comp):
            Z = []
            for row in range(100):
                Z.append([p[cat] for p in classif.predict_proba([[x, y] for x, y in zip(X[row], Y[row])])])
            plt.contourf(X, Y, np.array(Z), cmap = cmaps[cat], levels=np.linspace(0.5, 1, 10))
        plt.colorbar().set_label(label = "Probability", size=12)
        
        # plot the points
        sns.scatterplot(data = data, x = xlab, y = ylab, hue = "node", alpha = 0.8)
        
        # save fig
        plt.savefig("{}_classifier.png".format(classif_name), dpi=500)
        
        plt.close()
    
    # creates the classifier for degrees of relatedness
    def degree_classifier(self, classifier = "lda"):
        
        training_data = self.get_dataframe(["FS", "2nd", "3rd", "4th"], ["ibd1", "ibd2"])
        
        if classifier == "gmm":
            
            classif = GaussianMixture(n_components=4, means_init = [[0.5, 0.25], [0.5, 0], [0.25, 0], [0.125, 0]],
                                         weights_init = np.full(4, 0.25),
                                         covariance_type = "tied").fit(training_data[["ibd1", "ibd2"]].values)
            
        elif classifier == "lda":
            classif = LinearDiscriminantAnalysis().fit(training_data[["ibd1", "ibd2"]].values, training_data["node"].values)
            
        self.plot_classifier(classif, "degree", (0,0.75), (0,0.4), "ibd1", "ibd2", training_data, 4)
        
        return classif
        
    # creates the classifier for the n of segments
    def nsegs_classifier(self):
        
        # get the training data
        training_df = self.get_dataframe(["AV", "MGP", "MHS", "PGP", "PHS"], ["n", "ibd1", "ibd2"])
        
        # convert ibd1 and ibd2 to a cM kinship
        training_df["k"] = training_df.apply(lambda x: x.ibd2 + x.ibd1/2, axis = 1)
        
        # fit the model
        lda = LinearDiscriminantAnalysis().fit(training_df[["n", "k"]].values, training_df["node"].values)
        
        # plot classifier
        mean_k = training_df["k"].mean()
        
        # compute probabilities
        probs = lda.predict_proba([[n, mean_k] for n in range(20, 100)])
        
        # create df and plot
        fig, (ax1, ax2) = plt.subplots(1, 2)
        df = pd.DataFrame(probs)
        df[5] = np.arange(20,100)
        df = df.melt(id_vars = [5], value_vars=[0,1,2,3,4])
        sns.lineplot(data = df, ax = ax1, x = 5, y = "value", hue = "variable", palette = "tab10")
        plt.xlabel("n of segments")
        plt.ylabel("Probability")
        sns.histplot(data = training_df, ax = ax2, x = "n", hue = "node", bins = 60)
        
        fig.set_size_inches(12, 5)

        
        plt.savefig("nsegs_classifier.png", dpi=500)
        plt.close()
        
        return lda
    
    # creates the classifier for the haplotype score classifier
    def hap_classifier(self):
        
        # get the training data
        training_data = self.get_dataframe(["GP/AV", "HS", "PhaseError"], ["h1", "h2"])
        
        # train the classifier
        lda = LinearDiscriminantAnalysis().fit(training_data[["h1", "h2"]].values, training_data["node"].values)
        
        self.plot_classifier(lda, "hap", (0.5, 1), (0.5, 1), "h1", "h2", training_data, 3)
        
        return lda
    
    def train_classifiers(self, directory = ""):

        ### hap classifier
        hap_classif = self.hap_classifier()
        f = open(f"{directory}hap_classifier_{self.popln}.pkl", "wb")
        pkl.dump(hap_classif, f)
        f.close()

        ### n segs classifier
        nsegs_classif = self.nsegs_classifier()
        f = open(f"{directory}nsegs_classifier_{self.popln}.pkl", "wb")
        pkl.dump(nsegs_classif, f)
        f.close()

        ### degree classifier
        degree_classif = self.degree_classifier()
        f = open(f"{directory}degree_classifier_{self.popln}.pkl", "wb")
        pkl.dump(degree_classif, f)
        f.close()
        

'''
This function is for filling in the genetic position of a MAP file. It takes as input
(1) the map file that needs its genetic position column filled (input_map=)
(2) a genetic map file (map_file=)
(3) if (2) is not plink MAP format, the 0-indexed index of the chromosome number, physical position, genetic coordinate (in that order); otherwise optional (columns=)
(4) a file containing two columns: chromosome number and physical position; the outputted MAP file will ONLY contain these (sites=)
'''
def interpolate_map(input_map, map_file, **kwargs):
    ### Prep the map file
    cols = kwargs.get("columns", [0, 2, 3])

    map_file_df = pd.read_csv(map_file, delim_whitespace=True, header=None)
    map_file_df = map_file_df.rename({index: col for index, col in zip(cols, ["CHROM", "CM", "MB"])}, axis=1)

    ### Open the in map file
    in_map_df = pd.read_csv(input_map, delim_whitespace=True, header=None, names=["CHROM", "rsID", "CM", "MB"])

    ### Open the optional file containing the sites to keep
    sites_file = kwargs.get("sites", "")

    # the file exists, merge with the inputted map_df
    if os.path.exists(sites_file):
        sites_df = pd.read_csv(sites_file, delim_whitespace=True, header=None, names=["CHROM", "MB"])
        sites_df["keep"] = True

        in_map_df = in_map_df.merge(sites_df, on=["CHROM", "MB"], how="outer")
        in_map_df["keep"] = in_map_df["keep"].fillna(False)

    # no keep sites file, so keep all the  sites
    else:
        in_map_df["keep"] = True

    out_maps = []
    for chrom, chrom_df in in_map_df.groupby("CHROM"):
        # get the df of the genetic map df
        map_df = map_file_df[map_file_df.CHROM==chrom]

        # make sure it's sorted correctly
        chrom_df = chrom_df.sort_values("MB")

        # linear interpolation
        chrom_df["CM"] = np.interp(chrom_df["MB"], map_df["MB"], map_df["CM"])

        out_maps.append(chrom_df)

    # concat all the chromosomes together
    out_map_df = pd.concat(out_maps).fillna("NA")

    # write it out
    out_map_df[out_map_df.keep==True][["CHROM", "rsID", "CM", "MB"]].to_csv(input_map.replace(".map", "_interpolated.map"), header=False, index=False, sep=" ")



def run_phasedibd(input_vcf, input_map, **kwargs):

    import phasedibd as ibd

    haps = ibd.VcfHaplotypeAlignment(input_vcf, input_map)
    tpbwt = ibd.TPBWTAnalysis()
    
    phase_corr_mode=kwargs.get("use_phase_correction", False)  # SDS edited
    
    ibd_results = tpbwt.compute_ibd(haps, 
                                    use_phase_correction=phase_corr_mode,  # SDS edited
                                    L_f=kwargs.get("L_f", 3.0))
    
    print("Phase correction mode: " + str(phase_corr_mode))  # SDS added
    
    # get a list of the samples in the vcf
    with open(input_vcf) as vcf:
        for line in vcf:
            if "#CHROM" in line:
                samples = line.split()[9:]
                break

    # rename the column with the corresponding VCF ID
    ibd_results["id1"] = ibd_results["id1"].apply(lambda x: samples[x])
    ibd_results["id2"] = ibd_results["id2"].apply(lambda x: samples[x])

    # if output file is blank, return the dataframe, otherwise write it
    if kwargs.get("output", "") != "":
        ibd_results.to_csv(kwargs["output"], index=False, sep="\t")
    
    else:
        return ibd_results





        
class RemoveRelateds:

    def __init__(self):
        self.seed = np.random.choice(np.arange(20000))
        np.random.seed = self.seed

    # takes as input a king file
    # threshold_func is a func that takes as input PropIBD, IBD1Seg, IBD2Seg, InfType and returns True if the input is considered to be related
    def king_graph(self, king_file, threshold_func, other_args=[]):

        if type(king_file) == str:
            # read in king file
            king = pd.read_csv(king_file, delim_whitespace = True, dtype = {"ID1": str, "ID2": str})

        # we can also load an existing king df
        else:
            king = king_file

        # create graph structure with the kinship coeff
        self.kinG = nx.Graph()
        self.kinG.add_weighted_edges_from(king[["ID1", "ID2", "PropIBD"]].values)

        # build the kinship graph
        G = nx.Graph()
        G.add_nodes_from(self.kinG.nodes)
        king_related = king[king[["PropIBD", "IBD1Seg", "IBD2Seg", "InfType"]].apply(lambda x: threshold_func(*x, *other_args), axis = 1)]
        G.add_edges_from(king_related[["ID1", "ID2"]].values)
        import pdb; pdb.set_trace()

        return G

    def unrelated_family(self, g):

        # keep track of the unrelated nodes
        unrelated_nodes = list()

        # dict maps an ID to the number of close relatives
        degree_d = dict(g.degree())

        # each iteration removes a node from degree_d so will be executed len(degree_d) times
        while len(degree_d) > 0:

            # create function that returns the num of close relatives +- random noise for tie-breakers
            randmin = lambda x: np.random.normal(degree_d[x] , 1.5)

            # picks the node with the fewest close relatives
            node1 = min(degree_d, key = randmin)

            # add the node to unrelated_nodes but only if not already related
            add = True
            for node2 in unrelated_nodes:
                if g.has_edge(node1, node2):
                    add = False
            if add:
                unrelated_nodes.append(node1)

            # delete the node from degree_d regardless of if it's added or not
            del degree_d[node1]

        return unrelated_nodes

    def get_unrelateds(self, G):
        # object to store various components of the run
        # n_comp is number of distinct families, unrelateds holds the unrelateds, max k holds the highest kinship value of the set
        run = type('run', (object,), {"n_comp": 0, "unrelateds": []})

        # iterate through each "family" (clusters of relatives entirely unrelated)
        for i in nx.connected_components(G):
            g = G.subgraph(i)
            run.unrelateds += self.unrelated_family(g)
            run.n_comp += 1

        # return the run object
        return run

    # run it over multiple iterations using different seeds to get more individuals included
    # target is the min len of unrelateds for the loop to stop
    # max_iter is the max num of iterations before it stops
    def multiple_runs(self, G, target, max_iter = 10):

        # keep track of the runs
        run_list = []

        # only run it max_iter times
        for i in range(max_iter):

            # choose and set new seed
            seed = np.random.choice(np.arange(20000))
            np.random.seed = seed

            # run the unrelateds algorithm
            run = self.get_unrelateds(G)
            run.seed = seed

            # add the new run to run_list
            run_list.append(run)

            # if the most recent run exceeds our target, stio
            if len(run.unrelateds) >= target:
                print(f"Target of {target} relatives found")
                break

            print(f"Running iteration {i+1}...found a set of {len(run.unrelateds)}")

        # sort by length and get run with longest unrelateds list
        run_list.sort(key = lambda x: len(x.unrelateds), reverse = True)
        run = run_list[0]

        # set the class seed
        self.seed = run.seed

        return run


    def write_out(self, run, prefix = "unrelateds"):

        # write log file
        log = open(f"{prefix}.log", "w")
        log.write(f"{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n")
        log.write(f"Random seed: {self.seed}\n\n")
        log.write(f"Found {run.n_comp} distinct families/family at the given kinship threshold\n")
        log.write(f"Wrote a total of {len(run.unrelateds)} relatives to {prefix}.txt\n")
        log.write(f"Of these, the max kinship (or kinship equivalent) is {run.max_k}\n")
        log.close()

        # write out unrelateds
        out = open(f"{prefix}.txt", "w")
        out.write("\n".join(run.unrelateds) + "\n")
        out.close()

    def get_unrelated_set(self, G, **kwargs):
        # run iterations to find the largest set
        run = self.multiple_runs(G,
                                 target = kwargs.get("target", np.inf),
                                 max_iter = kwargs.get("max_iter", 25))
        
        run.max_k = max([d["weight"] for _,_,d in self.kinG.subgraph(run.unrelateds).edges(data=True)])

        if kwargs.get("prefix", "") != "":
            self.write_out(run, kwargs["prefix"])

        return run

    def get_unrelated_set_lda(self, king_file, path_to_lda, max_degree, **kwargs):

        if os.path.exists(path_to_lda):
            i = open(path_to_lda, "rb")
            degree_lda = pkl.load(i)
        else:
            degree_lda = path_to_lda.lda

        degrees = ["4th", "4th+", "3rd", "3rd+", "2nd", "2nd+", "PO", "FS", "1st"]
        related_degree = degrees[degrees.index(max_degree)+1:]

        king_df = pd.read_csv(king_file, delim_whitespace=True)

        king_df["predicted"] = degree_lda.predict(king_df[["IBD1Seg", "IBD2Seg"]].values)

        self.kinG = nx.Graph()
        self.kinG.add_weighted_edges_from(king_df[["ID1","ID2","PropIBD"]].values)

        G = nx.Graph()
        G.add_edges_from(king_df[king_df.predicted.isin(related_degree)][["ID1","ID2"]].values)

        self.get_unrelated_set(G, prefix=kwargs.get("prefix", ""), max_iter=kwargs.get("max_iter", 25), target=kwargs.get("target", np.inf))

    def king_unrelateds(self, king_file: str, max_k: float):

        # create the relatedness network from the king file
        G = self.king_graph(king_file, lambda propIBD, a, IBD2Seg, c: propIBD > 0.1 or IBD2Seg > 0.03)

        # run 10 iterations to find the largest set
        run = self.multiple_runs(G, target = np.inf, max_iter = 10)

        # write out the file
        outfile = open("sim_keep.txt", "w")
        _ = outfile.write("\n".join(run.unrelateds))
        
        print(f"Wrote out {len(run.unrelateds)} IDs to keep to 'sim_keep.txt'")


def quick_kinship(pair_df):
    # keep track of cm ibd1, ibd2
    ibd1, ibd2 = 0, 0

    # iterate through the chromosomes
    for _, chrom_df in pair_df.groupby("chromosome"):
        # add a column that is the cm of the ibd segments
        chrom_df["r"] = chrom_df.apply(lambda x: np.arange(ceil(x.start_cm), floor(x.end_cm)+1), axis=1)
        chrom_df["sites1"] = chrom_df.apply(lambda x: [(i, x.id1_haplotype) for i in x.r], axis=1)
        chrom_df["sites2"] = chrom_df.apply(lambda x: [(i, x.id2_haplotype) for i in x.r], axis=1)

        sites1 = [i[0] for i in set(it.chain(*chrom_df.sites1.values))]
        sites2 = [i[0] for i in set(it.chain(*chrom_df.sites2.values))]

        sites = Counter(sites1 + sites2)

        ibd1 += sum([1 for site, c in sites.items() if c < 4])
        ibd2 += sum([1 for site, c in sites.items() if c == 4])

    return ibd1, ibd2      

    
class Karyogram:
    def __init__(self, map_file, cm = True):
        if type(map_file) != list:
            map_file = [map_file]

        df = pd.DataFrame()
        for mapf in map_file:
            temp = pd.read_csv(mapf, delim_whitespace=True, header = None)
            df = pd.concat([df, temp])

        self.chrom_ends = {}
        self.max_x = 0
        for chrom, chrom_df in df.groupby(0):
            self.chrom_ends[chrom] = (min(chrom_df[2 if cm else 3]), max(chrom_df[2 if cm else 3])-min(chrom_df[2]))
            self.max_x = self.max_x if sum(self.chrom_ends[chrom]) < self.max_x else sum(self.chrom_ends[chrom])

        self.chrom_y = {(chrom, hap): (chrom - 1)*9 + 4*hap for chrom, hap in it.product(np.arange(1, 23), [0, 1])}

    def plot_segments(self, segments, **kwargs):

        # init the figure
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 20)

        # add the chromosome templates
        for chrom, hap in it.product(np.arange(1, 23), [0, 1]):
            rect = patches.Rectangle((self.chrom_ends[chrom][0], self.chrom_y[(chrom, hap)]),
                                    self.chrom_ends[chrom][1], 3, edgecolor = "black",
                                    facecolor = "darkgrey" if hap == 0 else "grey")
            ax.add_patch(rect)

        # add the segments
        for chrom, start, stop, hap in segments:
            facecolor = kwargs.get("hap0_color", "cornflowerblue") if hap == 0 else kwargs.get("hap1_color", "tomato")
            rect = patches.Rectangle((start, self.chrom_y[(chrom, hap)]), stop - start, 3,
                                    edgecolor = "black", facecolor = facecolor, alpha = 0.8)
            ax.add_patch(rect)

        # re-label the y ticks
        ax.set_yticks([self.chrom_y[(chrom, 0)] + 3.5 for chrom in range(1, 23)])
        ax.set_yticklabels([str(chrom) for chrom in range(1, 23)])

        # set axes limits, remove spines, modify ticks
        plt.xlim(0, self.max_x)
        plt.ylim(-2, self.chrom_y[(22, 1)] + 10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', which='major', labelsize=16)
        plt.tick_params(left = False)

        plt.savefig(f"{kwargs.get('file_name', 'karyogram')}.png", dpi = kwargs.get('dpi', 500))

'''This is a function that is designed to resolve half-sibs from FS.
Given that a pair shares one parent, is it more likely that they also
share their other parent (and are FS) or that their other parent
are different. It takes as arguments:
pedigree: networkx DiGraph where edges indicate parent --> offspring
pair_data: networkx Graph where edges contain attrs ibd1 and ibd2
classifier: an sklearn classifier that will spit out 2nd or FS; optional in which case it will train a model'''
def SiblingClassifier(pedigree, pair_data, dummy_n, classifier=None):
    # keep track of the sibling pairs
    sibling_pairs = set()

    # keep track of parentless nodes
    parentless = set()

    # iterate through all nodes
    for parent in pedigree.nodes:
        # get all children
        children = sorted(nx.descendants_at_distance(pedigree, parent, 1))

        # add all pairs to sibling_pairs
        sibling_pairs |= set(it.combinations(children, r=2))

        # see if there are any predecessors
        if len(set(pedigree.predecessors(parent))) == 0:
            parentless |= {parent}

    # iterate through all edges between parentless nodes
    for id1, id2, data in pair_data.subgraph(parentless).edges(data=True):
        # only take putative FS
        if 0.10 < data["ibd2"] < 0.50:
            sibling_pairs |= {(id1, id2)}

    # create a df containing the pairs and ther ibd values
    sib_df = pd.DataFrame(list(sibling_pairs), columns=["id1", "id2"])
    sib_df["ibd1"] = sib_df.apply(lambda x: pair_data.get_edge_data(x.id1, x.id2, {"ibd1": np.nan})["ibd1"], axis=1)
    sib_df["ibd2"] = sib_df.apply(lambda x: pair_data.get_edge_data(x.id1, x.id2, {"ibd2": np.nan})["ibd2"], axis=1)
    sib_df = sib_df.dropna().reset_index(drop=True)

    if sib_df.shape[0] == 0:
        return dummy_n, pedigree.copy()

    # now get the parents each pair has in common
    sib_df["parents"] = sib_df.apply(lambda x: set(pedigree.predecessors(x.id1)) & set(pedigree.predecessors(x.id2)), axis=1)

    # they are FS if they have two parents in common and HS if one parent in common AND there is one other parent ruling out their HS
    sib_df["other_parents"] = sib_df.apply(lambda x: len((set(pedigree.predecessors(x.id1)) | set(pedigree.predecessors(x.id2))) - x.parents) > 0, axis=1)
    sib_df["half"] = sib_df.apply(lambda x: "FS" if len(x.parents) == 2 else ("2nd" if x.other_parents else np.nan), axis=1)

    ### Now we need to either load the classifier or train one
    if classifier == None:
        sibClassif = LinearDiscriminantAnalysis().fit(sib_df.dropna()[["ibd1", "ibd2"]].values, sib_df.dropna()["half"].values)

    # classifier is supplied; use it
    else:
        sibClassif = classifier

    # classify the NaN rows, which are putative FS
    putFS = sib_df[sib_df["half"].isna()].copy()
    putFS["half"] = sibClassif.predict(putFS[["ibd1", "ibd2"]].values)

    # get a graph of all parentless FS pairs; makes it easier to find clusters of FS
    FS = nx.Graph()
    FS.add_weighted_edges_from(putFS[putFS.half=="FS"][["id1", "id2", "parents"]].values)

    # make a copy of the pedigree obj
    ped_copy = pedigree.copy()

    # iterate through families of FS
    for fam in nx.connected_components(FS):

        # get the family subgraph
        sub = FS.subgraph(fam)
        # the set of the parents 
        parents = set(it.chain(*[list(e["weight"]) for _,_,e in sub.edges(data=True)]))

        # no parents; create both
        if len(parents) == 0:
            # parent names
            father = f"dummy{dummy_n+1}"
            mother = f"dummy{dummy_n+2}"
            dummy_n += 2

            # add the edges and the nodes attrs
            ped_copy.add_edges_from([[i, j] for i,j in it.product([father, mother], list(fam))])
            ped_copy.nodes[father]["sex"] = "1"; ped_copy.nodes[father]["age"] = np.nan
            ped_copy.nodes[mother]["sex"] = "2"; ped_copy.nodes[mother]["age"] = np.nan

        # one parent; create the other parent
        elif len(parents) == 1:
            # get the parent
            parent1 = parents.pop()
            # name the new parent
            parent2 = f"dummy{dummy_n+1}"; dummy_n += 1

            # get the sex of the parent and get the sex of the new parent
            sex1 = pedigree.nodes[parent]["sex"]
            sex2 = {"0": "0", "1": "2", "2": "1"}[sex1]

            # add the edges and the node attrs of the new parent
            ped_copy.add_edges_from([[i, j] for i,j in it.product([parent2], list(fam))])
            ped_copy.nodes[parent2]["sex"] = sex2; ped_copy.nodes[parent2]["age"] = np.nan

    return dummy_n, ped_copy

def remove_missing(vcffile):

    with open(vcffile) as file:

        for lines in file:

            if "#" in lines:
                print(lines.strip())

            else:
                lines = lines.replace(".|.", "0|0")
                lines = lines.replace(".|0", "0|0")
                lines = lines.replace("0|.", "0|0")
                lines = lines.replace(".|1", "1|1")
                lines = lines.replace("1|.", "1|1")
                print(lines.strip())
                

def parse_args():
    parser = argparse.ArgumentParser()

    # args for interpolation
    parser.add_argument("-interpolate", action="store_true")
    parser.add_argument("-genetic_map", type=str)
    parser.add_argument("-input_map", type=str)
    parser.add_argument("-columns", nargs='+', default=[0, 2, 3], type=int)
    parser.add_argument("-sites", type=str, default="")

    # args for running phasedibd
    parser.add_argument("-phasedibd", action="store_true")
    parser.add_argument("-input_vcf", type=str)
    parser.add_argument("-output", type=str)
    parser.add_argument("-use_phase_correction", action="store_true")  # SDS added

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse_args()


    if args.interpolate:
        interpolate_map(input_map=args.input_map, map_file=args.genetic_map, columns=args.columns, sites=args.sites)

    if args.phasedibd:
        run_phasedibd(input_vcf=args.input_vcf,
                      input_map=args.input_map,
                      output=args.output,
                      use_phase_correction=args.use_phase_correction)  # SDS added


'''
Archived code

p = PedigreeNetwork("Himba_allPO.fam")
p.get_paths()

p = PedSims("")
p.subset_vcf("plink_keep.txt", "../ponderosa/plink_data/Himba_shapeit.chr1.vcf")
if __name__ == "__main__":
    if sys.argv[-1] == "plink_unrelateds":
        r = RemoveRelateds()
        r.plink_unrelateds(sys.argv[-3], float(sys.argv[-2]))

    if sys.argv[-1] == "concat_rel":
        p = PedSims("")
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(p.run_ibd_iteration, sys.argv[-4], sys.argv[-3], sys.argv[-2], chrom) for chrom in range(1, 23)]
            ibd_results = pd.concat([f.result() for f in concurrent.futures.as_completed(results)]).reset_index(drop=True)
        try:
            df = pd.read_feather("simulated_segments.f")
            df = pd.concat([df, ibd_results]).reset_index(drop=True)
        except:
            df = ibd_results

        df.to_feather("simulated_segments.f")

    if sys.argv[-1] == "missing":
        remove_missing(sys.argv[-2])

    if sys.argv[-1] == "pop_ibd":
        p = PedSims("")
        # map_file, rel, iter, pop, chrom
        map_file, rel, iter, pop, chrom = [sys.argv[i] for i in range(-6, -1)]

        # run ibd
        ibd_results = p.run_pop_ibd_iteration(map_file, rel, iter, pop, chrom)
        
        ibd_results.reset_index(drop=True).to_csv(f"{pop}/{pop}_chr{chrom}_{rel}{iter}.txt", sep="\t", index=False)

    if sys.argv[-1] == "rename_pop_ibd":
        segment_file, pop, rel, sim_iter = [sys.argv[i] for i in range(-5, -1)]

        p = PedSims("")
        p.rename_ibd_segments(segment_file, pop, rel, sim_iter)

    if sys.argv[-1] == "interpolate":
        interpolate_map(sys.argv[-3], sys.argv[-2])
    if sys.argv[-1] == "subset_vcf":
        p = PedSims("")
        p.subset_vcf(sys.argv[-3], sys.argv[-2])
        

class PedSims:

    # takes as input the path where the all sim data is
    def __init__(self, path):
        self.path = path
        self.ibd = pd.DataFrame()
        self.log = "simulated_pedigrees.log"
        out = open(self.log, "w")
        _ = out.write(f"Pedigree simulations with ped-sim\n{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n\n")
        out.close()

    def write_log(self, msg):
        out = open(self.log, "a")
        _ = out.write(msg)
        out.close()          

    def run_ibd_iteration(self, map_file, rel, iter, chrom):
        # converts indices to the actual iid
        def return_ids(vcff):
            with open(vcff) as vcf:
                    for lines in vcf:
                        if "#CHROM" in lines:
                            return lines.split()[9:]

        # dict that stores the IDs for each relationship type
        relative_ids = {"av": {"AV": ("g2-b2-i1", "g3-b1-i1"),
                                "FS": ("g2-b1-i1", "g2-b2-i1"),
                                "CO": ("g3-b1-i1", "g3-b2-i1"),
                                "CORM": ("g3-b2-i1", "g4-b1-i1")},
                        "mhs": {"MHS": ("g2-b1-i1", "g2-b2-i1"),
                                "MHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "MHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "MHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "phs": {"PHS": ("g2-b1-i1", "g2-b2-i1"),
                                "PHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "PHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "PHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "mgp": {"MGP": ("g1-b1-i1", "g3-b1-i1"),
                                "MGGP": ("g1-b1-s1", "g4-b1-i1")},
                        "pgp": {"PGP": ("g1-b1-i1", "g3-b1-i1"),
                                "PGGP": ("g1-b1-s1", "g4-b1-i1")}}

        def rename_id(cur_id, rel_name, index):
            cur_id = cur_id.split("_")[0]
            sim_iter = int("".join([i for i in cur_id if i.isnumeric()]))
            return f"i{iter}{rel_name}{sim_iter}_{index}"

        # hap = ibd.VcfHaplotypeAlignment(f"sim_chr{chrom}.vcf", map_file.replace("chr1", f"chr{chrom}"))
        hap = ibd.VcfHaplotypeAlignment(f"sim_chr{chrom}.vcf", map_file.replace("chr1", f"chr{chrom}"))
        tpbwt = ibd.TPBWTAnalysis()
        ibd_results = tpbwt.compute_ibd(hap, use_phase_correction=False)

        # convert the index IDs to their actual IDs
        convert = {index:i.split()[0] for index,i in enumerate(return_ids(f"sim_chr{chrom}.vcf"))}
        ibd_results["id1"] = ibd_results.apply(lambda x: convert[x.id1],axis=1)
        ibd_results["id2"] = ibd_results.apply(lambda x: convert[x.id2],axis=1)

        # extract the relative segments
        relative_segs = pd.DataFrame()

        for relative, (id1, id2) in relative_ids[rel].items():
            df = ibd_results[ibd_results.apply(lambda x: id1 in x.id1 and id2 in x.id2 and x.id1.split("_")[0] == x.id2.split("_")[0], axis = 1)].copy()
            df["id1"] = df["id1"].apply(lambda x: rename_id(x, relative, 1))
            df["id2"] = df["id2"].apply(lambda x: rename_id(x, relative, 2))
            df["relative"] = relative
            relative_segs = pd.concat([relative_segs, df.drop(["start", "end"], axis=1)])

        return relative_segs

    def run_pop_ibd_iteration(self, map_file, rel, iter, pop, chrom):
        import phasedibd as ibd

        # converts indices to the actual iid
        def return_ids(vcff):
            with open(vcff) as vcf:
                    for lines in vcf:
                        if "#CHROM" in lines:
                            return lines.split()[9:]

        # dict that stores the IDs for each relationship type
        relative_ids = {"av": {"AV": ("g2-b2-i1", "g3-b1-i1"),
                                "FS": ("g2-b1-i1", "g2-b2-i1"),
                                "CO": ("g3-b1-i1", "g3-b2-i1"),
                                "CORM": ("g3-b2-i1", "g4-b1-i1")},
                        "mhs": {"MHS": ("g2-b1-i1", "g2-b2-i1"),
                                "MHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "MHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "MHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "phs": {"PHS": ("g2-b1-i1", "g2-b2-i1"),
                                "PHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "PHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "PHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "mgp": {"MGP": ("g1-b1-i1", "g3-b1-i1"),
                                "MGGP": ("g1-b1-s1", "g4-b1-i1")},
                        "pgp": {"PGP": ("g1-b1-i1", "g3-b1-i1"),
                                "PGGP": ("g1-b1-s1", "g4-b1-i1")}}

        def rename_id(cur_id, rel_name, index):
            cur_id = cur_id.split("_")[0]
            sim_iter = int("".join([i for i in cur_id if i.isnumeric()]))
            return f"i{iter}{rel_name}{sim_iter}_{index}"
        
        vcff = f"{pop}/{rel}{iter}_chr{chrom}.vcf"

        print(vcff)

        # hap = ibd.VcfHaplotypeAlignment(f"sim_chr{chrom}.vcf", map_file.replace("chr1", f"chr{chrom}"))
        hap = ibd.VcfHaplotypeAlignment(vcff, map_file)
        tpbwt = ibd.TPBWTAnalysis()
        ibd_results = tpbwt.compute_ibd(hap, use_phase_correction=False)

        # convert the index IDs to their actual IDs
        convert = {index:i.split()[0] for index,i in enumerate(return_ids(vcff))}
        ibd_results["id1"] = ibd_results.apply(lambda x: convert[x.id1],axis=1)
        ibd_results["id2"] = ibd_results.apply(lambda x: convert[x.id2],axis=1)

        # extract the relative segments
        relative_segs = pd.DataFrame()

        for relative, (id1, id2) in relative_ids[rel].items():
            df = ibd_results[ibd_results.apply(lambda x: id1 in x.id1 and id2 in x.id2 and x.id1.split("_")[0] == x.id2.split("_")[0], axis = 1)].copy()
            # df["id1"] = df["id1"].apply(lambda x: rename_id(x, relative, 1))
            # df["id2"] = df["id2"].apply(lambda x: rename_id(x, relative, 2))
            df["relative"] = relative
            relative_segs = pd.concat([relative_segs, df.drop(["start", "end"], axis=1)])

        return ibd_results

    def rename_ibd_segments(self, segment_file, pop, rel, sim_iter):
        relative_ids = {"av": {"AV": ("g2-b2-i1", "g3-b1-i1"),
                                "FS": ("g2-b1-i1", "g2-b2-i1"),
                                "CO": ("g3-b1-i1", "g3-b2-i1"),
                                "CORM": ("g3-b2-i1", "g4-b1-i1")},
                        "mhs": {"MHS": ("g2-b1-i1", "g2-b2-i1"),
                                "MHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "MHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "MHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "phs": {"PHS": ("g2-b1-i1", "g2-b2-i1"),
                                "PHCO": ("g3-b1-i1" , "g3-b2-i1"),
                                "PHAV1": ("g2-b1-i1", "g3-b2-i1"),
                                "PHAV2": ("g2-b2-i1", "g3-b1-i1")},
                        "mgp": {"MGP": ("g1-b1-i1", "g3-b1-i1"),
                                "MGGP": ("g1-b1-s1", "g4-b1-i1")},
                        "pgp": {"PGP": ("g1-b1-i1", "g3-b1-i1"),
                                "PGGP": ("g1-b1-s1", "g4-b1-i1")}}
        
        def rename_id(rel_name, sim_iter, iiter, id_index):
            return f"{rel_name}{id_index}_{pop}_i{iiter}_sim{sim_iter}"

        keep_ids = {j:i for i,j in relative_ids[rel].items()}

        segments = pd.concat([pd.read_csv(segment_file.replace("chr1", f"chr{chrom}"), delim_whitespace=True) for chrom in range(1,23)])

        # id1, id2 should be ordered correctly
        ordered = segments["id1"] < segments["id2"]
        to_order = segments[~ordered]
        to_order[["id1","id1_haplotype","id2","id2_haplotype"]] = to_order[["id2","id2_haplotype","id1","id1_haplotype"]]

        # concat the ordered, not ordered pairs
        segments = pd.concat([segments[ordered], to_order])
        segments["pair"] = segments.apply(lambda x: (x.id1.split("_")[-1], x.id2.split("_")[-1]), axis=1)
        segments["relative"] = segments["pair"].apply(lambda x: keep_ids.get(x, np.nan))
        segments = segments.dropna(subset=["relative"], axis=0)

        # get the different sim iters and make sure they're from the same simuation
        segments["iter1"] = segments["id1"].apply(lambda x: int("".join([i for i in x.split("_")[1] if i.isnumeric()])))
        segments["iter2"] = segments["id2"].apply(lambda x: int("".join([i for i in x.split("_")[1] if i.isnumeric()])))
        segments = segments[segments.iter1 == segments.iter2]

        # rename the segments
        # segments["id1"] = segments.apply(lambda x: rename_id(x.relative, sim_iter, x.iter1, 1), axis=1)
        # segments["id2"] = segments.apply(lambda x: rename_id(x.relative, sim_iter, x.iter2, 2), axis=1)

        # rename the segments (alt)
        segments["id1"] = [rename_id(r, sim_iter, i, 1) for r, i in segments[["relative", "iter1"]].values]
        segments["id2"] = [rename_id(r, sim_iter, i, 2) for r, i in segments[["relative", "iter2"]].values]

        # add the pop col
        segments["pop"] = pop

        segments = segments[["chromosome", "id1", "id1_haplotype", "id2", "id2_haplotype", "start_cm", "end_cm", "start_bp", "end_bp", "relative", "pop"]].reset_index(drop=True)

        # write out
        segments.to_feather(f"segments/{rel}_{pop}_sim{sim_iter}_segments.f")


    # takes as input a pair's segments and returns IBD1 and IBD2 cM
    def ibd1_ibd2(self, seg_df):
        ibd1, ibd2 = 0, 0
        for chrom, chrom_df in seg_df.groupby("chromosome"):
            r = {}
            for _, row in chrom_df.iterrows():
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id1_haplotype"]])
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id2_haplotype"]+2])
            for (start, end), hap in r.items():
                l = end - start
                if 0 in hap and 1 in hap and 2 in hap and 3 in hap:
                    ibd2 += l
                else:
                    ibd1 += l
        return ibd1, ibd2
    
    def load_ibd(self):
        # load ibd segments
        ibd = pd.read_feather("simulated_segments.f")

        # initialize ibd graph
        ibdG = nx.Graph()

        # add eches to ibd graph
        for (id1, id2), pair_df in ibd.groupby(["id1", "id2"]):
            ibd1, ibd2 = self.ibd1_ibd2(pair_df)
            rel = pair_df["relative"].values[0]
            ibdG.add_edge(id1, id2, ibd = pair_df, ibd1 = ibd1, ibd2 = ibd2, rel = rel)

        # store the graph
        self.ibd = ibdG

        return ibdG

    def segment_stitcher(self, segment_list, max_gap = 1):
        regions = {}
        for start, stop in segment_list:
            overlapped = {start, stop}
            updated_regions = set()
            for r1, r2 in regions:
                if min(stop, r2) - max(start, r1) > -max_gap:
                    overlapped |= {r1, r2}
                else:
                    updated_regions |= {(r1, r2)}
            updated_regions |= {(min(overlapped), max(overlapped))}
            regions = updated_regions
        return regions

    def pair_segments(self, id1, id2, pair_df, max_gap = 1):
        pair_data = type('pair_data', (object,), {"n_segs": 0, "tot_cov": 0, "tot_all": 0, "hap": {id1: 0, id2: 0}})
        for chrom, chrom_df in pair_df.groupby("chromosome"):
            regions = self.segment_stitcher(chrom_df[["start_cm", "end_cm"]].values, max_gap)
            pair_data.n_segs += len(regions)
            pair_data.tot_cov += sum([stop - start for start, stop in regions])
            pair_data.tot_all += chrom_df["l"].sum()
            temp = {id1:{0:0, 1:0}, id2:{0:0, 1:0}}
            for _, row in chrom_df.iterrows():
                temp[row["id1"]][row["id1_haplotype"]] += row["l"]
                temp[row["id2"]][row["id2_haplotype"]] += row["l"]
            pair_data.hap[id1] += max([temp[id1][0], temp[id1][1]])
            pair_data.hap[id2] += max([temp[id2][0], temp[id2][1]])
        return pair_data

    # pair_df is a dataframe of a pair of relatives
    # mean_d is the mean distance between switch errors
    def introduce_phase_error(self, pair_df, mean_d):
        
        # given a mean distance between switch error returns a list of randomly drawn sites
        def generate_switches(mean_d, index):
            #start the switch at 0
            switches, last_switch = [], 0
            
            # longest chrom is 287 cM 
            while last_switch < 300:
                
                # add the new site to the previous site
                switches.append(np.random.exponential(mean_d) + last_switch)
                
                # previous site is now the new site
                last_switch = switches[-1]
                
            # return
            return [(i, index) for i in switches]
        
        # store the newly create segments
        new_segments = []
        
        for chrom, chrom_df in pair_df.groupby("chromosome"):
        
            # generate the switch locations
            s1 = generate_switches(mean_d, 0)
            s2 = generate_switches(mean_d, 1)
            switches = np.array(sorted(s1 + s2))
            switch_index, switches = switches[:,1], switches[:,0]

            # old segments
            segments = chrom_df[["start_cm", "end_cm", "id1_haplotype", "id2_haplotype", "id1", "id2", "chromosome"]].values

            # iterate through the segments
            for start, stop, hap1, hap2, id1, id2, chrom in segments:

                # get number of switches before the segment
                n_dict = {0: len(np.where(np.logical_and(switches<start, switch_index==0))[0]),
                        1: len(np.where(np.logical_and(switches<start, switch_index==1))[0])}


                # get the index of switches within the segment
                b = np.where(np.logical_and(switches>=start, switches<=stop))[0]

                # iterate through the switches and the switch index in the segment
                for s, index in zip(switches[b], switch_index[b]):
                    # add the broken segment as the current start --> s and add n, which is the number of preceding switches
                    new_segments.append([chrom, id1, id2, hap1, hap2, start, s, n_dict[0], n_dict[1]])

                    # new start
                    start = s

                    # increase the number of switches by 1 but only on the relevant switch
                    n_dict[index] += 1

                # add the final segment
                new_segments.append([chrom, id1, id2, hap1, hap2, start, stop, n_dict[0], n_dict[1]])

        pair_df = pd.DataFrame(new_segments, columns = ["chromosome", "id1", "id2", "id1_haplotype", "id2_haplotype", "start_cm", "end_cm", "n1", "n2"])
        pair_df["l"] = pair_df["end_cm"] - pair_df["start_cm"]
        pair_df = pair_df[pair_df.l >= 2.5]
        pair_df["id1_haplotype"] = pair_df[["id1_haplotype", "n1"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)
        pair_df["id2_haplotype"] = pair_df[["id2_haplotype", "n2"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)

        return pair_df


    def train_H_classifier(self, training_data, error_data, title, covariance_type = "tied"):
        error_data["drel"] = "No inference"

        pair_df = pd.concat([error_data, training_data])

        # get haplotype scores
        pair_df["h1"] = (pair_df["id1_hap"] / pair_df["tot_cov"]).apply(lambda x: min(x, 1))
        pair_df["h2"] = (pair_df["id2_hap"] / pair_df["tot_cov"]).apply(lambda x: min(x, 1))

        # get the haplotype scores to train
        X = pair_df[["h1", "h2"]].values.tolist() + pair_df[["h2", "h1"]].values.tolist()

        # get the means to initialize the gmm
        hs_mean = pair_df[pair_df.drel == "HS"]["h1"].mean()
        gpav1_mean, gpav2_mean = pair_df[pair_df.drel != "HS"][["h1", "h2"]].mean()
        init_means = [[gpav1_mean, gpav2_mean], [gpav2_mean, gpav1_mean], [hs_mean, hs_mean], [0.7, 0.7]]

        # init and fit the gmm
        gmm = GaussianMixture(n_components = 4, means_init = init_means, covariance_type = "tied")
        gmm.fit(X)
        gmm.weights_ = [0.25, 0.25, 0.25, 0.25]

        # plot the classifier
        x = np.linspace(0.65, 1, 50)
        y = np.linspace(0.65, 1, 50)
        X, Y = np.meshgrid(x, y)
        for cat in range(4):
            Z = []
            for row in range(50):
                Z.append([p[cat] for p in gmm.predict_proba([[x, y] for x, y in zip(X[row], Y[row])])])
            plt.contourf(X, Y, np.array(Z), cmap = ["Greens", "Greens", "Blues", "Reds"][cat], levels=np.linspace(0.5, 1, 10))

        plt.colorbar().set_label(label = "Probability", size=12)

        pair_df["shuffle"] = np.random.binomial(1, 0.5, pair_df.shape[0])
        pair_df["temp_h1"] = pair_df.apply(lambda x: x.h1 if x.shuffle else x.h2, axis=1)
        pair_df["temp_h2"] = pair_df.apply(lambda x: x.h2 if x.shuffle else x.h1, axis=1)
        pair_df["Relationship"] = pair_df["drel"]

        ax = sns.scatterplot(data = pair_df, x = "temp_h1", y = "temp_h2", hue = "Relationship", alpha = 0.5, palette = ["palegreen", "skyblue", "lightsalmon"], hue_order = ["GP/AV", "HS", "No inference"])

        sns.move_legend(ax, "lower center",
        bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False, fontsize = 10)

        plt.xlabel("Haplotype score 1", size = 12)
        plt.ylabel("Haplotype score 2", size = 12)
        
        plt.savefig(f"{title}_gmm.png", dpi = 500)
        plt.close()

        with open(f"{title}_gmm.pkl", "wb") as f:
            pickle.dump(gmm, f)

    def train_N_classifier(self, training_data, title):
    
        # initialize and train the classifier
        lda = LinearDiscriminantAnalysis()
        lda.fit(training_data[["tot_all", "n_segs"]].values, training_data["rel"].values)
        
        # for plotting the probabilities, plot at fixed k, which is set to the mean k value here
        mean_k = training_data["tot_all"].mean()
        
        # over a range of num of segments from 22 to 81
        X = [[mean_k, n] for n in range(22, 81)]
        
        # dataframe of the different probabilities
        pred = pd.DataFrame(lda.predict_proba(X), columns = lda.classes_)
        pred["n_segs"] = np.arange(22, 81)
        
        # the prob of being MGP or PGP is the prob of being GP, eg
        pred["GP"] = pred["MGP"] + pred["PGP"]
        pred["HS"] = pred["MHS"] + pred["PHS"]
        
        # change shape of df
        pred = pred.melt("n_segs", var_name = "Relationship", value_name = "Probability")
        
        # plot each relationship
        sns.lineplot(data = pred[~pred.Relationship.isin(["GP", "HS"])], x = "n_segs", y = "Probability", hue = "Relationship",
                    hue_order = ["PGP", "MGP", "PHS", "MHS", "AV"], palette = ["firebrick", "lightcoral", "mediumblue", "cornflowerblue", "mediumseagreen"]) 
        
        # plot only the higher order relationships (GP, HS, AV)
        ax = sns.lineplot(data = pred[pred.Relationship.isin(["GP", "HS", "AV"])], x = "n_segs", y = "Probability", hue = "Relationship",
                    hue_order = ["GP", "HS", "AV"], palette = ["tomato", "dodgerblue", "mediumseagreen"], legend = False, linewidth = 4, style = "Relationship", dashes = {"GP": (3, 3), "AV": (3, 3), "HS": (3, 3)})
        
        # plot formatting code
        sns.move_legend(ax, "lower center",
        bbox_to_anchor=(.5, 1), ncol=5, title=None, frameon=False, fontsize = 10)
        plt.xlabel("Number of IBD segments shared", fontsize = 12)
        plt.ylabel("Probability", fontsize = 12)
        
        plt.savefig(f"{title}_lda.png", dpi = 500)
        plt.close()

        with open(f"{title}_lda.pkl", "wb") as f:
            pickle.dump(lda, f)

    def generate_training_data(self, training_nodes, mean_d):

        if mean_d == None:
            data = [[node, self.famG.nodes[node]["data"]] for node in training_nodes]

        else:
            data = []
            for id1, id2 in training_nodes:
                ibd = self.ibd.get_edge_data(id1, id2)["ibd"]
                ibd_phase_error = self.introduce_phase_error(ibd, mean_d)
                data.append([(id1, id2), self.pair_segments(id1, id2, ibd_phase_error)])
        df_data = []
        for node, info in data:
            rel = self.famG.nodes[node]["rel"]
            id1, id2 = node
            df_data.append([info.hap[id1], info.hap[id2], info.tot_cov, info.tot_all, info.n_segs, rel, "HS" if "HS" in rel else "GP/AV"])

        pair_df = pd.DataFrame(df_data, columns=["id1_hap", "id2_hap", "tot_cov", "tot_all", "n_segs", "rel", "drel"])

        return pair_df

    def analyze_ibd(self, max_couple_k, max_inlaw_k):
        # add the pairs
        self.add_fam()

        # iterate through each pair and store data
        for id1, id2 in self.famG.nodes():
            ibd = self.ibd.get_edge_data(id1, id2)
            pair_data = self.pair_segments(id1, id2, pd.DataFrame(columns=["chromosome"]) if ibd == None else ibd["ibd"])
            self.famG.nodes[(id1, id2)]["k"] = pair_data.tot_all
            self.famG.nodes[(id1, id2)]["data"] = pair_data

        # only keep training pairs whose founders are unrelated
        keep_pairs = []
        for node1 in self.root_nodes:
            add = True
            for node2 in self.famG.successors(node1):
                k = self.famG.nodes[node2]["k"]
                rel = self.famG.nodes[node2]["rel"]
                if (rel == "couple" and k > max_couple_k) or (rel == "inlaws" and k > max_inlaw_k):
                    add = False
                    break
            if add:
                keep_pairs.append(node1)

        training_data = self.generate_training_data(keep_pairs, None)

        # introduce phase error
        error_data = self.generate_training_data(keep_pairs, 75)

        # look across different phase errors
        phase_errors = {"Perfect phase": training_data, 75: error_data}
        for mean_d in [5, 25, 50, 100, 125]:
            phase_errors[mean_d] = self.generate_training_data(keep_pairs, mean_d)
        
        self.train_H_classifier(training_data, error_data, "test")
        self.train_N_classifier(training_data, "test")

        return phase_errors

    def get_items(self):
        return self.famG, self.ibd, self.root_nodes

    def run_phasedibd(self, args):

        import phasedibd as ibd

        vcf_files = str(subprocess.run(f"ls {self.pop}_sims/*_{self.iter}.vcf", shell=True, stdout=subprocess.PIPE).stdout, encoding='utf-8').split()

        ### bgzip and index the files
        for index, vcf in enumerate(vcf_files):

            subprocess.run(f"bcftools view {vcf} -O z -o {vcf}.gz", shell=True)
            subprocess.run(f"bcftools index {vcf}.gz", shell=True)
            subprocess.run(f"rm {vcf}", shell=True)
            vcf_files[index] = f"{vcf}.gz"


        # hold all ibd segments
        ibd_segments = pd.DataFrame()

        for chrom in range(1, 23):

            map_df = pd.read_csv(args.map.replace("chr1", f"chr{chrom}"), delim_whitespace=True,
                                 names=["chrom", "rs", "cm", "mb"],
                                 dtype={"chrom": int, "rs": str, "cm":float, "mb":int})

            for vcf in vcf_files:

                tmp_file = vcf.split(".vcf.gz")[0] + f"_chr{chrom}.vcf"

                subprocess.run(f"bcftools view {vcf} --regions {chrom} -o {tmp_file}", shell=True)

                pos = str(subprocess.run(f'bcftools query -f "%POS " {tmp_file}', shell=True, stdout=subprocess.PIPE).stdout, encoding='utf-8').split()

                tmp_map = pd.DataFrame([int(i) for i in pos], columns=["mb"])
                tmp_map["cm"] = np.interp(tmp_map["mb"].values, map_df["mb"].values, map_df["cm"].values)
                tmp_map = tmp_map.merge(map_df[["chrom", "rs", "mb"]], on="mb", how="left")

                # create the file
                tmp_map[["chrom", "rs", "cm", "mb"]].to_csv(f"{tmp_file.replace('.vcf', '.map')}", index=False, header=False, sep=" ")

                # run ibd
                haplotypes = ibd.VcfHaplotypeAlignment(tmp_file, tmp_file.replace(".vcf", ".map"))
                tpbwt = ibd.TPBWTAnalysis()
                ibd_results = tpbwt.compute_ibd(haplotypes, use_phase_correction=False)

                # get the samples from the header
                head = str(subprocess.run(f'bcftools query -l {tmp_file}', shell=True, stdout=subprocess.PIPE).stdout, encoding='utf-8').split()
                ibd_results["id1"] = ibd_results["id1"].apply(lambda x: head[x])
                ibd_results["id2"] = ibd_results["id2"].apply(lambda x: head[x])

                # add to existing
                ibd_segments = pd.concat([ibd_segments, ibd_results])

        ibd_segments.to_csv(f"{self.path}ibd_segments_iter{self.iter}.txt", sep="\t", header=True, index=False)

    def analyze_ibd(self, args):

        id_files = str(subprocess.run(f"ls {self.pop}_sims/*_{self.iter}.ids", shell=True, stdout=subprocess.PIPE).stdout, encoding='utf-8').split()

        for id_file in id_files:

            rel = id_file.split("/")[-1].split("_")[0]

            ids = [tuple(i.split()) for i in open(id_file).readlines()]

            ids_df = pd.DataFrame([i if len(i)==2 else (np.nan, i[0])  for i in ids], columns=["sim_id", "id1"])

            fam_size = self.rel_founders[rel]

            rm_family = []
            for index in range(0, ids_df.shape[0], fam_size):

                fam_df = ids_df.iloc[index:index+fam_size]

                ks = [self.ibd_g.get_edge_data(i, j, {}).get("PropIBD", 0) for i,j in it.combinations(fam_df["id1"].values, r=2)]

                if max(ks) > self.k:
                    fam_name = fam_df.dropna()["sim_id"].values[0].split("_")[0]
                    rm_family.append(fam_name)
                    print(f"Removed family {fam_name} in iteration {self.iter}.")

        ibd_segments = pd.read_csv(f"{self.path}ibd_segments_iter{self.iter}.txt", delim_whitespace=True)
        ibd_segments["fid1"] = ibd_segments["id1"].apply(lambda x: x.split("_")[0])
        ibd_segments["fid2"] = ibd_segments["id2"].apply(lambda x: x.split("_")[0])
        ibd_segments = ibd_segments[ibd_segments.apply(lambda x: x.fid1 == x.fid2 and x.fid1 not in rm_family, axis=1)]

        segs = ProcessSegments(pd.DataFrame())
        for pair, pair_df in ibd_segments.groupby(["id1", "id2"]):

            segs.segs = pair_df
            ibd1, ibd2 = segs.get_ibd1_ibd2()

            print(ibd1, ibd2)

'''