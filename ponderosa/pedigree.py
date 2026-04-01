import yaml
from pathlib import Path
import networkx as nx
import numpy as np
import itertools as it
from typing import Tuple, List, Set, Dict, Any
from sklearn.utils import _safe_indexing
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.mixture import GaussianMixture
import copy

from .data_loading import Individuals, Pairs
from .config import PonderosaConfig
from .pedigree_utils import get_parental_paths, format_pedigree_paths, PATERNAL, MATERNAL

_PACKAGE_DIR = Path(__file__).parent
DEFAULT_PEDIGREE_CODES = _PACKAGE_DIR / "config" / "pedigree_codes.yaml"

class PedigreeHierarchy(nx.DiGraph):
    def __init__(self, rel_dict: dict):
        super().__init__()  # Initialize the DiGraph

        edges = [[data.get("parent", data["degree"]), r] for r, data in rel_dict.items()]
        self.add_edges_from(edges)

        degree_nodes = [node for node in self.nodes if self.in_degree(node) == 0]
        self.add_edges_from([("relatives", node) for node in degree_nodes])

        # Add attributes
        for node, node_data in rel_dict.items():
            node_attrs = node_data.copy()
            node_attrs[1] = node_attrs.get(1, [])
            node_attrs[2] = node_attrs.get(2, [])
            self.nodes[node].update(node_attrs)

    @classmethod 
    def from_yaml(cls, yaml_file: str = DEFAULT_PEDIGREE_CODES):
        with open(yaml_file, "r") as f:
            rel_dict = yaml.safe_load(f)
        return cls(rel_dict)

    def to_dict(self):
        return dict(self.nodes(data=True))

    def copy(self):
        return copy.deepcopy(self)


class PedigreeRegistry:

    def __init__(self, hierarchy: PedigreeHierarchy):

        # Keep track of the pairs in each node category
        self.nodes = {node: [] for node in hierarchy.nodes}

        # Keep track of the pairs added and errors
        self.pairs = []
        self.problem_pairs = []

        # Store the hierarchy
        self.hier = hierarchy.copy()
        
        # Store cache
        self._cache = dict()

    @classmethod
    def from_dict(cls, rel_dict: dict):
        hierarchy = PedigreeHierarchy(rel_dict)
        return cls(hierarchy)

    @classmethod
    def from_yaml(cls, yaml_file: str):
        hierarchy = PedigreeHierarchy.from_yaml(yaml_file)
        return cls(hierarchy)

    def add_pair(self, id1: str, id2: str, rel: str):
        self.nodes[rel].append((id1, id2))

    def add_pairs_from(self, pair_list: List[Tuple[str, str, str]]):
        for (id1, id2, rel) in pair_list:
            pair = self._order_pair(id1, id2)
            if pair in self.pairs:
                self.problem_pairs.append(pair)
                continue
            
            self.pairs.append(pair)
            self.add_pair(id1, id2, rel)

    def _order_pair(self, id1: str, id2: str) -> Tuple[str, str]:
        if id1 < id2:
            return (id1, id2)
        return (id2, id1)

    def _get_pairs(self, rel: str, ordered: bool = False) -> List[Tuple[str, str]]:
        if ordered:
            return [self._order_pair(*pair) for pair in self.nodes[rel]]

        return self.nodes[rel]

    def _get_pairs_from(self, rel_list: List[str], ordered: bool = False) -> List[Tuple[str, str]]:
        pairs = [self._get_pairs(rel, ordered) for rel in rel_list]

        return np.array(list(it.chain(*pairs)))

    def _get_descendant_nodes(self, rel: str):
        if rel in self._cache:
            descendant_nodes = self._cache[rel]
        else:
            descendant_nodes = list(nx.descendants(self.hier, rel)) + [rel]
            self._cache[rel] = descendant_nodes

        return descendant_nodes


    def get_relatives(self, rel: str, ordered: bool = False) -> List[Tuple[str, str]]:
        descendant_nodes = self._get_descendant_nodes(rel)

        return self._get_pairs_from(descendant_nodes, ordered)

    def _format_return(self, pair_list: np.array, rel_list: np.array, output_style: str = "flatten"):

        if output_style == "zip":
            return np.array(list(pair_list)), np.array(list(rel_list))

        elif output_style == "flatten":

            return np.array([[*pair, rel] for pair, rel in zip(pair_list, rel_list)])

    def get_relatives_from(self, rel_list: List[str], ordered: bool = False, output_style: str = "flatten") -> List[Tuple[str, str]]:
        pairs = [self.get_relatives(rel, ordered) for rel in rel_list]

        rel_list = [[rel]*len(i) for rel, i in zip(rel_list, pairs)]

        return self._format_return(it.chain(*pairs), it.chain(*rel_list), output_style=output_style)

    def get_relationships(self, pair_list: list) -> List[str]:

        out_pairs = {}
        for nodes, pairs in self.nodes.items():
            for id1, id2 in pairs:
                out_pairs[(id1, id2)] = nodes
                out_pairs[(id2, id1)] = nodes

        return [out_pairs.get(tuple(pair), "unknown") for pair in pair_list]

class Relationship:

    def __init__(self, rel_dict: dict):

        p1, p2 = get_parental_paths(rel_dict)

        self.sex_specific = rel_dict["sex"]
        self.type3 = len(p1) > 0 and len(p2) > 0

        matrix = format_pedigree_paths(p1 + p2)

        self.mat = matrix if self.sex_specific else matrix[:,1:]

    # returns boolean if it is the given relationship
    def is_relationship(self, mat) -> bool:

        mat = format_pedigree_paths(mat)


        type3 = np.ptp(mat[:,0]) > 0 # Gets the range of the parents

        mat = mat if self.sex_specific else mat[:,1:]

        # not the same shape == can't be the relationship
        if self.mat.shape != mat.shape:
            return False

        # returns True if the relationship matches
        return np.array_equal(self.mat, mat) and (type3 == self.type3)
        

class PedigreeCodes:


    def __init__(self, node_data_dict: dict):

        # load each relationship, create the matrix
        self.codes = []
        for name, data in node_data_dict.items():
            if data.get("path", False) == False:
                continue
            r = Relationship(data)
            self.codes.append([name, r])

        # sort such that the sex-specific relationships come first
        self.codes.sort(key=lambda x: int(x[1].sex_specific), reverse=True)

        self.unknown_codes = {}

    @classmethod
    def from_hierarchy(cls, hierarchy: PedigreeHierarchy):

        data_dict = hierarchy.to_dict()

        return cls(data_dict)

    @classmethod
    def from_yaml(cls, yaml_file: str):
        hierarchy = PedigreeHierarchy.from_yaml(yaml_file)

        data_dict = hierarchy.to_dict()

        return cls(data_dict)


    def _compute_kinship(self, p1: list, p2: list) -> Tuple[float, float]:
        k1 = sum([0.5**(len(path)-2) for path in p1])
        k2 = sum([0.5**(len(path)-2) for path in p2])

        return k1*(1-k2) + k2*(1-k1), k1*k2

    def add_unknown_code(self, code_matrix: np.ndarray, ibd1: float, ibd2: float, pair: Tuple[str, str] = None):

        matrix_str = str(code_matrix)

        if matrix_str not in self.unknown_codes:
            self.unknown_codes[matrix_str] = [ibd1, ibd2, []]

        if pair:
            self.unknown_codes[matrix_str][2].append(pair)

    def get_unknown_codes(self):
        pass

    # returns the relationship
    def determine_relationship(self, path_dict):

        p1, p2 = get_parental_paths(path_dict)

        ibd1, ibd2 = self._compute_kinship(p1, p2)

        mat = format_pedigree_paths(p1 + p2)

        # the pair are the same generation
        same_gen = mat[:,1:].sum() == 0

        # in the genetically older individual
        if mat[:,1:].sum() < 0:
            return "nan", ibd1, ibd2, mat, same_gen

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


# TODO: look for FS with no parents
class PedigreeGraph(nx.DiGraph):

    def __init__(self, po_list: List[Tuple[str, str]], # Lists child --> parent
                       default_missing_parent: List[str],
                       sexes_list: List[Tuple[str, int]]):
        """ Creates the parent-offspring di-graph """

        super().__init__()

        assert type(po_list) is np.ndarray

        self.add_edges_from(po_list[:,[1,0]])

        for iid, sex in sexes_list:
            assert type(sex) is int or type(sex) is np.int64
            self.nodes[iid]["SEX"] = sex

        self.remove_nodes_from(default_missing_parent)

    @classmethod
    def load_from_individual_data(cls, individuals: Individuals):

        po_list = individuals.retrieve_data("FATHER", "MOTHER", output_style="expand")

        sexes_list = individuals.retrieve_data("SEX", output_style="expand")

        return cls(po_list, individuals.get_default_missing_parent(), sexes_list)

    def add_parents_from(self, edge_list: List[Tuple[str, str, dict]]):

        updated_edge_list = []

        for parent_iid, child_iid, *data in edge_list:

            if len(data) > 0: # Attribute has been added

                data = data[0]

                parent_sex = data[parent_iid]

                if self._parent_of_sex(child_iid, parent_sex) \
                    or parent_sex == 0: # Already has a parent of that sex or parent sex is missing
                    continue

                self.add_node(parent_iid, SEX=parent_sex)

            updated_edge_list.append([parent_iid, child_iid])

        self.add_edges_from(updated_edge_list)

    def _sort_pair(self, id1: str, id2: str) -> Tuple[str, str]:
        if id1 < id2:
            return (id1, id2)
        return (id2, id1)

    def _get_sex(self, iid: str) -> int:
        return self.nodes[iid]["SEX"]

    def _parent_of_sex(self, iid: str, sex: int) -> bool:
        assert sex in [1, 2]

        parents = self._get_parents(iid)

        for parent in parents:
            if self._get_sex(parent) == sex:
                return True
        return False

    def has_mother(self, iid: str) -> bool:
        return self._parent_of_sex(iid, 2)

    def has_father(self, iid: str) -> bool:
        return self._parent_of_sex(iid, 1)

    def _get_parents(self, iid: str) -> Set[str]:
        parents = self.predecessors(iid)

        return set(parents)


class RelationshipPathFinder:

    def __init__(self, pedigree: PedigreeGraph):

        self.pedigree = pedigree
        self.nodes = self.pedigree.nodes

    @classmethod
    def from_edge_list(cls, po_list: List[Tuple[str, str]],
                       default_missing_parent: List[str],
                       sexes_list: List[Tuple[str, int]]):

        pedigree = PedigreeGraph(po_list, default_missing_parent, sexes_list)

        return cls(pedigree)

    def _get_successors(self, id1: str):
        return self.pedigree.successors(id1)

    def _get_predecessors(self, id1: str):
        return self.pedigree.predecessors(id1)

    def _get_paths(self, cur_relative: str, path_ids: List[str], path_dirs: List[int], paths: List[Tuple[List[str], List[int]]]) -> None:
        # init the next set of relatives
        next_set = []

        # past the first iteration, so we can get down nodes, but only down nodes that are not in the path
        if len(path_dirs) > 1:
            next_set += [(nxt_relative,-1) for nxt_relative in self._get_successors(cur_relative) if nxt_relative not in path_ids]

        # we're still moving up, so we can get the up nodes
        if path_dirs[-1] > 0:
            next_set += [(nxt_relative, 1) for nxt_relative in self._get_predecessors(cur_relative)]

        # we can't keep going; base case
        if len(next_set) == 0:
            paths.append((path_ids, path_dirs))
            return paths

        # iterate through the new set of relatives
        for nxt_relative, direction in next_set:
            paths = self._get_paths(nxt_relative, path_ids + [nxt_relative], path_dirs + [direction], paths)
        
        return paths

    def _merge_paths(self, paths: List[Tuple[List[str], List[int]]]) -> Dict[str, Dict[int, List[Tuple]]]:
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
                parent_sex = self.nodes[subpath[0][0]]["SEX"]
                # add to the rel pairs dictionary
                rel_pairs[id2][int(parent_sex)] |= {tuple(path[1:index+1])}

        out_paths = {}
        for id2, paths in rel_pairs.items():
            tmp_path = {}
            tmp_path[1] = [[step[1] for step in path] for path in paths.get(1, set())]
            tmp_path[2] = [[step[1] for step in path] for path in paths.get(2, set())]
            out_paths[id2] = tmp_path

        return out_paths

    def _get_focal_paths(self, focal: str) -> dict:

        path_list = self._get_paths(focal, [focal], [1], [])

        return self._merge_paths(path_list)

    def find_focal_relationships(self, focal: str, codes: PedigreeCodes) -> List[Tuple[str, str, str]]:

        focal_relatives = self._get_focal_paths(focal)

        pairs = []

        for id2, path_dict in focal_relatives.items():

            if focal == id2:
                continue

            # get the relationship 
            rname, e_ibd1, e_ibd2, mat, same_gen = codes.determine_relationship(path_dict)

            # don't want to add if rname is nan or they are same generation and id2 > id1
            if rname == "nan" or (same_gen and focal > id2):
                continue

            if rname == "unknown":
                codes.add_unknown_code(mat, e_ibd1, e_ibd2, pair=(focal, id2))
                continue

            pairs.append([focal, id2, rname])

        return pairs

    def find_all_relationships(self, registry: PedigreeRegistry, codes: PedigreeCodes):

        for focal in self.pedigree.nodes:

            relative_pairs = self.find_focal_relationships(focal, codes)

            for id1, id2, rname in relative_pairs:

                registry.add_pair(id1, id2, rname)

        return registry

    


# TODO: more sophisticated _get_fs_sets
class Siblings:

    def __init__(self, sibling_pairs: np.ndarray,
                       sibling_labels: np.ndarray,
                       sibling_ibd_data: np.ndarray,
                       min_training_pairs: int = 5):

        self.HS = np.array([])
        self.FS = np.array([])
        self.predHS = np.array([])
        self.predFS = np.array([])

        # Mask out pairs with NaN IBD values
        sibling_ibd_data, sibling_pairs, sibling_labels = self._nan_mask_data(sibling_ibd_data, sibling_pairs, sibling_labels)

        if sibling_ibd_data.shape[0] == 0:
            return
            
        self.HS = sibling_pairs[sibling_labels == "HS"]
        self.FS = sibling_pairs[sibling_labels == "FS"]

        # Get boolean masks
        resolved_mask = sibling_labels != "unresolved"
        unresolved_mask = sibling_labels == "unresolved"

        X = sibling_ibd_data[resolved_mask]
        y = sibling_labels[resolved_mask]

        classif, mapping_func = self._train_sibling_classifier(sibling_ibd_data, sibling_labels, resolved_mask, min_training_pairs)

        X_unresolved = sibling_ibd_data[unresolved_mask]
        unresolved_pairs = sibling_pairs[unresolved_mask]

        if unresolved_pairs.shape[0] == 0:
            self.predHS = np.array([])
            self.predFS = np.array([])

        else:
            unresolved_pairs, pred_sibship = self._predict_sibship(classif, mapping_func, X_unresolved, unresolved_pairs)
            
            self.predHS = unresolved_pairs[pred_sibship == "HS"]
            self.predFS = unresolved_pairs[pred_sibship == "FS"]

            self.dummy = 1

        self.classif = classif

    @classmethod
    def from_pedigree_data(cls, pedigree: PedigreeGraph, pairs: Pairs, min_training_pairs: int = 5):

        sibling_pairs, sibling_labels = cls._find_siblings(pedigree)

        sibling_ibd_data = pairs.get_pair_data_from(sibling_pairs, "IBD1", "IBD2", output_style="flatten")

        return cls(sibling_pairs, sibling_labels, sibling_ibd_data, min_training_pairs)

    def pair_is(self, id1: str, id2: str, rel: str):
        pair1, pair2 = (id1, id2), (id2, id1)
        if rel == "HS":
            if pair1 in self.HS or pair2 in self.HS:
                return True
            if pair1 in self.predHS or pair2 in self.predHS:
                return True

    @staticmethod
    def _find_siblings(pedigree: PedigreeGraph) -> Tuple[List[str], List[str]]:
        sib_pairs = set()

        for iid in pedigree.nodes:
            children = pedigree.successors(iid) # All the children

            sib_pairs |= set(pedigree._sort_pair(*pair) for pair in it.combinations(children, r=2))

        pairs, sibling_labels = [], []

        for id1, id2 in sib_pairs:
            id1_parents = pedigree._get_parents(id1)
            id2_parents = pedigree._get_parents(id2)

            sibling_class = "unresolved"

            if len(id1_parents) == 2 and len(id1_parents & id2_parents) == 2:
                pairs.append((id1, id2))
                sibling_labels.append("FS")
            if (len(id1_parents) == 2 or len(id2_parents) == 2) and len(id1_parents & id2_parents) == 1:
                pairs.append((id1, id2))
                sibling_labels.append("HS")

        return np.array(pairs, dtype=object), np.array(sibling_labels, dtype=str)

    def _nan_mask_data(self, arr: np.array, *args):

        if arr.shape[0] == 0:
            return tuple([arr, *args])

        nan_mask = ~np.isnan(arr).any(axis=1)

        cleaned_arrs = [arr[nan_mask]]

        for arg in args:

            cleaned_arrs.append(_safe_indexing(arg, nan_mask))

        return tuple(cleaned_arrs)

    def _train_classif(self, X: np.ndarray, y: np.array = np.array([])):

        if y.shape[0] > 0:
            assert X.shape[0] > 0
            
            classif = LinearDiscriminantAnalysis().fit(X, y)

            mapping_func = lambda x: x
        
        else:
            classif = GaussianMixture(n_components=2,
                              means_init=[[0.5, 0], [0.5, 0.25]],
                              covariance_type="spherical").fit(X)

            mapping_func = lambda x: "HS" if x == 0 else "FS"
        
        return classif, mapping_func


    def _train_sibling_classifier(self, X: np.ndarray, y: np.array, resolved_mask: np.array, min_training_pairs: int = 5):

        # No training pairs
        if resolved_mask.any() == False:

            return self._train_classif(X)

        X_clean = X[resolved_mask]
        y_clean = y[resolved_mask]

        unique_classes, class_counts = np.unique(y_clean, return_counts=True)

        if min(class_counts) >= min_training_pairs and len(unique_classes) == 2:

            return self._train_classif(X_clean, y_clean)

        return self._train_classif(X)


    def _predict_sibship(self, classif, mapping_func, X: np.array, pairs: np.array):

        X_clean, pairs_clean = self._nan_mask_data(X, pairs)

        y_pred = classif.predict(X_clean)

        return pairs_clean, np.array([mapping_func(i) for i in y_pred])

    def _get_fs_sets(self) -> List[Set[str]]:

        assert hasattr(self, "predFS")

        fs_sets = []

        for (id1, id2) in self.predFS:
            added = False
            for fs_set in fs_sets:
                if id1 in fs_set or id2 in fs_set:
                    fs_set |= {id1, id2}
                    added = True
                    break
            if not added:
                fs_sets.append({id1, id2})

        return fs_sets

    def _get_dummy_id(self):
        self.dummy += 1
        return f"Missing{self.dummy-1}"

    def get_dummy_parent_edges(self) -> List[Tuple[str, str, dict]]:
        def new_edge(parent, child, parent_sex):
            return [parent, child, {parent: {"SEX": parent_sex}}]

        dummy = 1

        edge_list = []
        for fs_set in self._get_fs_sets():
            for sex in [1, 2]:
                parent_iid = self._get_dummy_id()
                for iid in fs_set:
                    edge_list.append(new_edge(parent_iid, iid, sex))

        return edge_list


def build_pedigree(individuals: Individuals, pairs: Pairs, hierarchy: PedigreeHierarchy) -> Tuple[PedigreeGraph, PedigreeRegistry]:

    # Initialize the pedigree options
    pedigree = PedigreeGraph.load_from_individual_data(individuals)

    # Finds siblings and resolves ambuiguous relationships
    siblings = Siblings.from_pedigree_data(pedigree, pairs)

    # From the resolved siblings, get the edges for the missing parents
    updated_po_edges = siblings.get_dummy_parent_edges()
    # Add missing edges
    pedigree.add_edges_from(updated_po_edges)

    # Initialize the relationship registry to keep track of relative pairs
    registry = PedigreeRegistry(hierarchy)
    
    # Add relative pairs to the registry
    codes = PedigreeCodes.from_hierarchy(hierarchy)

    # Find all the relationships
    path_finder = RelationshipPathFinder(pedigree)
    registry = path_finder.find_all_relationships(registry, codes)

    return registry