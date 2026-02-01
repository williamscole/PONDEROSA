import numpy as np
import networkx as nx
import pandas as pd

from .pedigree import PedigreeHierarchy

def process_phase_error(X: np.ndarray, prob_threshold: float = 0.5) -> np.ndarray:

    assert np.all(np.isclose(X.sum(axis=1),1))

    output = np.full((X.shape[0], 2), np.nan)

    output[:,0] = np.where(X[:,2] < prob_threshold, X[:,0], output[:,0])
    output[:,1] = np.where(X[:,2] < prob_threshold, X[:,1], output[:,1])

    return output / output.sum(axis=1)[:, np.newaxis]


class MatrixHierarchy:
    """
    Docstring for MatrixHierarchy
    
    :var Args: Description
    :var sample_idx: Description
    :vartype sample_idx: Index
    :var Returns: Description
    """
    def __init__(self, directed_edge_list: list, index_to_pair: dict, methods: list):

        self.index_to_pair = index_to_pair
        self.n_pairs = len(index_to_pair)

        diG = nx.DiGraph()
        diG.add_edges_from(directed_edge_list)

        self.nodes = list(diG.nodes)

        self.node_to_idx = {node: idx for idx, node in enumerate(self.nodes)}

        self.method_to_idx = {method: idx for idx, method in enumerate(list(methods) + ["aggregated"])}

        # Creates a matrix that is mxnx3, whre m=no. of pairs, n=no. of relationships, 3 corresponds to posterior, conditional, LDA method
        self.matrix = np.full((self.n_pairs, len(self.nodes), 3), np.nan)
        self.relatives_idx = self.node_to_idx["relatives"]
        self.matrix[:, self.relatives_idx, 0] = 1
        self.matrix[:, self.relatives_idx, 1] = 1

        D = dict()
        for child in nx.bfs_tree(diG, "relatives"):

            idx = self.node_to_idx[child]

            D[idx] = []

            parent = list(diG.predecessors(child))

            if len(parent) > 0:

                pidx = self.node_to_idx[parent[0]]

                D[pidx].append(self.node_to_idx[child])

        self.levels = {i: np.array(j) for i, j in D.items()}

        self.g = diG

    @classmethod
    def from_hierarchy(cls, hierarchy: PedigreeHierarchy, index_to_pair: dict, methods: list):
        edge_list = list(hierarchy.edges())
        return cls(edge_list, index_to_pair, methods)

    def numeric_array(self, input_arr, nodes: bool, asfloat: bool):
        if type(input_arr) == str:
            output = self.node_to_idx[input_arr] if nodes else self.method_to_idx[input_arr]
            return np.float32(output) if asfloat else output
        if nodes:
            output_arr = np.array([self.node_to_idx[i] for i in input_arr])
        else:
            output_arr = np.array([self.method_to_idx[i] for i in input_arr])
        return output_arr.astype(np.float32) if asfloat else output_arr

    # Sets an attribute of the indices of the pairs whose probability of the specified node is >p
    def set_subset(self, subset_name, p, posterior=False):
        assert subset_name in self.nodes

        pidx = 0 if posterior else 1

        subset_idx = self.node_to_idx[subset_name]

        subset = np.where(self.matrix[:,subset_idx,pidx] > p)[0]

        setattr(self, subset_name, subset)

    def complex_subset(self, subset_name, index_arr):

        setattr(self, subset_name, index_arr)

    # Add the probabilities
    def add_probs(self, prob_arr, labels, method, subset_name=None):
        label_arr = self.numeric_array(labels, nodes=True, asfloat=False)
        method_idx = float(self.method_to_idx[method])

        if subset_name:
            subset_arr = getattr(self, subset_name)
            self.matrix[subset_arr[:, np.newaxis], label_arr, 1] = prob_arr
            self.matrix[subset_arr[:, np.newaxis], label_arr, 2] = method_idx

        else:
            self.matrix[:,label_arr,1] = prob_arr
            self.matrix[:,label_arr,2] = method_idx

    def add_priors(self):
        pass
        

    def _fill_nan(self):

        # Iterate through from leaf node to root
        for node in reversed(list(nx.topological_sort(self.g))):

            children_str = list(nx.descendants_at_distance(self.g, node, 1))

            if not children_str:
                continue

            node_idx = self.node_to_idx[node]


            children_idx = self.numeric_array(children_str, nodes=True, asfloat=False)

            needs_filling = np.isnan(self.matrix[:, node_idx, 1])

            self.matrix[:, node_idx, 2] = np.where(
                                                needs_filling,
                                                self.method_to_idx["aggregated"],
                                                self.matrix[:, node_idx, 2]
                                            )

            # Set the prob to the sum of the children IF empty
            self.matrix[:, node_idx, 1] = np.where(
                                                needs_filling,
                                                self.matrix[:, children_idx, 1].sum(axis=1),
                                                self.matrix[:, node_idx, 1]
                                            )

    def compute_probs(self):

        self._fill_nan()

        root_node = True
        for parent in nx.bfs_tree(self.g, "relatives"):
            children_str = nx.descendants_at_distance(self.g, parent, 1)
            if len(children_str) == 0:
                continue

            children = self.numeric_array(children_str, nodes=True, asfloat=False)
            
            # Normalize conditional probabilities
            p_arr = np.nansum(self.matrix[:, children, 1], axis=1)[:, np.newaxis] 
            self.matrix[:,children,1] = self.matrix[:,children,1] / p_arr
            
            # Get parent posterior probability
            if root_node:
                parent_p = np.ones(self.n_pairs)[:, np.newaxis]
                root_node = False
            else:
                parent_idx = self.node_to_idx[parent]
                parent_p = self.matrix[:,parent_idx,0][:, np.newaxis]

            # Update posterior probabilities
            self.matrix[:,children,0] = self.matrix[:,children,1] * parent_p

            
    def _operate_on(self, nodes, pidx, subset_name=None):
        node_idx = self.numeric_array(nodes, nodes=True, asfloat=False)

        if subset_name:
            subset_arr = getattr(self, subset_name)
            return self.matrix[subset_arr[:, np.newaxis], node_idx, pidx]
        
        return self.matrix[:,node_idx,pidx]

    def _posterior_proba(self, node_arr, subset_name=None):

        node_idx = self.numeric_array(node_arr, nodes=True, asfloat=False)

        if subset_name:
            subset_arr = getattr(self, subset_name)
            return self.matrix[subset_arr,:][np.arange(len(node_idx)), node_idx, 0]

        return self.matrix[np.arange(len(node_idx)), node_idx, 0]

    def _predict_proba(self, nodes, subset_name=None):

        return self._operate_on(nodes, 0, subset_name)

    def most_likely_among(self, nodes, asint, subset_name=None):

        mat = self._operate_on(nodes, 0, subset_name)

        most_likely = np.nanargmax(mat, axis=1)

        arr = np.array([nodes[i] for i in most_likely])

        if asint:
            return np.array([self.node_to_idx[i] for i in arr])
        
        return arr

    def top2(self, nodes, asint, subset_name=None):

        mat = self._operate_on(nodes, 0, subset_name)

        t2 = np.argpartition(mat, -2, axis=1)[:,-2:]

        arr = [(nodes[i], nodes[j]) for i, j in t2]

        if asint:
            return np.array([[self.node_to_idx[i], self.node_to_idx[j]] for i,j in arr])

        return arr

    def _most_probable_helper(self, arr):

        p = 1; new_node = 0

        while p >= self.min_p:

            node = new_node

            children = self.levels[node]

            if children.shape[0] == 0:
                break

            new_node = children[np.nanargmax(arr[children])]

            p = arr[new_node]

        return node

    def most_probable(self, min_p, asint=False):

        self.min_p = min_p

        node_arr = np.apply_along_axis(self._most_probable_helper, axis=1, arr=self.matrix[:,:,0])

        prob_arr = self.matrix[np.arange(self.n_pairs), node_arr, 0]

        if asint:
            return node_arr, prob_arr
        
        node_arr = [self.nodes[i] for i in node_arr]

        return node_arr, prob_arr

    def most_likely_degree(self, asint=False):

        degree_idx = self.levels[self.relatives_idx]

        most_likely_idx = np.nanargmax(self.matrix[:,degree_idx,0], axis=1)

        if asint:
            return most_likely_idx

        return [self.nodes[degree_idx[i]] for i in most_likely_idx]

    def to_dataframe(self, min_p: float = 0.50) -> pd.DataFrame:

        id_df = pd.DataFrame(index=range(self.n_pairs))

        pairs = np.vstack(id_df.index.map(self.index_to_pair))

        id_df["id1"] = pairs[:,0]
        id_df["id2"] = pairs[:,1]
        id_df["degree"] = self.most_likely_degree()

        most_probable, prob = self.most_probable(min_p)
        id_df["pred_rel"] = most_probable
        id_df["prob"] = prob

        df = pd.DataFrame(self.matrix[:,:,0]).replace(np.nan, 0)

        df.rename({j:i for i,j in self.node_to_idx.items()}, axis=1, inplace=True)
        df = df.drop("relatives", axis=1)
        probability_columns = list(df.columns)

        df = id_df.merge(df, right_index=True, left_index=True)

        df.attrs["prob_columns"] = probability_columns

        return df

    def visualize_hierarchy(self, sample_idx: int) -> str:
        """
        Generate ASCII tree visualization of hierarchy for a specific sample.
        
        Args:
            sample_idx: Index of the sample/pair to visualize
            
        Returns:
            String containing ASCII tree representation
        """
        def build_tree_lines(node_name, prefix="", is_last=True):
            """Recursively build tree lines for a node and its children"""
            lines = []
            
            # Get node info
            node_idx = self.node_to_idx[node_name]
            posterior = self.matrix[sample_idx, node_idx, 0]
            conditional = self.matrix[sample_idx, node_idx, 1]
            method_idx = self.matrix[sample_idx, node_idx, 2]
            
            # Get method name (reverse lookup)
            method_name = "unknown"
            for method, idx in self.method_to_idx.items():
                if idx == method_idx:
                    method_name = method
                    break
            
            # Format node line
            connector = "└── " if is_last else "├── "
            prob_info = f"P:{posterior:.3f} C:{conditional:.3f}"
            
            if not np.isnan(method_idx) and method_name != "unknown":
                node_line = f"{prefix}{connector}{node_name} ({prob_info}) [{method_name}]"
            else:
                node_line = f"{prefix}{connector}{node_name} ({prob_info})"
                
            lines.append(node_line)
            
            # Get children
            children = list(nx.descendants_at_distance(self.g, node_name, 1))
            
            if children:
                # Sort children for consistent display
                children.sort()
                
                for i, child in enumerate(children):
                    is_last_child = (i == len(children) - 1)
                    
                    # Determine prefix for children
                    if is_last:
                        child_prefix = prefix + "    "
                    else:
                        child_prefix = prefix + "│   "
                    
                    child_lines = build_tree_lines(child, child_prefix, is_last_child)
                    lines.extend(child_lines)
            
            return lines
        
        # Start building from root
        tree_lines = build_tree_lines("relatives")
        
        # Add header
        header = f"Hierarchy for sample {sample_idx}:"
        result = [header, "=" * len(header)] + tree_lines
        
        return "\n".join(result)