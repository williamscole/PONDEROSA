import numpy as np
import yaml
from pathlib import Path
import pandas as pd

class AgePriors:

    def __init__(self, prior_dict: dict):

        self.hard_cutoffs = {}

        for node, prior_info in prior_dict.items():
            # Hard cut off for age is given for the relationship
            if "min_gap" in prior_info or "max_gap" in prior_info:
                # If no min_gap specified, min_gap be as small/negative as possible
                min_gap = prior_info.get("min_gap", -np.inf)
                # If no max_gap specified, max_gap can be as large as possible
                max_gap = prior_info.get("max_gap", np.inf)

                # Function returns True if the age gap is too small or too large
                self.hard_cutoffs[node] = lambda age_gap: age_gap < min_gap or age_gap > max_gap


    def _hard_cut_zeros(self, age_gap: float):

        if pd.isna(age_gap) or age_gap is None:
            return []
        
        zero_prob_nodes = []

        for node, node_func in self.hard_cutoffs.items():
            if node_func(age_gap):
                zero_prob_nodes.append(node)

        return zero_prob_nodes
    
    def get_hard_cut_zeros(self, age_gap_arr: np.ndarray):
        n_pairs = age_gap_arr.shape[0]
        out_dict = {node: np.zeros(n_pairs)+1 for node in self.hard_cutoffs.keys()}

        for pair_idx, age_gap in enumerate(age_gap_arr):
            for node in self._hard_cut_zeros(age_gap):
                # Set the probability to zero for this pair
                out_dict[node][pair_idx] = 0

        return out_dict
    
    @classmethod
    def from_yaml(cls, yaml_file: Path):

        assert yaml_file.exists()

        with open(yaml_file, 'r') as file:
            prior_data = yaml.safe_load(file)
        
        return cls(prior_data["age_priors"])



    
        

    



