from pathlib import Path
from typing import Dict, Tuple, List
import numpy as np
import pandas as pd

from .prediction import MatrixHierarchy
from .pedigree import PedigreeHierarchy
from .data_loading import load_truth
        
def create_truth_matrix_hierarchy(truth_df: pd.DataFrame, 
                                   hierarchy: PedigreeHierarchy,
                                   pair_dict: Dict[int, Tuple[str, str]]) -> MatrixHierarchy:
    """
    Create a MatrixHierarchy from ground truth where each relationship has P=1.0
    """
    
    # Create MatrixHierarchy
    truth_mhier = MatrixHierarchy.from_hierarchy(
        hierarchy=hierarchy,
        index_to_pair=pair_dict,
        methods=["truth"]
    )

    truth_df = truth_df[truth_df[["ID1", "ID2"]].apply(tuple, axis=1).isin(pair_dict.values())]

    truth_mhier._add_one_probs(truth_df.truth.values)
    
    return truth_mhier

class EvaluationResults:

    def __init__(self,
                 truth_dict: Dict[Tuple[str, str], Tuple[str]],  # maps ID pair to relationship tuple
                 infer_dict: Dict[Tuple[str, str], Tuple[str]],
                 level_names: List[str] = None
                 ):
                
        n_levels = len(list(truth_dict.values())[0])

        if level_names is None:
            level_names = [f"Level {i}" for i in range(1, n_levels+1)]
        else:
            assert len(level_names) == n_levels

        results = np.full((len(truth_dict), n_levels), np.nan)

        for pair_idx, (pair, truth) in enumerate(truth_dict.items()):

            inferred = infer_dict.get(pair, None)

            if inferred:
                # Must have the same level of relationship!
                assert len(truth) == n_levels and \
                       len(inferred) == n_levels
                
                for idx in range(n_levels):
                    results[pair_idx, idx] = int(truth[idx] == inferred[idx])

        self.results = results
        self.n_levels = n_levels
        self.level_names = level_names
        self.n_pairs = len(truth_dict)

    def _calculate_level_accuracy(self, level_idx: int) -> Dict[str, float]:
        """Calculate accuracy metrics for a specific level, ignoring NaN values."""
        level_results = self.results[:, level_idx]
        
        # Filter out NaN values
        valid_mask = ~np.isnan(level_results)
        valid_results = level_results[valid_mask]
        
        if len(valid_results) == 0:
            return {
                'n_evaluated': 0,
                'n_correct': 0,
                'n_incorrect': 0,
                'accuracy': np.nan
            }
        
        n_correct = np.sum(valid_results == 1)
        n_incorrect = np.sum(valid_results == 0)
        accuracy = n_correct / len(valid_results)
        
        return {
            'n_evaluated': len(valid_results),
            'n_correct': int(n_correct),
            'n_incorrect': int(n_incorrect),
            'accuracy': accuracy
        }
    
    def get_overall_metrics(self) -> Dict[str, Dict]:
        """Get metrics for all levels."""
        metrics = {}
        for level_idx in range(self.n_levels):
            level_name = self.level_names[level_idx]
            metrics[level_name] = self._calculate_level_accuracy(level_idx)
        return metrics
    
    def __str__(self) -> str:
        """Return formatted evaluation report as string."""
        lines = []
        lines.append("=" * 60)
        lines.append("PONDEROSA Evaluation Report")
        lines.append("=" * 60)
        lines.append(f"\nTotal Pairs in Truth: {self.n_pairs}")
        
        metrics = self.get_overall_metrics()
        
        for level_name, level_metrics in metrics.items():
            lines.append("\n" + "-" * 60)
            lines.append(f"{level_name}:")
            lines.append("-" * 60)
            
            if level_metrics['n_evaluated'] == 0:
                lines.append("  No pairs evaluated at this level")
            else:
                lines.append(f"  Pairs Evaluated: {level_metrics['n_evaluated']}")
                lines.append(f"  Correct:         {level_metrics['n_correct']}")
                lines.append(f"  Incorrect:       {level_metrics['n_incorrect']}")
                lines.append(f"  Accuracy:        {level_metrics['accuracy']:.2%}")
        
        lines.append("\n" + "=" * 60)
        
        return "\n".join(lines)
    
    @classmethod
    def from_ponderosa(cls,
                    truth_file: Path,
                    pred_mhier: MatrixHierarchy,
                    hierarchy: PedigreeHierarchy,
                    level_nodes: List[List[str]] = None,
                    level_names: List[str] = None):
        """
        Create EvaluationResults from PONDEROSA outputs.
        
        Compares ground truth relationships to predictions at multiple hierarchical levels.
        For each level, uses most_likely_among to get the predicted node at that level.
        
        Args:
            truth_file: Path to truth file with columns [ID1, ID2, truth]
            pred_mhier: MatrixHierarchy with predicted relationships
            hierarchy: PedigreeHierarchy defining relationship structure
            level_nodes: Nodes to evaluate at each level.
                        e.g., [["1st", "2nd", "3rd", "4th"], ["PHS", "MHS", "PGP", "MGP", ...]]
                        If None, uses default levels (degree, then specific)
            level_names: Names for each level (e.g., ["Degree", "Specific Relationship"])
                        If None, generates names like "Level_0", "Level_1"
        
        Returns:
            EvaluationResults object with accuracy metrics at each hierarchy level
        
        Example:
            >>> eval_results = EvaluationResults.from_ponderosa(
            ...     truth_file=Path("truth.txt"),
            ...     pred_mhier=matrix_hierarchy,
            ...     hierarchy=hierarchy,
            ...     level_names=["Degree", "Specific"]
            ... )
            >>> print(eval_results)
        """
        from .data_loading import load_truth
        
        # Load truth data
        truth_df = load_truth(truth_file)
        
        # Determine levels to evaluate if not provided
        if level_nodes is None:
            level_nodes = cls._get_default_levels(hierarchy)
        
        n_levels = len(level_nodes)
        
        # Generate level names if not provided
        if level_names is None:
            level_names = [f"Level_{i}" for i in range(n_levels)]
        else:
            assert len(level_names) == n_levels, \
                f"level_names length ({len(level_names)}) must match level_nodes length ({n_levels})"
        
        # Build truth_dict: (id1, id2) -> tuple of relationships at each level
        truth_dict = {}
        for _, row in truth_df.iterrows():
            pair = (row['ID1'], row['ID2'])
            relationship = row['truth']
            
            # Get the hierarchical path for this relationship
            try:
                rel_path = cls._get_relationship_hierarchy(relationship, hierarchy)
            except ValueError as e:
                # Skip relationships not in hierarchy
                continue
            
            # For each level, find which node matches
            rel_at_levels = []
            for level in level_nodes:
                matching_node = None
                for node in level:
                    if node in rel_path:
                        matching_node = node
                        break
                rel_at_levels.append(matching_node)
            
            truth_dict[pair] = tuple(rel_at_levels)
        
        # Build infer_dict: (id1, id2) -> tuple of predictions at each level
        infer_dict = {}
        for idx, pair in pred_mhier.index_to_pair.items():
            pred_at_levels = []
            
            for level in level_nodes:
                # Use most_likely_among to get prediction at this specific level
                pred_array = pred_mhier.most_likely_among(level, asint=False)
                pred_at_levels.append(pred_array[idx])
            
            infer_dict[pair] = tuple(pred_at_levels)
        
        return cls(truth_dict, infer_dict, level_names)

    @staticmethod
    def _get_relationship_hierarchy(relationship: str, hierarchy: PedigreeHierarchy) -> List[str]:
        """
        Get the hierarchical path from root to relationship.
        
        Args:
            relationship: Specific relationship (e.g., "PHS")
            hierarchy: PedigreeHierarchy
        
        Returns:
            List of nodes from most general to most specific
            e.g., ["2nd", "HS", "PHS"] for PHS
            e.g., ["1st", "PO"] for PO
        
        Raises:
            ValueError: If relationship not found in hierarchy or no path exists
        """
        import networkx as nx
        
        if relationship not in hierarchy.nodes:
            raise ValueError(f"Relationship '{relationship}' not found in hierarchy")
        
        # Get path from 'relatives' to this relationship
        try:
            path = nx.shortest_path(hierarchy, 'relatives', relationship)
        except nx.NetworkXNoPath:
            raise ValueError(f"No path from 'relatives' to '{relationship}'")
        
        # Remove 'relatives' from path (it's the root, not a meaningful level)
        path = path[1:]
        
        return path

    @staticmethod
    def _get_default_levels(hierarchy: PedigreeHierarchy) -> List[List[str]]:
        """
        Get default evaluation levels from hierarchy.
        
        Returns two levels:
        1. Degree level: Direct children of 'relatives' (e.g., ["1st", "2nd", "3rd", "4th"])
        2. Specific level: All leaf nodes (relationships with no children)
        
        Args:
            hierarchy: PedigreeHierarchy
        
        Returns:
            List of lists, where each inner list contains nodes at that level
            e.g., [["1st", "2nd", "3rd", "4th"], ["PO", "FS", "PHS", "MHS", ...]]
        
        Example:
            >>> levels = EvaluationResults._get_default_levels(hierarchy)
            >>> degree_level, specific_level = levels
            >>> print(degree_level)  # ["1st", "2nd", "3rd", "4th"]
        """
        import networkx as nx
        
        # Level 1: Degree nodes (direct children of 'relatives')
        degree_nodes = sorted(list(hierarchy.successors('relatives')))
        
        # Level 2: Leaf nodes (specific relationships with no children)
        leaf_nodes = sorted([
            node for node in hierarchy.nodes 
            if hierarchy.out_degree(node) == 0 and node != 'relatives'
        ])
        
        return [degree_nodes, leaf_nodes]