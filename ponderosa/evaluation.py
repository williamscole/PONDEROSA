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
    # Create reverse mapping: (id1, id2) -> index
    pair_to_index = {pair: idx for idx, pair in pair_dict.items()}
    
    # Create MatrixHierarchy
    truth_mhier = MatrixHierarchy.from_hierarchy(
        hierarchy=hierarchy,
        index_to_pair=pair_dict,
        methods=["truth"]
    )

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
                       truth_file: Path):
        """Logic here - to be implemented"""
        raise NotImplementedError("from_ponderosa not yet implemented")