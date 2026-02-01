import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from ponderosa.evaluation import create_truth_matrix_hierarchy, EvaluationResults
from ponderosa.pedigree import PedigreeHierarchy


# Simple hierarchy for testing (1st and 2nd degree only)
SIMPLE_HIERARCHY = {
    'PO': {'sex': False, 'degree': '1st', 'path': True, 1: [[1]]},
    'FS': {'sex': False, 'degree': '1st', 'path': True, 1: [[1, -1]], 2: [[1, -1]]},
    'HS': {'sex': False, 'degree': '2nd', 'path': True, 1: [[1, -1]]},
    'PHS': {'sex': True, 'degree': '2nd', 'parent': 'HS', 'path': True, 1: [[1, -1]]},
    'MHS': {'sex': True, 'degree': '2nd', 'parent': 'HS', 'path': True, 2: [[1, -1]]},
    'GP': {'sex': False, 'degree': '2nd', 'path': True, 1: [[1, 1]]},
    'PGP': {'sex': True, 'degree': '2nd', 'parent': 'GP', 'path': True, 1: [[1, 1]]},
    'MGP': {'sex': True, 'degree': '2nd', 'parent': 'GP', 'path': True, 2: [[1, 1]]},
}


class TestCreateTruthMatrixHierarchy:
    
    def test_basic_truth_hierarchy(self):
        """Test creating a truth MatrixHierarchy from simple ground truth."""
        
        # Create simple truth data
        truth_data = {
            'ID1': ['ID1', 'ID3', 'ID5'],
            'ID2': ['ID2', 'ID4', 'ID6'],
            'truth': ['PHS', 'PO', 'FS']
        }
        truth_df = pd.DataFrame(truth_data)
        
        # Create hierarchy from dict
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        # Create pair dict (index to pair mapping)
        all_pairs = [
            ('ID1', 'ID2'),
            ('ID3', 'ID4'),
            ('ID5', 'ID6'),
        ]
        pair_dict = {i: pair for i, pair in enumerate(all_pairs)}
        
        # Create truth matrix hierarchy
        truth_mhier = create_truth_matrix_hierarchy(
            truth_df=truth_df,
            hierarchy=hierarchy,
            pair_dict=pair_dict
        )
        
        # Verify it was created
        assert truth_mhier is not None
        
        # Check that probabilities were set correctly
        pred, prob = truth_mhier.most_probable(1)
        
        # Should predict the correct relationships
        assert pred[0] == 'PHS'  # First pair
        assert pred[1] == 'PO'   # Second pair
        assert pred[2] == 'FS'   # Third pair
        
        # Probabilities should be 1.0 for ground truth
        assert prob[0] == 1.0
        assert prob[1] == 1.0
        assert prob[2] == 1.0
        
        print("✓ Basic truth hierarchy test passed")
    
    def test_truth_hierarchy_propagation(self):
        """Test that probabilities propagate up the hierarchy."""
        
        truth_data = {
            'ID1': ['ID1'],
            'ID2': ['ID2'],
            'truth': ['PHS']  # Paternal half-sibling
        }
        truth_df = pd.DataFrame(truth_data)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        pair_dict = {0: ('ID1', 'ID2')}
        
        truth_mhier = create_truth_matrix_hierarchy(
            truth_df=truth_df,
            hierarchy=hierarchy,
            pair_dict=pair_dict
        )
        
        # Get predictions at different levels
        pred, prob = truth_mhier.most_probable(0.5)
        
        # Should correctly identify as PHS
        assert pred[0] == 'PHS'
        assert prob[0] == 1.0
        
        # PHS should propagate to HS (parent node)
        # (Exact method depends on your MatrixHierarchy interface)
        # You might need to add methods to check probabilities at specific nodes
        
        print("✓ Hierarchy propagation test passed")
    
    def test_missing_pairs(self):
        """Test handling pairs in truth but not in pair_dict."""
        
        truth_data = {
            'ID1': ['ID1', 'ID99'],  # ID99 pair doesn't exist in pair_dict
            'ID2': ['ID2', 'ID100'],
            'truth': ['PHS', 'FS']
        }
        truth_df = pd.DataFrame(truth_data)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        pair_dict = {0: ('ID1', 'ID2')}  # Only one pair
        
        # Should handle missing pairs gracefully (skip or warn)
        truth_mhier = create_truth_matrix_hierarchy(
            truth_df=truth_df,
            hierarchy=hierarchy,
            pair_dict=pair_dict
        )
        
        assert truth_mhier is not None
        
        # Should still work for the pair that exists
        pred, prob = truth_mhier.most_probable(0.5)
        assert len(pred) == 1  # Only one pair in pair_dict
        assert pred[0] == 'PHS'
        
        print("✓ Missing pairs test passed")
    
    def test_multiple_relationships_same_level(self):
        """Test with multiple different relationships at the same hierarchical level."""
        
        truth_data = {
            'ID1': ['ID1', 'ID3', 'ID5', 'ID7'],
            'ID2': ['ID2', 'ID4', 'ID6', 'ID8'],
            'truth': ['PHS', 'MHS', 'PGP', 'MGP']  # All 2nd degree
        }
        truth_df = pd.DataFrame(truth_data)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        pair_dict = {i: (f'ID{2*i+1}', f'ID{2*i+2}') for i in range(4)}
        
        truth_mhier = create_truth_matrix_hierarchy(
            truth_df=truth_df,
            hierarchy=hierarchy,
            pair_dict=pair_dict
        )
        
        pred, prob = truth_mhier.most_probable(0.5)
        
        # All should be predicted correctly
        assert pred[0] == 'PHS'
        assert pred[1] == 'MHS'
        assert pred[2] == 'PGP'
        assert pred[3] == 'MGP'
        
        # All should have probability 1.0
        assert np.all(prob == 1.0)
        
        print("✓ Multiple relationships test passed")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])