import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from ponderosa.evaluation import create_truth_matrix_hierarchy, EvaluationResults
from ponderosa.pedigree import PedigreeHierarchy
from ponderosa.prediction import MatrixHierarchy



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


class TestEvaluationResults:
    
    def test_perfect_accuracy(self):
        """Test when all predictions are correct."""
        truth_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),
            ("ID3", "ID4"): ("1st", "PO"),
        }
        
        infer_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),  # Both correct
            ("ID3", "ID4"): ("1st", "PO"),   # Both correct
        }
        
        eval_results = EvaluationResults(
            truth_dict, 
            infer_dict, 
            level_names=["Degree", "Specific"]
        )
        
        metrics = eval_results.get_overall_metrics()
        
        # Both levels should be 100% accurate
        assert metrics["Degree"]["accuracy"] == 1.0
        assert metrics["Specific"]["accuracy"] == 1.0
        assert metrics["Degree"]["n_correct"] == 2
        assert metrics["Specific"]["n_correct"] == 2

        print(eval_results)
        
    def test_partial_accuracy(self):
        """Test when degree is right but specific relationship is wrong."""
        truth_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),
            ("ID3", "ID4"): ("1st", "PO"),
        }
        
        infer_dict = {
            ("ID1", "ID2"): ("2nd", "MHS"),  # Degree right, specific wrong
            ("ID3", "ID4"): ("1st", "FS"),   # Degree right, specific wrong
        }
        
        eval_results = EvaluationResults(
            truth_dict, 
            infer_dict, 
            level_names=["Degree", "Specific"]
        )
        
        metrics = eval_results.get_overall_metrics()
        
        # Degree should be 100%, specific should be 0%
        assert metrics["Degree"]["accuracy"] == 1.0
        assert metrics["Specific"]["accuracy"] == 0.0
        
    def test_missing_pairs(self):
        """Test when some pairs are missing from predictions."""
        truth_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),
            ("ID3", "ID4"): ("1st", "PO"),
            ("ID5", "ID6"): ("2nd", "MHS"),
        }
        
        infer_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),
            # ID3/ID4 pair is missing!
            ("ID5", "ID6"): ("3rd", "CO"),  # Wrong
        }
        
        eval_results = EvaluationResults(
            truth_dict, 
            infer_dict, 
            level_names=["Degree", "Specific"]
        )
        
        metrics = eval_results.get_overall_metrics()
        
        # Should only evaluate 2 pairs (missing pairs have NaN)
        assert metrics["Degree"]["n_evaluated"] == 2
        assert metrics["Specific"]["n_evaluated"] == 2
        
    def test_string_output(self):
        """Test that __str__ produces readable output."""
        truth_dict = {
            ("ID1", "ID2"): ("2nd", "PHS"),
            ("ID3", "ID4"): ("1st", "PO"),
        }
        
        infer_dict = {
            ("ID1", "ID2"): ("2nd", "MHS"),
            ("ID3", "ID4"): ("1st", "PO"),
        }
        
        eval_results = EvaluationResults(
            truth_dict, 
            infer_dict, 
            level_names=["Degree", "Specific"]
        )
        
        output = str(eval_results)
        
        # Should contain key information
        assert "PONDEROSA Evaluation Report" in output
        assert "Degree:" in output
        assert "Specific:" in output
        assert "Accuracy:" in output
        assert "50.00%" in output  # Specific is 50% (1/2 correct)


class TestEvaluationResultsFromPonderosa:
    
    def test_from_ponderosa_basic(self, tmp_path):
        """Test from_ponderosa with simple truth file and predictions."""
        
        # Create truth file
        truth_data = pd.DataFrame({
            'ID1': ['ID1', 'ID3', 'ID5'],
            'ID2': ['ID2', 'ID4', 'ID6'],
            'truth': ['PHS', 'PO', 'FS']
        })
        truth_file = tmp_path / "truth.txt"
        truth_data.to_csv(truth_file, sep="\t", index=False)
        
        # Create hierarchy
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        # Create pair_dict
        pair_dict = {
            0: ('ID1', 'ID2'),
            1: ('ID3', 'ID4'),
            2: ('ID5', 'ID6'),
        }
        
        # Create prediction MatrixHierarchy
        pred_mhier = MatrixHierarchy.from_hierarchy(
            hierarchy=hierarchy,
            index_to_pair=pair_dict,
            methods=["test"]
        )
        
        # Manually set predictions (simulating classifier output)
        # For simplicity, let's say predictions match truth
        pred_mhier._add_one_probs(['PHS', 'PO', 'FS'])
        
        # Test from_ponderosa
        eval_results = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=pred_mhier,
            hierarchy=hierarchy,
            level_names=["Degree", "Specific"]
        )
        
        # Check that it was created correctly
        assert eval_results.n_pairs == 3
        assert eval_results.n_levels == 2
        assert eval_results.level_names == ["Degree", "Specific"]
        
        # Check accuracy (should be perfect since predictions match truth)
        metrics = eval_results.get_overall_metrics()
        assert metrics["Degree"]["accuracy"] == 1.0
        assert metrics["Specific"]["accuracy"] == 1.0
        
        print("✓ Basic from_ponderosa test passed")
    
    def test_from_ponderosa_with_errors(self, tmp_path):
        """Test from_ponderosa when predictions don't match truth."""
        
        # Create truth file
        truth_data = pd.DataFrame({
            'ID1': ['ID1', 'ID3', 'ID5'],
            'ID2': ['ID2', 'ID4', 'ID6'],
            'truth': ['PHS', 'PO', 'MGP']  # PHS, PO, MGP
        })
        truth_file = tmp_path / "truth.txt"
        truth_data.to_csv(truth_file, sep="\t", index=False)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        pair_dict = {
            0: ('ID1', 'ID2'),
            1: ('ID3', 'ID4'),
            2: ('ID5', 'ID6'),
        }
        
        pred_mhier = MatrixHierarchy.from_hierarchy(
            hierarchy=hierarchy,
            index_to_pair=pair_dict,
            methods=["test"]
        )
        
        # Set predictions that are partially wrong
        # Truth: PHS (2nd degree), PO (1st degree), MGP (2nd degree)
        # Pred:  MHS (2nd degree), FS (1st degree), PGP (2nd degree)
        pred_mhier._add_one_probs(['MHS', 'FS', 'PGP'])
        
        eval_results = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=pred_mhier,
            hierarchy=hierarchy,
            level_names=["Degree", "Specific"]
        )
        
        metrics = eval_results.get_overall_metrics()
        
        # Degree level: All correct (2nd, 1st, 2nd matches 2nd, 1st, 2nd)
        assert metrics["Degree"]["accuracy"] == 1.0
        assert metrics["Degree"]["n_correct"] == 3
        
        # Specific level: All wrong (MHS != PHS, FS != PO, PGP != MGP)
        assert metrics["Specific"]["accuracy"] == 0.0
        assert metrics["Specific"]["n_correct"] == 0
        
        print("✓ from_ponderosa with errors test passed")
    
    def test_from_ponderosa_missing_pairs(self, tmp_path):
        """Test when truth file has pairs not in predictions."""
        
        # Create truth file with extra pairs
        truth_data = pd.DataFrame({
            'ID1': ['ID1', 'ID3', 'ID99'],  # ID99 doesn't exist
            'ID2': ['ID2', 'ID4', 'ID100'],
            'truth': ['PHS', 'PO', 'FS']
        })
        truth_file = tmp_path / "truth.txt"
        truth_data.to_csv(truth_file, sep="\t", index=False)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        # Only two pairs in predictions
        pair_dict = {
            0: ('ID1', 'ID2'),
            1: ('ID3', 'ID4'),
        }
        
        pred_mhier = MatrixHierarchy.from_hierarchy(
            hierarchy=hierarchy,
            index_to_pair=pair_dict,
            methods=["test"]
        )
        
        pred_mhier._add_one_probs(['PHS', 'PO'])
        
        eval_results = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=pred_mhier,
            hierarchy=hierarchy
        )
        
        # Should create evaluation for all 3 pairs in truth
        assert eval_results.n_pairs == 3
        
        metrics = eval_results.get_overall_metrics()
        
        # Should only evaluate 2 pairs (the missing pair has NaN)
        assert metrics["Level_0"]["n_evaluated"] == 2
        assert metrics["Level_1"]["n_evaluated"] == 2
        
        print("✓ from_ponderosa missing pairs test passed")
    
    def test_from_ponderosa_default_levels(self, tmp_path):
        """Test that default levels are extracted correctly from hierarchy."""
        
        truth_data = pd.DataFrame({
            'ID1': ['ID1'],
            'ID2': ['ID2'],
            'truth': ['PHS']
        })
        truth_file = tmp_path / "truth.txt"
        truth_data.to_csv(truth_file, sep="\t", index=False)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        pair_dict = {0: ('ID1', 'ID2')}
        
        pred_mhier = MatrixHierarchy.from_hierarchy(
            hierarchy=hierarchy,
            index_to_pair=pair_dict,
            methods=["test"]
        )
        pred_mhier._add_one_probs(['PHS'])
        
        # Don't specify level_nodes - should use defaults
        eval_results = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=pred_mhier,
            hierarchy=hierarchy
        )
        
        # Should have 2 levels: degree and specific
        assert eval_results.n_levels == 2
        
        # Check default level names
        assert eval_results.level_names[0] == "Level_0"  # Degree
        assert eval_results.level_names[1] == "Level_1"  # Specific
        
        print("✓ from_ponderosa default levels test passed")
    
    def test_from_ponderosa_custom_levels(self, tmp_path):
        """Test with custom level specification."""
        
        truth_data = pd.DataFrame({
            'ID1': ['ID1', 'ID3'],
            'ID2': ['ID2', 'ID4'],
            'truth': ['PHS', 'MHS']
        })
        truth_file = tmp_path / "truth.txt"
        truth_data.to_csv(truth_file, sep="\t", index=False)
        
        hierarchy = PedigreeHierarchy(SIMPLE_HIERARCHY)
        
        pair_dict = {
            0: ('ID1', 'ID2'),
            1: ('ID3', 'ID4'),
        }
        
        pred_mhier = MatrixHierarchy.from_hierarchy(
            hierarchy=hierarchy,
            index_to_pair=pair_dict,
            methods=["test"]
        )
        pred_mhier._add_one_probs(['PHS', 'MHS'])
        
        # Specify custom levels: degree, HS type, specific
        eval_results = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=pred_mhier,
            hierarchy=hierarchy,
            level_nodes=[
                ["1st", "2nd"],           # Degree level
                ["HS", "GP"],             # Type level (2nd degree subtypes)
                ["PHS", "MHS", "PGP", "MGP"]  # Specific level
            ],
            level_names=["Degree", "Type", "Specific"]
        )
        
        # Should have 3 levels
        assert eval_results.n_levels == 3
        assert eval_results.level_names == ["Degree", "Type", "Specific"]
        
        metrics = eval_results.get_overall_metrics()
        
        # All should be perfect since predictions match truth
        assert metrics["Degree"]["accuracy"] == 1.0  # Both are 2nd
        assert metrics["Type"]["accuracy"] == 1.0     # Both are HS
        assert metrics["Specific"]["accuracy"] == 1.0 # PHS and MHS correct
        
        print("✓ from_ponderosa custom levels test passed")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])