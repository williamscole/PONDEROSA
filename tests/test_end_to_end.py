
from tests.generate_files import GeneratePairs

import pytest
from pathlib import Path
import shutil
import subprocess
import yaml
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

@pytest.fixture(scope="module")
def test_data_dir(tmp_path_factory, request):
    """
    Create a persistent test directory that:
    1. Persists on test failure (for debugging)
    2. Reuses the same directory on re-runs
    3. Cleans up only on successful completion
    """
    # Use pytest's cache to store the directory path
    cache = request.config.cache
    cache_key = "test_end_to_end/data_dir"
    
    # Try to get cached directory path
    cached_dir = cache.get(cache_key, None)
    
    if cached_dir and Path(cached_dir).exists():
        # Reuse existing directory
        data_dir = Path(cached_dir)
        print(f"\n♻️  Reusing cached test directory: {data_dir}")
    else:
        # Create new directory
        data_dir = tmp_path_factory.mktemp("ponderosa_test")
        cache.set(cache_key, str(data_dir))
        print(f"\n📁 Created new test directory: {data_dir}")
        
        # Generate test data
        pairs = GeneratePairs(n_pairs=100, simple_segments=False)
        pairs.generate_some(["po", "fs", "hs", "gp", "av", "co"])
        pairs.write_out(str(data_dir / "test"), n_chrom=22)

        # DEBUG: Print what files were actually created
        print(f"\n📋 Files created in {data_dir}:")
        for f in sorted(data_dir.glob("*")):
            print(f"  - {f.name}")
    
    yield data_dir
    
    # Cleanup only if all tests passed
    if request.session.testsfailed == 0:
        shutil.rmtree(data_dir, ignore_errors=True)
        cache.set(cache_key, None)
        print(f"\n🗑️  Cleaned up test directory: {data_dir}")
    else:
        print(f"\n💾 Test failed - keeping directory for debugging: {data_dir}")


class TestEndToEnd:
    """End-to-end tests for PONDEROSA using generated test data."""
    
    def test_1_direct_api_call(self, test_data_dir):
        """Test 1: Call run_ponderosa directly with config object."""
        from ponderosa import run_ponderosa, PonderosaConfig
        from ponderosa.evaluation import EvaluationResults
        from ponderosa.pedigree import PedigreeHierarchy
        
        config_file = test_data_dir / "test_config.yaml"
        assert config_file.exists(), f"Config file not found: {config_file}"
        
        # Load config and validate
        config = PonderosaConfig.from_yaml(config_file)
        config.validate()
        
        # Run PONDEROSA
        results = run_ponderosa(config)

        print(f"\n=== PEDIGREE REGISTRY ===")
        for node, pairs_list in results.registry.nodes.items():
            if pairs_list:
                print(f"{node}: {len(pairs_list)} pairs")
        
        # Verify outputs exist
        expected_files = [
            f"{config.output.output}_pairs.txt",
            f"{config.output.output}.mhier.pkl",
        ]
        
        for expected_file in expected_files:
            output_file = Path(expected_file)
            assert output_file.exists(), f"Expected output file not found: {output_file}"
        
        # EVALUATE ACCURACY
        truth_file = test_data_dir / "test_truth.txt"
        assert truth_file.exists(), f"Truth file not found: {truth_file}"
        
        # Load the saved matrix hierarchy
        with open(f"{config.output.output}.mhier.pkl", "rb") as f:
            matrix_hierarchy = pickle.load(f)
        
        # Load hierarchy
        hierarchy = PedigreeHierarchy.from_yaml()

        # Evaluate
        evaluation = EvaluationResults.from_ponderosa(
            truth_file=truth_file,
            pred_mhier=matrix_hierarchy,
            hierarchy=hierarchy,
            level_names=["Degree", "Specific"]
        )

        with open(f"{config.output.output}.mhier.pkl", "rb") as f:
            matrix_hierarchy = pickle.load(f)

        # Get predictions
        pred_df = matrix_hierarchy.to_dataframe(min_p=0.0)

        # Load pairs to get IBD features
        from ponderosa.data_loading import load_pairs
        pairs = load_pairs(config.files, config.algorithm)

        # Get IBD1 and IBD2 for all pairs
        ibd_data = []
        for idx, (id1, id2) in matrix_hierarchy.index_to_pair.items():
            ibd1, ibd2 = pairs.get_pair_data_from(
                np.array([(id1, id2)]),  # <-- Convert to numpy array
                "IBD1", "IBD2", 
                output_style="flatten"
            )[0]
            
            # Get truth relationship
            truth_rel = None
            for (tid1, tid2), trel in evaluation.truth_dict.items():
                if (tid1, tid2) == (id1, id2) or (tid2, tid1) == (id1, id2):
                    truth_rel = trel[0] if len(trel) > 0 else None  # Degree level
                    break
            
            ibd_data.append({
                'ibd1': ibd1,
                'ibd2': ibd2,
                'truth': truth_rel,
                'pred': pred_df.iloc[idx]['degree']
            })

        # Convert to dataframe
        plot_df = pd.DataFrame(ibd_data)

        # Create scatter plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: Colored by truth
        for rel in plot_df['truth'].unique():
            if pd.notna(rel):
                subset = plot_df[plot_df['truth'] == rel]
                ax1.scatter(subset['ibd1'], subset['ibd2'], label=rel, alpha=0.6, s=50)
        ax1.set_xlabel('IBD1 (cM)')
        ax1.set_ylabel('IBD2 (cM)')
        ax1.set_title('Ground Truth Relationships')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Colored by prediction
        for rel in plot_df['pred'].unique():
            if pd.notna(rel):
                subset = plot_df[plot_df['pred'] == rel]
                ax2.scatter(subset['ibd1'], subset['ibd2'], label=rel, alpha=0.6, s=50)
        ax2.set_xlabel('IBD1 (cM)')
        ax2.set_ylabel('IBD2 (cM)')
        ax2.set_title('Predicted Relationships')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('ibd_scatter.png', dpi=150)
        print("\n✓ Saved scatter plot to ibd_scatter.png")
        plt.show()

        import pytest; pytest.set_trace()
        
        # Print evaluation
        print("\n" + str(evaluation))
        
        # Assert 100% accuracy at both levels
        metrics = evaluation.get_overall_metrics()
        
        degree_acc = metrics["Degree"]["accuracy"]
        specific_acc = metrics["Specific"]["accuracy"]
        
        assert degree_acc == 1.0, f"Expected 100% degree accuracy, got {degree_acc:.2%}"
        assert specific_acc == 1.0, f"Expected 100% specific accuracy, got {specific_acc:.2%}"
        
        print(f"✓ Test 1 passed: Direct API call successful with 100% accuracy!")
    
    def test_2_cli_with_config(self, test_data_dir):
        """Test 2: Run PONDEROSA from CLI using config file."""
        config_file = test_data_dir / "test_config.yaml"
        assert config_file.exists(), f"Config file not found: {config_file}"
        
        # Run PONDEROSA via subprocess
        result = subprocess.run(
            ["python", "-m", "ponderosa", "--config", str(config_file)],
            capture_output=True,
            text=True
        )
        
        # Check exit code
        assert result.returncode == 0, f"PONDEROSA failed with exit code {result.returncode}\nStderr: {result.stderr}"
        
        # Load config to get output prefix
        with open(config_file) as f:
            config_dict = yaml.safe_load(f)
        output_prefix = config_dict["output"]["output"]
        
        # Verify outputs exist
        expected_files = [
            f"{output_prefix}_pairs.txt",
            f"{output_prefix}.mhier.pkl",
        ]
        
        for expected_file in expected_files:
            output_file = Path(expected_file)
            assert output_file.exists(), f"Expected output file not found: {output_file}"
        
        print(f"✓ Test 2 passed: CLI with config successful")
    
    def test_3_cli_without_config(self, test_data_dir):
        """Test 3: Run PONDEROSA from CLI with individual arguments (no config file)."""
        config_file = test_data_dir / "test_config.yaml"
        assert config_file.exists(), f"Config file not found: {config_file}"
        
        # Parse config to extract individual arguments
        with open(config_file) as f:
            config_dict = yaml.safe_load(f)
        
        files = config_dict["files"]
        algorithm = config_dict.get("algorithm", {})
        output = config_dict.get("output", {})
        
        # Build CLI command
        cmd = [
            "python", "-m", "ponderosa",
            "--ibd", files["ibd"],
            "--fam", files["fam"],
            "--ibd-caller", files["ibd_caller"],
            "--output", output["output"] + "_cli",  # Different output to avoid conflicts
        ]
        
        # Add optional arguments
        if files.get("ages"):
            cmd.extend(["--ages", files["ages"]])
        if files.get("mapf"):
            cmd.extend(["--mapf", files["mapf"]])
        if algorithm.get("min_segment_length"):
            cmd.extend(["--min-segment-length", str(algorithm["min_segment_length"])])
        if algorithm.get("min_total_ibd"):
            cmd.extend(["--min-total-ibd", str(algorithm["min_total_ibd"])])
        if algorithm.get("max_gap"):
            cmd.extend(["--max-gap", str(algorithm["max_gap"])])
        if output.get("verbose"):
            cmd.append("--verbose")
        
        # Run PONDEROSA via subprocess
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        # Check exit code
        assert result.returncode == 0, f"PONDEROSA failed with exit code {result.returncode}\nStderr: {result.stderr}"
        
        # Verify outputs exist
        output_prefix = output["output"] + "_cli"
        expected_files = [
            f"{output_prefix}_pairs.txt",
            f"{output_prefix}.mhier.pkl",
        ]
        
        for expected_file in expected_files:
            output_file = Path(expected_file)
            assert output_file.exists(), f"Expected output file not found: {output_file}"
        
        print(f"✓ Test 3 passed: CLI without config successful")

if __name__ == "__main__":

    pairs = GeneratePairs(n_pairs=100, simple_segments=False)

    pairs.generate_some(["po", "fs", "hs", "gp", "av", "co"])

    pairs.write_out("del_now/test", n_chrom=22)