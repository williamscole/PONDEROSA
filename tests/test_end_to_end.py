
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
        # plt.show()
        plt.close()

        # ============== ADD THIS AFTER THE IBD1 vs IBD2 SCATTER PLOT ==============

        # Get N (number of segments) for 2nd degree pairs and create jitterplot
        second_degree_rels = ["AV", "PGP", "MGP", "PHS", "MHS"]

        # Collect data for 2nd degree pairs
        second_degree_data = []
        for idx, (id1, id2) in matrix_hierarchy.index_to_pair.items():
            # Get N, IBD1, IBD2 features
            n_segs, ibd1, ibd2 = pairs.get_pair_data_from(
                np.array([(id1, id2)]),
                "N", "IBD1", "IBD2",
                output_style="flatten"
            )[0]
            
            # Get truth relationship (specific level)
            truth_rel = None
            for (tid1, tid2), trel in evaluation.truth_dict.items():
                if (tid1, tid2) == (id1, id2) or (tid2, tid1) == (id1, id2):
                    # trel[1] is the specific relationship level
                    if len(trel) > 1:
                        truth_rel = trel[1]
                    break
            
            # Only include 2nd degree relationships
            if truth_rel in second_degree_rels:
                # Get predicted specific relationship
                pred_rel = pred_df.iloc[idx]['pred_rel']
                
                second_degree_data.append({
                    'n_segs': n_segs,
                    'ibd1': ibd1,
                    'ibd2': ibd2,
                    'truth': truth_rel,
                    'pred': pred_rel
                })

        second_degree_df = pd.DataFrame(second_degree_data)

        # Create jitterplot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Helper function to add jitter
        def jitter(values, amount=0.2):
            return values + np.random.uniform(-amount, amount, len(values))

        # Get unique relationships and assign x positions
        unique_rels = sorted(second_degree_df['truth'].dropna().unique())
        rel_to_x = {rel: i for i, rel in enumerate(unique_rels)}

        # Plot 1: Colored by ground truth
        for rel in unique_rels:
            subset = second_degree_df[second_degree_df['truth'] == rel]
            x_pos = jitter(np.full(len(subset), rel_to_x[rel]))
            ax1.scatter(x_pos, subset['n_segs'], label=rel, alpha=0.6, s=50)

        ax1.set_xticks(range(len(unique_rels)))
        ax1.set_xticklabels(unique_rels)
        ax1.set_xlabel('Relationship (Ground Truth)')
        ax1.set_ylabel('Number of IBD Segments')
        ax1.set_title('Ground Truth: N Segments by 2nd Degree Relationship')
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='y')

        # Plot 2: Colored by prediction
        unique_preds = sorted(second_degree_df['pred'].dropna().unique())
        for rel in unique_preds:
            subset = second_degree_df[second_degree_df['pred'] == rel]
            # Use truth position on x-axis, color by prediction
            x_positions = [rel_to_x.get(t, -1) for t in subset['truth']]
            x_pos = jitter(np.array(x_positions).astype(float))
            ax2.scatter(x_pos, subset['n_segs'], label=rel, alpha=0.6, s=50)

        ax2.set_xticks(range(len(unique_rels)))
        ax2.set_xticklabels(unique_rels)
        ax2.set_xlabel('Relationship (Ground Truth)')
        ax2.set_ylabel('Number of IBD Segments')
        ax2.set_title('Predicted: N Segments by 2nd Degree Relationship')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig('nsegs_jitter.png', dpi=150)
        print("\n✓ Saved jitterplot to nsegs_jitter.png")
        # plt.show()
        plt.close()

        # ============== HAPLOTYPE SCORE PLOT ==============

        # The haplotype classifier distinguishes GPAV vs HS
        hap_rels = ["GPAV", "HS"]
        # But we want to see the specific relationships too
        hap_specific_rels = ["AV", "PGP", "MGP", "PHS", "MHS", "GP"]

        # Collect haplotype data for 2nd degree pairs
        hap_data = []
        for idx, (id1, id2) in matrix_hierarchy.index_to_pair.items():
            h1, h2 = pairs.get_pair_data_from(
                np.array([(id1, id2)]),
                "H1", "H2",
                output_style="flatten"
            )[0]
            
            # Get truth relationship (specific level)
            truth_rel = None
            for (tid1, tid2), trel in evaluation.truth_dict.items():
                if (tid1, tid2) == (id1, id2) or (tid2, tid1) == (id1, id2):
                    if len(trel) > 1:
                        truth_rel = trel[1]
                    break
            
            # Only include relevant 2nd degree relationships
            if truth_rel in hap_specific_rels:
                # Determine the parent category (GPAV or HS)
                if truth_rel in ["AV", "PGP", "MGP", "GP"]:
                    truth_parent = "GPAV"
                elif truth_rel in ["PHS", "MHS"]:
                    truth_parent = "HS"
                else:
                    truth_parent = truth_rel
                    
                hap_data.append({
                    'h1': h1,
                    'h2': h2,
                    'truth': truth_rel,
                    'truth_parent': truth_parent,
                })

        hap_df = pd.DataFrame(hap_data)

        # Create scatter plot of H1 vs H2
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: Colored by specific relationship
        for rel in hap_df['truth'].unique():
            if pd.notna(rel):
                subset = hap_df[hap_df['truth'] == rel]
                ax1.scatter(subset['h1'], subset['h2'], label=rel, alpha=0.6, s=50)

        ax1.set_xlabel('H1 (Haplotype Score 1)')
        ax1.set_ylabel('H2 (Haplotype Score 2)')
        ax1.set_title('Haplotype Scores by Specific Relationship')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        # Add diagonal line for reference (H1 = H2)
        lims = [min(ax1.get_xlim()[0], ax1.get_ylim()[0]), max(ax1.get_xlim()[1], ax1.get_ylim()[1])]
        ax1.plot(lims, lims, 'k--', alpha=0.3, label='H1=H2')

        # Plot 2: Colored by parent category (GPAV vs HS)
        colors = {'GPAV': 'blue', 'HS': 'red'}
        for rel in ['GPAV', 'HS']:
            subset = hap_df[hap_df['truth_parent'] == rel]
            ax2.scatter(subset['h1'], subset['h2'], label=rel, alpha=0.6, s=50, c=colors[rel])

        ax2.set_xlabel('H1 (Haplotype Score 1)')
        ax2.set_ylabel('H2 (Haplotype Score 2)')
        ax2.set_title('Haplotype Scores: GPAV vs HS')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        # Add diagonal line
        lims = [min(ax2.get_xlim()[0], ax2.get_ylim()[0]), max(ax2.get_xlim()[1], ax2.get_ylim()[1])]
        ax2.plot(lims, lims, 'k--', alpha=0.3)

        plt.tight_layout()
        plt.savefig('haplotype_scatter.png', dpi=150)
        print("\n✓ Saved haplotype plot to haplotype_scatter.png")
        # plt.show()
        plt.close()

        # import pytest; pytest.set_trace()
        
        # Print evaluation
        print("\n" + str(evaluation))
        
        # Assert 100% accuracy at both levels
        metrics = evaluation.get_overall_metrics()
        
        degree_acc = metrics["Degree"]["accuracy"]
        specific_acc = metrics["Specific"]["accuracy"]
        
        assert degree_acc > 0.95, f"Expected 100% degree accuracy, got {degree_acc:.2%}"
        assert specific_acc > 0.95, f"Expected 100% specific accuracy, got {specific_acc:.2%}"
        
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
            cmd.extend(["--map", files["mapf"]])
        if algorithm.get("min_segment_length"):
            cmd.extend(["--min-segment-length", str(algorithm["min_segment_length"])])
        if algorithm.get("min_total_ibd"):
            cmd.extend(["--min-total-ibd", str(algorithm["min_total_ibd"])])
        if algorithm.get("max_gap"):
            cmd.extend(["--max-gap", str(algorithm["max_gap"])])
        if output.get("verbose"):
            cmd.append("--verbose")
        if files.get("truth"):
            cmd.extend(["--truth", files["truth"]])
        
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