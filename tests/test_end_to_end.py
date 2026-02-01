
from tests.generate_files import GeneratePairs

import pytest
from pathlib import Path
import shutil
import subprocess
import yaml


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
        
        config_file = test_data_dir / "test_config.yaml"
        assert config_file.exists(), f"Config file not found: {config_file}"
        
        # Load config and run PONDEROSA
        config = PonderosaConfig.from_yaml(config_file)
        results = run_ponderosa(config)
        
        # Verify outputs exist
        output_prefix = config.output.output
        expected_files = [
            f"{output_prefix}_pairs.txt",
            f"{output_prefix}.mhier.pkl",
        ]
        
        for expected_file in expected_files:
            output_file = Path(expected_file)
            assert output_file.exists(), f"Expected output file not found: {output_file}"
        
        print(f"✓ Test 1 passed: Direct API call successful")
        print(results.summary())
    
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