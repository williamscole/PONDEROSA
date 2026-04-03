"""
Script for managing/organizing the simulation
"""

from pathlib import Path
import subprocess
import gzip

from ..config import PonderosaConfig
from ..core import run_ponderosa
from .config import SimulationConfig
from .setup import simulation_workspace
from .founders import pedsim_dryrun, calculate_relatedness, create_founders_file
from .pedsim import PedSim


def simulate(config: SimulationConfig) -> PonderosaConfig:
    """
    Main simulation workflow.
    
    Returns a PonderosaConfig for running PONDEROSA analysis.
    """
    
    with simulation_workspace(config) as temp_dir:
        
        vcf_ext = ".vcf.gz" if str(config.pedsim.vcf_file).endswith('.gz') else ".vcf"

        pedsim = PedSim(
            vcf_file=str(temp_dir / f"input{vcf_ext}"),  # Use dynamic extension
            def_file=str(temp_dir / "pedigree.def"),
            intf_file=str(temp_dir / "interference.tsv"),
            map_file=str(temp_dir / "input.map"),
            output=str(temp_dir / "simulation"),
            executable_path=str(config.pedsim.pedsim_executable),
        )
        
        # 2. Dry run to get family structures
        dry_run_families = pedsim_dryrun(pedsim)
        
        # 3. Calculate relatedness from KING file
        relatedness_dict = calculate_relatedness(
            ibd_file=[str(config.king_file)],
            ibd_caller="king",  # or from config
            genetic_map_file_list=[]  # Add map files if needed
        )
        
        # 4. Get VCF samples
        vcf_samples = _get_vcf_samples(vcf_file=Path(pedsim.get_input("vcf")))
        
        # 5. Create founders mapping
        founders_df = create_founders_file(
            vcf_samples=vcf_samples,
            dry_run_families=dry_run_families,
            relatedness_dict=relatedness_dict,
            max_k=config.training.max_kinship,
            n_sim=config.training.n_pairs_per_relationship
        )
        
        # 6. Write founders file
        founders_file = temp_dir / "founders.txt"
        founders_df.to_csv(founders_file, sep="\t", index=False, header=False)
        
        # 7. Update PedSim with founders and execute
        pedsim.update_flag("--set_founders", str(founders_file))
        pedsim.execute()
        
        # 8. Run IBD caller on simulated output
        _run_ibd_caller(config.ibd_caller, temp_dir)
        
        # 9. Create PonderosaConfig from simulation results
        ponderosa_config = _create_ponderosa_config(temp_dir, config)
        ponderosa_config.algorithm.train_only = True
        ponderosa_config.validate()

        # 10. Run PONDEROSA
        run_ponderosa(ponderosa_config, train_only=True)
        
        return Path(f"{config.output_path}.classif.pkl")


def _get_vcf_samples(vcf_file: Path) -> list[str]:
    """Extract sample IDs from VCF header."""
    import gzip
    
    # Determine how to open the file
    if str(vcf_file).endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt', encoding='utf-8')
    else:
        open_func = lambda f: open(f, 'r', encoding='utf-8')
    
    with open_func(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Header line found - extract sample IDs (columns 9 onward)
                columns = line.strip().split('\t')
                return columns[9:]
    
    # If we get here, no header line was found
    raise ValueError(f"No #CHROM header line found in VCF file: {vcf_file}")


def _run_ibd_caller(ibd_caller: str, temp_dir: Path):
    """
    Run the user's IBD calling script on simulated output.
    
    The script receives the temp directory as its only argument.
    It is expected to find the simulated VCF at {temp_dir}/simulation.vcf(.gz)
    and write IBD output with the prefix {temp_dir}/simulation.
    
    Args:
        ibd_caller: Path to the user's bash script
        temp_dir: Simulation workspace directory
    """
    script_path = Path(ibd_caller).resolve()
    
    if not script_path.exists():
        raise FileNotFoundError(f"IBD caller script not found: {script_path}")
    
    result = subprocess.run(
        ["bash", str(script_path), str(temp_dir)],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(
            f"IBD caller script failed (exit code {result.returncode}):\n"
            f"stderr: {result.stderr}\n"
            f"stdout: {result.stdout}"
        )


def _create_ponderosa_config(temp_dir: Path, sim_config: SimulationConfig) -> PonderosaConfig:
    """
    Create PonderosaConfig from simulation output.
    
    Expects the following files in temp_dir after steps 7-8:
        - simulation-everyone.fam   (from ped-sim)
        - simulation.ibd.gz         (from IBD caller)
        - genetic.map               (from IBD caller script)
    """
    
    # Locate the fam file from ped-sim output
    fam_file = temp_dir / "simulation-everyone.fam"
    if not fam_file.exists():
        raise FileNotFoundError(f"Ped-sim fam file not found: {fam_file}")
    
    # Locate IBD output
    ibd_file = temp_dir / "simulation.ibd.gz"
    if not ibd_file.exists():
        raise FileNotFoundError(
            f"IBD output not found: {ibd_file}\n"
            f"Check that the IBD caller script wrote output with prefix: {temp_dir / 'simulation'}"
        )
    
    # Derive the IBD format from the script name (e.g. "hap-ibd.sh" -> "hap-ibd")
    ibd_caller = Path(sim_config.ibd_caller).stem
    
    config_dict = {
        "files": {
            "ibd": str(ibd_file),
            "fam": str(fam_file),
            "ibd_caller": ibd_caller,
        },
        "algorithm": {
            "min_segment_length": 3.0,
            "min_total_ibd": 50.0,
        },
        "output": {
            "output": sim_config.output_path,
            "write_training": True,
        }
    }
    
    # Add map file if the IBD caller script produced one
    map_file = temp_dir / "genetic.map"
    if map_file.exists():
        config_dict["files"]["mapf"] = str(map_file)
    
    return PonderosaConfig.from_dict(config_dict)