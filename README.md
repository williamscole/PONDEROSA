## Simulation Pipeline

*(In testing — not yet available for general use)*

PONDEROSA includes a simulation pipeline for generating population-specific training data using [Ped-sim](https://github.com/williamslab/ped-sim). The simulation workflow simulates relative pairs from real genotype data, detects IBD segments, and trains classifiers on the resulting IBD sharing patterns. The output is a trained classifier pickle file (`.classif.pkl`) that can be passed directly to PONDEROSA via `--training` when analyzing your real data.

### Overview

The simulation pipeline:

1. Reads your founder genotypes (VCF) and kinship estimates (KING output)
2. Selects unrelated founders from your data based on a kinship threshold
3. Simulates pedigrees using ped-sim with those founders
4. Runs IBD detection on the simulated genotypes (via a user-provided bash script)
5. Trains PONDEROSA classifiers on the simulated IBD sharing patterns
6. Outputs a `.classif.pkl` file for use with `--training`

### Prerequisites

In addition to the standard PONDEROSA dependencies, the simulation pipeline requires:

- [Ped-sim](https://github.com/williamslab/ped-sim) compiled and available locally
- A KING kinship output file (`.seg` format) for your founder samples
- A VCF file containing founder genotypes
- An IBD calling tool (e.g., [hap-ibd](https://github.com/browning-lab/hap-ibd)) and a bash script to run it (see below)

### Running a Simulation

```bash
python -m ponderosa simulate simulation_config.yaml
```

### Simulation Configuration

Create a YAML file with the following structure:

```yaml
# Required: ped-sim configuration
pedsim:
  pedsim_path: "/path/to/ped-sim"           # Directory containing ped-sim executable
  vcf_file: "/path/to/founders.vcf.gz"      # Input VCF with founder genotypes
  # Optional ped-sim settings (auto-detected if not provided):
  # simmap_file: "/path/to/refined_mf.simmap"
  # interference_file: "/path/to/nu_p_campbell.tsv"
  # def_file: "/path/to/custom_pedigrees.def"
  # random_seed: 42

# Required: KING output for founder selection
king_file: "/path/to/king.seg"

# Optional: training parameters
training:
  n_pairs_per_relationship: 100   # Pairs to simulate per relationship type (default: 100)
  max_kinship: 0.05               # Max kinship between founders in a family (default: 0.05)

# Required: path to your IBD caller bash script
# IMPORTANT: The script name (minus .sh) must match the IBD format name
# that PONDEROSA recognizes (e.g., "hap-ibd.sh" or "phasedibd.sh")
ibd_caller: "/path/to/hap-ibd.sh"

# Output prefix for the trained classifier pickle
# The simulation will produce: {output_path}.classif.pkl
output_path: "my_simulation_output"

# Whether to clean up temporary files after simulation (default: true)
cleanup_temp: true
```

### Simulation Configuration Reference

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `pedsim.pedsim_path` | Path | Yes | — | Directory containing the ped-sim executable |
| `pedsim.vcf_file` | Path | Yes | — | VCF file with founder genotypes (`.vcf` or `.vcf.gz`) |
| `pedsim.simmap_file` | Path | No | Auto-detected | Sex-specific recombination map |
| `pedsim.interference_file` | Path | No | Auto-detected | Crossover interference parameters |
| `pedsim.def_file` | Path | No | Built-in default | Pedigree definition file |
| `pedsim.random_seed` | Int | No | Random | Seed for reproducible simulations |
| `king_file` | Path | Yes | — | KING `.seg` output for founder kinship |
| `training.n_pairs_per_relationship` | Int | No | 100 | Simulated pairs per relationship type |
| `training.max_kinship` | Float | No | 0.05 | Maximum kinship allowed between founders |
| `ibd_caller` | Path | Yes | — | Path to IBD caller bash script |
| `output_path` | String | No | `"ponderosa_simulation"` | Output prefix for classifier pickle |
| `cleanup_temp` | Bool | No | `true` | Delete temporary files when done |

### Writing an IBD Caller Script

The simulation pipeline calls a user-provided bash script to run IBD detection on the simulated VCF. This keeps PONDEROSA agnostic to the IBD calling tool and its dependencies — your script handles everything (jar paths, map files, memory, etc.).

**Naming convention:** The script filename (minus `.sh`) must match the IBD caller format that PONDEROSA recognizes. For example, name your script `hap-ibd.sh` if it wraps hap-ibd, or `phasedibd.sh` if it wraps phasedibd.

**Contract:** PONDEROSA calls `bash your_script.sh <temp_dir>`, where `<temp_dir>` is the simulation workspace. Your script should:

1. Find the simulated VCF at `<temp_dir>/simulation.vcf` or `<temp_dir>/simulation.vcf.gz`
2. Write IBD output with the prefix `<temp_dir>/simulation` (e.g., `simulation.ibd.gz`)
3. Write a genetic map file as `<temp_dir>/genetic.map` in plink format (chr, rsid, cM, bp)
4. Exit with a non-zero code on failure

**Example `hap-ibd.sh` template:**

```bash
#!/bin/bash

# =========== Do not change this code ===========
# Script name must match the IBD caller format (e.g. hap-ibd.sh, phasedibd.sh)
# PONDEROSA uses the filename (minus .sh) to determine the IBD output format.

dirpath=${1}
prefix=${dirpath}/simulation

# Find the simulated VCF (ped-sim may output .vcf or .vcf.gz)
if [[ -f "${prefix}.vcf.gz" ]]; then
    vcffile="${prefix}.vcf.gz"
elif [[ -f "${prefix}.vcf" ]]; then
    vcffile="${prefix}.vcf"
else
    echo "ERROR: No simulated VCF found at ${prefix}.vcf or ${prefix}.vcf.gz"
    exit 1
fi
# ================================================

# ========== Edit these to your config ==========

basedir="/path/to/your/tools"
hapibd_jar="$basedir/hap-ibd.jar"

chr_map_dir="$basedir/plink.GRCh37.map"

# Memory in gigs
mem=8

# No. of cpus to use
nthreads=1

# ================================================

# ======== Checking that files exist =============
if [[ ! -f $hapibd_jar ]]; then
    echo "ERROR: hap-ibd jar not found: $hapibd_jar"
    exit 1
fi

if [[ ! -d $chr_map_dir ]]; then
    echo "ERROR: map directory not found: $chr_map_dir"
    exit 1
fi
# ================================================

# Create the concatenated map file
# PONDEROSA expects this file as "genetic.map" in plink format (chr, rsid, cM, bp)
cat "${chr_map_dir}/plink.chr1.GRCh37.map" > "${dirpath}/genetic.map"
for i in {2..22}; do
    cat "${chr_map_dir}/plink.chr${i}.GRCh37.map" >> "${dirpath}/genetic.map"
done

java -Xmx${mem}g -jar $hapibd_jar \
    gt=$vcffile \
    map=${dirpath}/genetic.map \
    nthreads=$nthreads \
    out=$prefix
```

### Using the Trained Classifiers

After the simulation completes, use the output pickle with PONDEROSA on your real data:

```bash
python -m ponderosa \
    --ibd your_real_data.ibd \
    --fam your_samples.fam \
    --ibd-caller hap-ibd \
    --map your_genetic.map \
    --training my_simulation_output.classif.pkl \
    --output my_results
```

When `--training` is provided, PONDEROSA skips classifier training and uses the pre-trained models from the simulation. This is useful when your study population differs from the default training data, such as in endogamous populations where IBD sharing patterns may not match standard expectations.