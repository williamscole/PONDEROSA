# PONDEROSA

PONDEROSA (Parent OffspriNg peDigree infErence RObuSt to endogAmy) is a Python tool for inferring genetic relationships between individuals using Identity By Descent (IBD) segments. The tool uses machine learning classifiers trained on IBD sharing patterns to distinguish between different degrees of biological relationships.

## Overview

PONDEROSA analyzes IBD segments to infer relationships between pairs of individuals. It supports multiple IBD calling software outputs and provides flexible configuration through both command-line arguments and YAML configuration files.

## Installation

### Prerequisites

PONDEROSA requires Python 3.8+ and the following dependencies:
- pandas
- polars  
- numpy
- scikit-learn
- pyyaml
- networkx
- matplotlib
- seaborn

### Environment Setup

The recommended installation method is using conda with the provided environment file:

```bash
conda env create -f environment.yml
conda activate ponderosa_v2
```

## Usage

### Basic Command Line Usage

```bash
python -m ponderosa [options]
```

### Using Configuration Files

PONDEROSA supports YAML configuration files for complex analyses:

```bash
python -m ponderosa --config config.yaml
```

## Command Line Arguments

### File Arguments

These arguments specify input and output files:

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `--config` | Path | No | YAML configuration file |
| `--ibd` | Path | Yes | IBD segments file |
| `--fam` | Path | Yes | PLINK FAM file with individual information |
| `--ibd-caller` | Choice | Yes | IBD calling software: `phasedibd`, `hap-ibd` |
| `--map` | Path | Yes* | Genetic map file for coordinate conversion |
| `--ages` | Path | No | File containing age information for individuals |
| `--priors` | Path | No | File specifying relationship priors (e.g., age-based priors) |
| `--populations` | Path | No | Population assignment file *(in testing — not yet available)* |
| `--training` | Path | No | Directory containing pre-trained models |

*Unless your IBD caller outputs the segment length in cM


### Algorithm Arguments

These control the relationship inference algorithm:

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--min-segment-length` | Float | 3.0 | Minimum IBD segment length in centiMorgans (cM) |
| `--min-total-ibd` | Float | 50.0 | Minimum total IBD sharing in cM for a pair to be analyzed |
| `--max-gap` | Float | 1.0 | Maximum gap in cM for stitching adjacent segments |
| `--population` | String | "pop1" | Population identifier for analysis |
| `--genome-length` | Float | 3545.0 | Total genome length in cM |

### Output Arguments

These control output format and verbosity:

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--output` | String | "ponderosa_results" | Output file prefix |
| `--min-probability` | Float | 0.5 | Minimum probability threshold for reporting relationships |
| `--create-plots` | Flag | False | Generate visualization plots *(in testing — not yet available)* |
| `--verbose`, `-v` | Count | 0 | Increase verbosity (can be used multiple times: `-v`, `-vv`, `-vvv`) |
| `--write-training` | Flag | False | Write trained classifiers to pickle file |
| `--debug` | Flag | False | Show full error tracebacks |

## Configuration File Format

PONDEROSA supports YAML configuration files with the following structure:

```yaml
# File inputs
files:
  ibd: "path/to/ibd_segments.txt"
  fam: "path/to/individuals.fam"
  ibd_caller: "hap-ibd"
  ages: "path/to/ages.txt"                    # Optional
  mapf: "path/to/genetic.map"                 # Optional
  priors: "path/to/priors.yaml"              # Optional
  populations: "path/to/populations.txt"      # Optional (in testing — not yet available)
  training: "path/to/trained_models/"         # Optional

# Algorithm parameters  
algorithm:
  min_segment_length: 3.0
  min_total_ibd: 50.0
  max_gap: 1.0
  population: "pop1"
  genome_length: 3545.0

# Output settings
output:
  output: "my_analysis_results"
  min_probability: 0.5
  write_readable: true
  verbose: true
  write_training: false
```

## Simulation Pipeline

*(In testing — not yet available)*

PONDEROSA includes a simulation pipeline for generating training data using [Ped-sim](https://github.com/williamslab/ped-sim). The simulation workflow simulates relative pairs from real genotype data and trains classifiers on the resulting IBD sharing patterns.

### Simulation Configuration

```yaml
pedsim:
  pedsim_path: "/path/to/ped-sim"
  vcf_file: "/path/to/founders.vcf.gz"

king_file: "/path/to/king.seg"

training:
  n_pairs_per_relationship: 100
  max_kinship: 0.05

# Optional: processing settings
ibd_caller: "hap-ibd.sh"
output_path: "simulation_output"
cleanup_temp: true
```

## Input File Formats

### IBD Segments File

The format depends on the IBD caller specified:

#### phasedibd Format
```
id1    id2    chromosome    start_cm    end_cm    id1_haplotype    id2_haplotype
IND1   IND2   1             10.5        25.3      0                0
IND1   IND3   1             30.1        45.7      1                0
```

#### hap-ibd Format  
```
id1    id1_haplotype    id2    id2_haplotype    chromosome    start_bp    end_bp    length_cm
IND1   1                IND2   2                1             1000000     2500000   15.2
IND1   2                IND3   1                1             3000000     4200000   12.8
IND2   1                IND4   2                2             5000000     7800000   18.5
```

*Actual hap-ibd file should have no header

### FAM File (PLINK Format)
```
fam_id  ind_id  father  mother  sex  phenotype
FAM1    IND1    0       0       1    -9
FAM1    IND2    0       0       2    -9
FAM2    IND3    0       0       1    -9

```
*The actual fam file should have no header

### Ages File (Optional)
```
individual_id    age or year of birth
IND1            45
IND2            67
IND3            32
```
*The actual age file should have no header

### Genetic Map File (Optional, for hap-ibd)
```
chromosome    position_bp    position_cm
1             1000000        0.5
1             2000000        1.2
1             3000000        1.8
```
*The actual map file should have no header

### Priors File (Optional)
The priors file allows you to specify age-based relationship constraints. The format is a tab-separated or space-separated file with three columns:

```
rel    operator    age_gap
MHS    >           25
GP     <=          30
```

**Column Descriptions:**
- `rel`: Relationship abbreviation (e.g., MHS for maternal half-siblings, GP for grandparent-grandchild)
- `operator`: Comparison operator (`>`, `<`, `=`, `>=`, `<=`)  
- `age_gap`: Age difference threshold in years

**Example Usage:**
In this example, if two 2nd degree individuals have a >25 year age gap, P(MHS) would be set to 0 and the other probabilities rescaled.

## Output Files

PONDEROSA generates several output files with the specified prefix:

- `{prefix}_pairs.txt`: Processed pair-wise IBD statistics and relationship predictions
- `{prefix}.mhier.pkl`: Relationship-likelihood tree data structure (binary pickle file)
- `{prefix}.classif.pkl`: Trained classifiers (binary pickle file)
- `{prefix}_evaluation.txt`: Evaluation report (if truth file is provided)
- `{prefix}_plots/`: Visualization directory (if `--create-plots` specified) *(in testing — not yet available)*

## Relationship Categories

PONDEROSA classifies relationships into hierarchical categories:

- **1st Degree**: Parent-Child (PC), Full Siblings (FS)
- **2nd Degree**: Grandparent-Grandchild (PGP/MGP), Aunt/Uncle-Niece/Nephew (AV), Half-Siblings (PHS/MHS)
- **3rd Degree**
- **Unrelated**: Distant or no biological relationship

## Performance Considerations

- **Single vs Multiple Files**: Single IBD files are more memory-efficient as they allow Polars to optimize filtering during file reading. Multiple IBD files require loading all data into memory before applying genome-wide filters like `--min-total-ibd`.
- **Memory Usage**: Large datasets may require substantial RAM. Consider filtering pairs with `--min-total-ibd` to reduce memory requirements.
- **Processing Time**: Runtime scales quadratically with the number of individuals. Use `--verbose` to monitor progress.
- **Accuracy**: Results improve with higher quality IBD calls and appropriate parameter tuning for your dataset.

## Troubleshooting

### Common Issues

1. **File Not Found Errors**: Ensure all input files exist and paths are correct
2. **Memory Errors**: Increase `--min-total-ibd` to reduce the number of pairs analyzed
3. **No Results**: Lower `--min-probability` threshold or check IBD segment quality
4. **Format Errors**: Verify IBD file format matches the specified `--ibd-caller`

### Debug Mode

Use `--debug` to see full error tracebacks:

```bash
python -m ponderosa --debug --config config.yaml
```

### Verbose Output

Use multiple `-v` flags for detailed progress information:

```bash
python -m ponderosa -vvv --config config.yaml
```

## Citation

If you use PONDEROSA in your research, please cite:

Williams CM, Scelza BA, Slack SD, Font-Porterias N, Al-Hindi DR, Mathias RA, Watson H, Barnes KC, Lange E, Johnson RK, Gignoux CR, Ramachandran S, Henn BM. A rapid accurate approach to inferring pedigrees in endogamous populations. *Genetics*. 2025;230(4):iyaf094. doi: [10.1093/genetics/iyaf094](https://doi.org/10.1093/genetics/iyaf094)

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).

PONDEROSA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html) for more details.

## Support

For questions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com).