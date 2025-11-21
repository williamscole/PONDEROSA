# PONDEROSA

### Setting up the virtual environment and installing phasedibd

*I have tested this with Python 3.10.9. Earlier versions of Python may work, but phasedibd needs at least Python 3.8 to install.*

1. Create the virtual environment: `python -m venv ponderosa_venv`
2. Activate the virtual environment: `source ponderosa_venv/bin/activate`
3. Install required packages: `python -m pip install -r requirements.txt` (`requirements.txt` is found in this directory)
4. Clone `phasedibd`: `git clone https://github.com/23andMe/phasedibd`
5. `cd phasedibd`
6. `make`**
7. `python setup.py install`
8. `python tests/unit_tests.py`

**On the HPC that I use, I need to replace the third line of the `phasedibd` Makefile with `python -m pip install cython` (currently it reads `pip install cython --user`, which may cause permissions issues).

### Setting up a conda environment

1. `conda env create -f environment.yaml` will create a conda environment called `ponderosa` using Python 3.10.9 and will install all necessary dependencies.
2. Activate the virtual environment: `conda activate ponderosa`
3. To set up `phasedibd`, proceed using steps 4-8 above.

### Creating a map file

Running `python pedigree_tools.py -interpolate [options]` will create a Plink-formatted MAP file with the genetic position filled in. It takes as input a genetic map and interpolates the cM coordinate of a physical position.

The following arguments are accepted:
- `-input_map`: the MAP file in which the genetic position needs to be filled in (no header; assumes the columns are: chromosome, SNP ID, cM position filled with 0, physical position)
- `-genetic_map`: the file containing genetic map coordinates (no header; requires the following columns: chromosome, cM position, physical position)
- `-columns`: default is to assume that `genetic_map` is a Plink MAP-formatted file. Alternatively, you can use this optional argument to specify the 0-indexed index of the chromosome, cM position, and physical position columns (in that order). E.g., `-columns 2 3 4` indicates that the chromosome, cM, and Mb columns are the 3rd, 4th, and 5th columns, respectively.
- `-sites`: if only a subset of loci are to be output. This file has no headers and expects two columns: chromosome and physical position. This is **optional** and omitting it will include all sites in `input_map`.

The output is the file: `[input_map]_interpolated.map`

### Running phasedibd

Running `python pedigree_tools.py -phasedibd [options]` will run `phasedibd`.

The following arguments are accepted:
- `-input_vcf`: uncompressed VCF for a single chromosome.
- `-input_map`: Plink MAP file. `phasedibd` requires that it contains *exactly* the sites in `input_vcf`; see above `-interpolate` to create this file.
- `-outfile`: the name of the file to write the IBD segments to.

### Running PONDEROSA

*Required file arguments*
- `--ibd`: phasedibd IBD output file. If all chromosomes are in the same file, simply provide the file name. Otherwise, provide the file name for chromosome 1 (assumes that "chr1" is in the file name).
- `--fam`: PLINK-formatted FAM file.
- `--king`: KING .seg file.

*Optional file arguments*
- `--ages`: Age file where the first column is the individual ID and the second column is the age.
- `--map`: PLINK-formatted MAP file (*highly recommended*). Sites from all chromosomes can be in the same file; otherwise provide the file for chromosome 1.
- `--populations`: For running Ponderosa on a subset of samples, this file contains individual IDs (column 1) and a population label (column 2); you may specify which population you'd like run with `--population`.
- `--yaml`: You may provide all of these arguments in a YAML file instead of on the command-line interface.
- `--pedigree_codes`: Instructs PONDEROSA which relationships to look for in the dataset.

*Other arguments*
- `--output`: Output file prefix. *Default: Ponderosa*
- `--min_p`: Minimum probability required for the relationship output. E.g., if P(2nd)=0.98 and P(MHS)=0.93, setting `min_p=0.95` would report the pair as `2nd`, but `min_p=0.9` would report the pair as `MHS`. *Default 0.50*
- `--population`: Used with `--populations`, specifies the population to run PONDEROSA on.
- `--assess`: For assessing the performance of Ponderosa on the known pairs.
- `--training`: Full path and file name to the degree classifier. Assumes that the the haplotype score classifier (hap) and number of IBD segments classifier (nsegs) have the same prefix/suffix.  

### Plotting IBD segments

```
from pedigree_tools import Karyogram

# either a list of .map files or a single map file that contains all chromosomes
map_files = [f"chr_{chrom} for chrom in range(1, 23)]

# initialize the object, set cm to True if you're plotting centimogran positions, cm = False plots Mb positions
kgram = Karyogram(map_files, cm = True)

# segments is a list of segments with the following info: [chromosome, start, stop, haplotype index (0 or 1)]
segments = [[1, 32.1, 45.6, 0], [2, 45.5, 123.4, 1]]

# plot segments; optional keyword arguments include
# file_name [default: "karyogram"]: writes the output as [file_name].png
# hap0_color [default: "cornflowerblue"]: the color of IBD segments on the 0 haplotype
# hap1_color [default: "tomato"]: the color of IBD segments on the 1 haplotype
# dpi [default: 500]: dpi of the plot

kgram.plot_segments(segments, file_name = "my_karyogram", hap0_color = "skyblue")
```
