# PONDEROSA

### **Introduction**  
PONDEROSA (_**P**arent **O**ffspri**N**g pe**D**igree Inf**E**rence **RO**bu**S**t to Endog**A**my_) is an algorithm designed to assist in pedigree construction. PONDEROSA works well in datasets with high-quality long-range phasing. We have found that this can be better achieved in endogamous populations. Even in datasets with poor phasing, PONDEROSA can still distinguish avuncular from grandparent-grandchildren and maternal half-siblings from paternal half-siblings and will work well in datasets with few half-siblings or datasets with a narrow age range (which may be unlikely to have avuncular or grandparent-grandchildren pairs). PONDEROSA works best in datasets with existing pedigree structure, which is necessary for training the machine-learning classifiers. PONDEROSA will work out this existing pedigree structure from tracing parent-offspring lineages; therefore, every parent-offspring pair as inferred by KING must be present in the .fam file.  

Please note that PONDEROSA is designed to _assist_ pedigree construction and further steps by the user are required to construct the pedigree. We hope to change this in future versions of PONDEROSA. For now, PONDEROSA largely infers relationships in a vacuum (i.e. without considering the context of the pedigree). Relationship inference should be double-checked against the existing pedigree structure.  

### **Requirements**
**python3 or higher**  
**scikit-learn** and its dependencies, including numpy and pandas. We recommend running python directly from anaconda3, which has all the packages needed to run PONDEROSA.   

### **Running PONDEROSA**. 
Running PONDEROSA from the command line is easy, requiring only a parameter file.    
`python3.6 PONDEROSA.py [par_file]`
### **Parameter file**
The template for the parameter file is provided (**par_file.txt**).  
##### _Run type_ 
PONDEROSA has three different run types. Only one can be True; the other two must be False.  
| Run type | Description |
| -------- | ----------- |
| **po_only** | If selected, PONDEROSA will compute haplotype scores for PO pairs. Using age first, and then haplotype scores (if age is unavailable), this run type will output all PO pairs oriented as parent-child. We suggest running this step to create the .fam file necessary for other run types. |  
| **ped_only** | PONDEROSA will output all pairwise relationships present in the .fam file provided.|  
| **run_all** | Will run the entirety of PONDEROSA.|

##### _Required files_  
Each run type uses and requires different files. See below.  
| Run type | Required files | Optional files |
| -------- | -------------- | -------------- |
| **po_only** | king_file, map_file, match_file | ped_file, age_file |   
| **ped_only** | king_file, fam_file | |  
| **run_all** | king_file, map_file, match_file, fam_file | ped_file, age_file, hap_file |

##### _File descriptions_
The following files can be used by PONDEROSA. They must be formatted correctly; see sample files provided in Sample/.
| Flag | Description |
| ---- | ----------- |
|**king_file** | KING .seg file (or any .seg-formatted IBD file).
|**map_file** | PLINK-formatted .map file. This should be the full path/file name for chromosome 1. The map files for the other 21 autosomes must be in the same directory. Note that PONDEROSA expects a .map file for each chromosome but only one **--map** flag (see example par file). This .map file **must** be the same .map file used to generate IBD segments.
|**fam_file** | PLINK-formatted .fam file. All PO present in the KING file must be present in the .fam file. If age data is unavailable/unreliable and the parent/offspring cannot be distinguished in the pair, PONDEROSA can be run and the haplotype scores of the individuals can be used to make the distinction. 
|**match_file** | GERMLINE-formatted match file for chromosome 1. Again, PONDEROSA expects a .match file for each chromosome but only one **--match** flag (see example par file). If GERMLINE file, must be generated with GERMLINE’s --haploid flag (we suggest GERMLINE v1.5.3). iLASH .match files can also be used and will be detected by PONDEROSA.
|**ped_file** | PLINK-formatted .ped file used by PONDEROSA to stitch IBD segments together. If no .ped file is supplied, PONDEROSA stitches together two segments that are within 1 cM of each other. If .ped file is supplied, PONDEROSA only stitches two segments that are within 1 cM (can be changed with **--cm_gap** flag) of each other and have, at most, one discordant homozygote (can be changed with **--disc_homoz** flag). This flag can add considerable computational time. |
|**age_file** | Age file where the first column corresponds to the individual ID and the second column corresponds to the age. Note that not all individuals need an age. |
|**hap_file** | If PONDEROSA has already been run (either with **po_only** or **run_all**), supplying the haplotype score file here will skip the haplotype score calculation step. Will drastically reduce computation time. |

##### _Parameters_  
| Flag | Description |
| ---- | ----------- |
|**out** | Output file prefix.|
|**num_chr** | Number of autosomes.|
|**cm_gap** | Maximum gap in cM between IBD segments for them to be considered a single segment (see **ped_file** for more detail).|
|**disc_homoz** | Maximum number of discordant homozygotes between two IBD segments in order for them to be considered the same IBD segment (see **ped_file** for more detail). Only use if **ped file** is provided.|  
|**likelihood** | Minimum likelihood (0.5 - 1) required for a pair to be inferred as a 2nd degree pair. We recommend being more conservative here.|  
|**mhs_gap** | Maximum age-gap for maternal half-siblings. If you do not want PONDEROSA to consider age here, use an arbitrarily large age gap (e.g. 100).|
|**po_gap** | Minimum age-gap for parent-offspring. If you do not want PONDEROSA to consider age here, use 0 for this flag.|
|**gp_gap** | Minimum age-gap for a grandparent-grandchild pair. Note that if you do not want PONDEROSA to consider age, use 0 for this flag.|
|**trust_fs** | If True, PONDEROSA will assume that all KING-inferred FS with IBD2 > 0.15 are true FS. Recommended when pedigree data is sparse. |

