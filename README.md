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
The template for the parameter file is provided (_par_file.txt_).  
#### _Run type_ 
PONDEROSA has three different run types. Only one can be True; the other two must be False.  
| Run type | Description |
| -------- | ----------- |
| **po_only** | If selected, PONDEROSA will compute haplotype scores for PO pairs. Using age first, and then haplotype scores (if age is unavailable), this run type will output all PO pairs oriented as parent-child. We suggest running this step to create the .fam file necessary for other run types. |  
| **ped_only** | PONDEROSA will output all pairwise relationships present in the .fam file provided.|  
| **run_all** | Will do the above but will also infer unresolved second degree relationships.|

#### _File requirements_  
Each run type uses and requires different files. See below.  
| Run type | Required files | Optional files |
| -------- | -------------- | -------------- |
| **po_only** | king_file, map_file, match_file | ped_file, age_file |   
| **ped_only** | king_file, fam_file | |  
| **run_all** | king_file, map_file, match_file, fam_file | ped_file, age_file, hap_file |

#### _File descriptions_
The following files can be used by PONDEROSA. They must be formatted correctly; see sample files provided in Sample/. The file name for optional files that are not supplied should be "None".
| Flag | Description |
| ---- | ----------- |
|**king_file** | KING .seg file (or any .seg-formatted IBD file).
|**map_file** | PLINK-formatted .map file. This should be the full path/file name for chromosome 1. The map files for the other 21 autosomes must be in the same directory. Note that PONDEROSA expects a .map file for each chromosome, but the user need only supply the name of the first chromosome (see example par file). This .map file **must** be the same .map file used to generate IBD segments.
|**fam_file** | PLINK-formatted .fam file. All PO present in the KING file should be present in the .fam file. If age data is unavailable/unreliable and the parent/offspring cannot be distinguished in the pair, PONDEROSA should be run with **po_only**, which will orient pairs into parent-offspring.
|**match_file** | GERMLINE-formatted match file for chromosome 1. Again, PONDEROSA expects a .match file for each chromosome but only one **--match** flag (see example par file). If GERMLINE file, must be generated with GERMLINE’s --haploid flag (we suggest GERMLINE v1.5.3). iLASH .match files can also be used and will be detected by PONDEROSA.
|**ped_file** | PLINK-formatted .ped file used by PONDEROSA to stitch IBD segments together. If no .ped file is supplied, PONDEROSA stitches together two segments that are within 1 cM of each other. If .ped file is supplied, PONDEROSA only stitches two segments that are within 1 cM (can be changed with **cm_gap**) of each other and have, at most, one discordant homozygote (can be changed with **disc_homoz**). This flag can add considerable computational time and is generally not recommdended. |
|**age_file** | Age file where the first column corresponds to the individual ID and the second column corresponds to the age. Note that not all individuals need an age, and vice versa. |
|**hap_file** | PONDEROSA will create a hap file if it has already been run with **po_only** or **run_all**. It can be supplied here, and PONDEROSA will skip the haplotype score calculation step. Will drastically reduce computation time. |

#### _Parameters_  
| Flag | Description |
| ---- | ----------- |
|**out** | Output file prefix.|
|**num_chr** | Number of autosomes.|
|**cm_gap** | Maximum gap in cM between IBD segments for them to be considered a single segment (see **ped_file** for more detail).|
|**disc_homoz** | Maximum number of discordant homozygotes between two IBD segments in order for them to be considered the same IBD segment (see **ped_file** for more detail). Will ignore if ped file is not provided. |  
|**likelihood** | Minimum likelihood (0.5 - 1) required for a pair to be inferred as a 2nd degree pair. We recommend being more conservative here.|  
|**mhs_gap** | Maximum age-gap for maternal half-siblings. If you do not want PONDEROSA to consider age here, use an arbitrarily large age gap (e.g. 100).|
|**po_gap** | Minimum age-gap for parent-offspring. If you do not want PONDEROSA to consider age here, use 0 for this flag.|
|**gp_gap** | Minimum age-gap for a grandparent-grandchild pair. Note that if you do not want PONDEROSA to consider age, use 0 for this flag.|
|**trust_fs** | If True, PONDEROSA will assume that all KING-inferred FS with IBD2 > 0.15 are true FS. Recommended when pedigree data is sparse. |

### **Output files.** 

#### _[out].log_ 

Provides information about the PONDEROSA run, including supplied parameters and files, run time, and any errors.  

##### _Error messages_
| Error code | Description |
| ---------- | ----------- |
| 01 | PONDEROSA is attempting to assign putative 2nd degree relatives to a pedigree relationship, but there are not enough training pairs of either AV, GP, MHS, PHS. The dataset is too sparse to train the classifier. | 
| 02 | To maximize the number of relative pairs, all PO pairs present in KING should be present in the .fam file provided. Pairs here are present in KING but not in the .fam file. They should be added to the .fam, but PONDEROSA will continue running. |
| 03 | The following sets of individuals are full siblings but have different parents listed in the .fam file. PONDEROSA has ignored this, but the user should double check. For example, if individual A and B are FS but the father ID of A is different than the father ID of B. |
| 04 | PONDEROSA is attempting to classify ambiguous sibships with a classifier, but there are not enough FS or 2nd degree relationships to train the classifier. Try running again with **trust_fs** set True, which will skip this step. |
| 05 | KING has a low IBD2 threshold for pairs to be considered FS. Pairs listed here have abnormally low IBD2 values for a FS pair, but have been inferred as such by KING. These are unlikely to be real FS pairs, and PONDEROSA will ignore them. |

#### _[out]\_PO.txt_
Pairwise data for parent-offspring present in KING file.  
| Column | Description |
| ------ | ----------- |
| PAIR_ID | Unique pair ID for the PO pair |
| IID1 | IID for individual 1 |
| IID2 | IID for individual 2 |
| H1 | Haplotype score for individual 1 |
| H2 | Haplotype score for individual 2 |
| AGE1 | Age of individual 1 |
| AGE2 | Age of individual 2 |
| CHILD | IID of the inferred child |
| PARENT | IID of the inferred parent |
| METHOD | AGE if age data was used to orient; H if haplotype scores used to orient. |
| STRENGTH | If METHOD is AGE, the difference in age of the PO pair. If METHOD is H, the difference in haplotype scores of the PO pair; PO pairs with a small difference in haplotype scores (i.e. close to 0) cannot be as confidentally oriented. | 

#### _[out]\_pairs.txt_
Relative pairs present in the .fam file.  
| Column | Description |
| ------ | ----------- |
| PAIR_ID | Unique pair ID for the relative pair |
| IID1 | The IID of the genetically younger individual in a pair (if applicable). For example, a child in a PO pair or niece of an avuncular pair. |
| IID2 | the IID of the genetically older individual. |
| GTD | True if both individuals in pair are genotyped. |
| IBD1 | KING IBD1 value. |
| IBD2 | KING IBD2 value. |
| PI_HAT | Total proportion IBD shared (equal to 0.5*\IBD1 + IBD2). | 
| KINGINF | KING-inferred degree of relatedness | 
| REL | Pedigree relationship of the pair as inferred by PONDEROSA. |
| DEGREE | Degree of relatedness of the pedigree relationship. |

#### _[out]\_second.txt_
Pairwise data of putative second degree relatives.  
| Column | Description |
| ------ | ----------- |
| PAIR_ID | Unique pair ID for the relative pair |
| YOUNGER | The genetically younger individual in the pair, if applicable. | 
| OLDER | The genetically older individual in the pair. |
| METHOD | AGE if age data used to orient younger-older; H is haplotype scores used to orient younger-older. |
| REL | Inferred second degree relative type. |
| SECOND_PROB | Probability of pair being second degree related. |
| PROB | Probability of the inferred second degree relative type. |
| HSR | Haplotype score ratio of the pair. |
| N | Number of IBD segments shared. |
| AV | Probability of pair being AV. |
| GP | Probability of pair being GP. |
| MHS | Probability of pair being MHS. |
| PHS | Probability of pair being PHS. |
| AV_ERROR | True if 1) Pair is inferred as AV _and_ 2) Age data disagrees with haplotype data (i.e. the older individual appears to be the niece/nephew). Note that it is possible that a niece/nephew is older than their uncle/aunt. |
