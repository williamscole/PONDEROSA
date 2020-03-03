# PONDEROSA

### **Introduction**  
PONDEROSA (_**P**arent **O**ffspri**N**g pe**D**igree Inf**E**rence **RO**bu**S**t to Endog**A**my_) is an algorithm designed to assist in pedigree construction. PONDEROSA works well in datasets with high-quality long-range phasing. We have found that this can be better achieved in endogamous populations. Even in datasets with poor phasing, PONDEROSA can still distinguish avuncular from grandparent-grandchildren and maternal half-siblings from paternal half-siblings and will work well in datasets with few half-siblings or datasets with a narrow age range (which may be unlikely to have avuncular or grandparent-grandchildren pairs). PONDEROSA works best in datasets with existing pedigree structure, which is necessary for training the machine-learning classifiers. PONDEROSA will work out this existing pedigree structure from tracing parent-offspring lineages; therefore, every parent-offspring pair as inferred by KING must be present in the .fam file.  

Please note that PONDEROSA is designed to _assist_ pedigree construction and further steps by the user are required to construct the pedigree. We hope to change this in future versions of PONDEROSA. For now, PONDEROSA largely infers relationships in a vacuum (i.e. without considering the context of the pedigree). Relationship inference should be double-checked against the existing pedigree structure.  

### **Requirements**
**python3 or higher**  
**scikit-learn** and its dependencies (install here: https://scikit-learn.org/stable/install.html). Alternatively, we recommend running python directly from anaconda3, which has all the packages needed to run PONDEROSA.    

### **Input**  
##### _Required arguments_  
| Flag | Description |
| ---- | ----------- |
|**king** | KING .seg file (or any .seg-formatted IBD file). |
|**map** | PLINK-formatted .map file. The chromosome number should be replaced with “%s”. Note that PONDEROSA expects a .map file for each chromosome but only one **--map** flag (see example script). This .map file must be the same .map file used to generate IBD segments. |
|**fam** | PLINK-formatted .fam file. All PO present in the KING file must be present in the .fam file. If age data is unavailable/unreliable and the parent/offspring cannot be distinguished in the pair, PONDEROSA can be run and the haplotype scores of the individuals can be used to make the distinction. |
|**match** | GERMLINE-formatted match file where the chromosome number is replaced with “%s”. Again, PONDEROSA expects a .match file for each chromosome but only one **--match** flag (see example script). If GERMLINE file, must be generated with GERMLINE’s --haploid flag (we suggest GERMLINE v1.5.3). iLASH .match files can also be used, but PONDEROSA’s **\-\-ilash** flag must be used.|

##### _Optional arguments_  
| Flag | Description |
| ---- | ----------- |
|**out** | Output file prefix. _Default: “PONDEROSA”_ |
|**chr** | Number of autosomes. Change only for non-human samples. _Default: 22_ |
|**ilash** | For use if .match file is in iLASH format. |
|**haps** | If PONDEROSA has already been run, supplying the haplotype score file here will skip the haplotype score calculation step. |
|**second_train** | If PONDEROSA has already been run, the user can supply the .training file from that run to train the 2nd degree classifier of another PONDEROSA run (with a different dataset that has few 2nd degree training pairs). |
|**age** | Age file where the first column corresponds to the individual ID and the second column corresponds to the age. Note that not all individuals need an age. |
|**gp_gap** | Minimum age-gap for a grandparent-grandchild pair. Note that if you do not want PONDEROSA to consider age, use 0 for this flag. _Default: 30_ |
|**mhs_gap** | Maximum age-gap for maternal half-siblings. If you do not want PONDEROSA to consider age here, use an arbitrarily large age gap (e.g. 100). _Default: 30_ |
|**po_gap** | Minimum age-gap for parent-offspring. If you do not want PONDEROSA to consider age here, use 0 for this flag. _Default: 15_ |
|**ped** | PLINK-formatted .ped file used by PONDEROSA to stitch IBD segments together. If no .ped file is supplied, PONDEROSA stitches together two segments that are within 1 cM of each other. If .ped file is supplied, PONDEROSA only stitches two segments that are within 1 cM (can be changed with **--cm_gap** flag) of each other and have, at most, one discordant homozygote (can be changed with **--disc_homoz** flag).|
|**cm_gap** | Maximum gap in cM between IBD segments for them to be considered a single segment (see **--ped** flag for more detail). _Default: 1_ |
|**disc_homoz** | Maximum number of discordant homozygotes between two IBD segments in order for them to be considered the same IBD segment (see **--ped** flag for more detail). Only use if **--ped file** is used. _Default: 1_ |  
|**likelihood** | Minimum likelihood (0.5 - 1) required for a pair to be inferred as a 2nd degree pair. We recommend being more conservative here. _Default: 0.80_ |  

##### _Example_
`python3.6 PONDEROSA.py --map Himba_chr%s.map --ped Himba_chr%s.ped --king king.seg –fam Himba.fam --match Himba_chr%s.match --out Himba`

