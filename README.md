# PONDEROSA

### **Introduction**  
PONDEROSA (_**P**arent **O**ffspri**N**g pe**D**igree Inf**E**rence **RO**bu**S**t to endog**A**my_) is an algorithm designed to assist in pedigree construction. PONDEROSA is highly sensitive to phase quality and will see reduced performance in datasets with any of the following: admixture, sparse pedigrees (especially sparse parent-offspring data), samples from different populations, etc. While we have found PONDEROSA to work best in endogamous populations, PONDEROSA can work well in any population as long as high phase quality can be achieved. However, running PONDEROSA with a dataset with suboptimal phase quality should only affect avuncular/grandparent-grandchild vs. half-sibling distinction and may still be useful for datasets with few half-siblings or datasets with a narrow age range (which may be unlikely to have avuncular or grandparent-grandchildren pairs).  

Please note that PONDEROSA is designed to _assist_ pedigree construction and further steps by the user are required to construct the pedigree. We hope to change this in future versions of PONDEROSA. For now, PONDEROSA largely infers relationships in a vacuum (i.e. without considering the context of the pedigree). Relationship inference should be double-checked against the existing pedigree structure.  

### **Requirements**
**python3.6**  
**scikit-learn** and its dependencies (install here: https://scikit-learn.org/stable/install.html)  

### **Input**  
##### _Required arguments_  
| Flag | Description |
| ---- | ----------- |
|**king** | KING .seg file (or any .seg-formatted IBD file). |
|**map** | PLINK-formatted .map file. The chromosome number should be replaced with “%s”. This .map file must be the same .map file used to generate IBD segments. |
|**fam** | PLINK-formatted .fam file. All PO present in the KING file must be present in the .fam file. If age data is unavailable/unreliable and the parent/offspring cannot be distinguished in the pair, PONDEROSA can be run and the haplotype scores of the individuals can be used to make the distinction. |
|**match** | GERMLINE-formatted match file where the chromosome number is replaced with “%s”. If GERMLINE file, must be generated with GERMLINE’s --haploid flag (we suggest GERMLINE v1.5.3). iLASH .match files can also be used, but PONDEROSA’s **\-\-ilash** flag must be used.|

##### _Optional arguments_  
| Flag | Description |
| ---- | ----------- |
|**out** | Output file prefix. _Default: “PONDEROSA”_ |
|**ilash** | For use if .match file is in iLASH format. |
|**haps** | If PONDEROSA has already been run, supplying the haplotype score file here will skip the haplotype score calculation step. |
|**age** | Age file where the first column corresponds to the individual ID and the second column corresponds to the age. Note that not all individuals need an age. |
|**gp_age** | Minimum age-gap for a grandparent-grandchild pair. _Default: 30_ |
|**mhs_age** | Maximum age-gap for maternal half-siblings. _Default: 30_ |
|**ped** | PLINK-formatted .ped file used by PONDEROSA to stitch IBD segments together. If no .ped file is supplied, PONDEROSA stitches together two segments that are within 1 cM of each other. If .ped file is supplied, PONDEROSA only stitches two segments that are within 1 cM of each other _and_ have, at most, one discordant homozygote.|  

##### _Example_
`python3.6 PONDEROSA.py --map Himba_chr%s.map --ped Himba_chr%s.ped --king king.seg –fam Himba.fam --match Himba_chr%s.match --out Himba`

