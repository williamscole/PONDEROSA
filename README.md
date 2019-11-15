# PONDEROSA

### **Introduction**  
PONDEROSA (_**P**arent **O**ffspri**N**g pe**D**igree Inf**E**rence **RO**bu**S**t to endog**A**my_) is an algorithm designed to assist in pedigree construction. PONDEROSA is highly sensitive to phase quality and will see reduced performance in datasets with any of the following: admixture, sparse pedigrees (especially sparse parent-offspring data), samples from different populations, etc. While we have found PONDEROSA to work best in endogamous populations, PONDEROSA can work well in any population as long as high phase quality can be achieved. However, running PONDEROSA with a dataset with suboptimal phase quality should only affect avuncular/grandparent-grandchild vs. half-sibling distinction and may still be useful for datasets with few half-siblings or datasets with a narrow age range (which may be unlikely to have avuncular or grandparent-grandchildren pairs).  

Please note that PONDEROSA is designed to assist pedigree construction and further steps by the user are required to construct the pedigree. We hope to change this in future versions of PONDEROSA. For now, PONDEROSA largely infers relationships in a vacuum (i.e. without considering the context of the pedigree). Relationship inference should be double-checked against the existing pedigree structure.  

### **Requirements**
**python3.6**  
**scikit-learn** and its dependencies (install here: https://scikit-learn.org/stable/install.html)