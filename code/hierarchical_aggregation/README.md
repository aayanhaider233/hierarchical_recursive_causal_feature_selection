# Hierarchical aggregation



This stage performs DMR detection, CpG-to-DMR aggregation, DMR-to-gene mapping, DMR-to-gene aggregation, and WGCNA.





### Input



* Batch corrected, M-scale CpG methylation matrix (\*\_cpg_mval\_matrix\_corrected.csv)
* Training metadata set (train\_metadata.csv)
* Illumina 450K array manifest





### Output



* data/intermediate/dmr\_groups.rds
* data/intermediate/dmr\_cpg\_map.csv
* data/intermediate/dmr\_gene\_map.csv
* data/intermediate/gene\_methylation\_matrix.csv
* data/intermediate/gene\_module\_assignments.csv
* data/intermediate/gene\_module\_kme\_matrix.csv
* data/methylation/train\_dmr\_matrix.csv 
* data/methylation/test\_dmr\_matrix.csv
* data/methylation/test\_cpg\_mval\_corrected.csv
* data/causal\_methylation\_inputs/me\_matrix.csv





### Execution



Run the stage driver:



```bash

python run_aggregation.py

```



Or run the full project pipeline from code/:



```bash

python main.py

```





### Structure



|**Order**|**Module**|**Description**|
|-|-|-|
|01|dmr\_detection.R|Detect DMRs from CpG methylation and disease status data|
|02|dmr\_aggregation.R|Construct aggregated DMR methylation matrix from CpG methylation matrix|
|03|dmr\_to\_gene\_mapping.py|Construct DMR-gene map|
|04|gene\_aggregation.py|Construct aggregated gene methylation matrix from DMR methylation matrix|
|05|wgcna.R|Perform WGCNA to obtain ME matrix, gene-module assignments, and gene-module membership matrix|




### Workflow



1. DMR detection
2. CpG-to-DMR aggregation 
3. DMR-to-gene assignment
4. DMR-to-gene aggregation 
5. WGCNA
6. Save aggregated methylation datasets and relevant intermediate data





#### Notes



* DMR detection is performed on the training set only. The resulting DMR definitions are reused for aggregation of both the training and test sets to avoid information leakage.
* Only autosomal probes (chromosomes 1–22) are considered during DMR detection.
* CpG-to-DMR aggregation uses the mean M-value of all CpGs belonging to a detected DMR.
* DMR-to-gene mapping follows a hierarchical strategy: constituent CpG annotation → promoter overlap (TSS −2000 bp to +500 bp) → nearest gene based on the minimum distance between a TSS and the DMR's endpoints.
* DMR-to-gene aggregation uses the mean methylation value across all DMRs assigned to the same gene.
* WGCNA is performed on the gene-level methylation matrix using an unsigned network with a minimum module size of 30 and module merging at eigengene correlation r ≥ 0.75 (mergeCutHeight = 0.25).
* The final ME matrix is intended as the primary methylation input for downstream causal discovery and causal effect estimation.


