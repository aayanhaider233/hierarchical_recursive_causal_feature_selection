# Preprocessing



This stage performs data harmonisation, quality control, epigenetic age estimation, probe filtering, dataset concatenation, M-value conversion, batch effect correction, train-test partitioning, and epigenetic age acceleration computation.





### Input



* Raw metadata files (\*\_series\_matrix.txt.gz)
* Raw methylation matrices (\*.csv.gz)
* Illumina annotation files in annotations/methylation/





### Output



* data/intermediate/batch\_info.csv 
* data/methylation/cpg\_mval\_matrix\_full\_corrected.csv 
* data/methylation/train\_cpg\_mval\_corrected.csv
* data/methylation/test\_cpg\_mval\_corrected.csv
* data/metadata/train\_metadata.csv
* data/metadata/test\_metadata.csv`





### Execution



Run the stage driver:



```bash

python run_preprocessing.py

```



Or run the full project pipeline from code/:



```bash

python main.py

```





### Structure



|**Order**|**Module**|**Description**|
|-|-|-|
|01|metadata\_harmonisation.py|Parse and harmonise GEO metadata|
|02|probe\_quality\_control.py|Initial probe filter, QC and, detection p-value filtering|
|03|epigenetic\_age.py|Generate epigenetic age estimates (Horvath, Hannum, and GrimAge)|
|04|probe\_filter.py|Remove sex-chromosome, cross-reactive, and multi-mapping probes|
|05|concatenate\_datasets.py|Merge cohorts and generate batch information|
|06|methylation\_beta\_to\_M.py|Convert methylation beta-values to M-values|
|07|batch\_effect\_correction.py|Batch effect correction|
|08|partition\_data.py|Stratified train-test split|
|09|epigenetic\_age\_acceleration.py|Compute epigenetic age acceleration (EAA)|
|10|run\_preprocessing.py|Execute the complete preprocessing workflow|





### Workflow



1. Metadata harmonization
2. Non-CpG probe removal and quality control
3. Epigenetic age estimation
4. Annotation-based probe filtering
5. Cohort concatenation and batch information generation
6. Beta-to-M value conversion
7. Batch effect correction
8. Train-test partitioning
9. Epigenetic age acceleration computation
10. Save processed datasets





#### Notes



* Epigenetic ages are computed before annotation-based probe filtering to preserve clock-specific CpG sites.
* Batch correction is performed on the full M-value matrix before train/test partitioning.
* EAA models are fitted on the training set only and applied to the test set to avoid information leakage.
* Epigenetic age features are replaced by their corresponding EAA features.


