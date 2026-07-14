# **Metadata files**



**train\_metadata.csv** and **test\_metadata.csv** contain covariates used throughout the study.





### **Variables**



|**Name**|**Description**|
|-|-|
|sample\_id|Unique identifier constructed from the original barcode and batch pair.|
|age|Chronological age in years.|
|sex|0 = Male, 1 = Female|
|disease|0 = Control, 1 = Schizophrenia|
|horvath2013\_EAA|Horvath epigenetic age acceleration|
|hannum\_EAA|Hannum epigenetic age acceleration|
|grimage2\_EAA|GrimAge epigenetic age acceleration|





#### Notes



* Batch and barcode are not retained after sample\_id construction.
* All feature-selection methods were applied to the methylation features; covariates remained unchanged.

