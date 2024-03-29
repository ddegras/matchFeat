# matchFeat
Fast Algorithms for One-To-One Feature Matching

Statistical methods to match feature vectors between multiple datasets in a one-to-one fashion. Applications include object tracking, video surveillance, remote sensing as well as multilevel modeling. Given a fixed number of classes/distributions, for each dataset, exactly one vector of each class is observed without label. The goal is to label the feature vectors using each label exactly once so to produce the best match across datasets, e.g. by minimizing the variability within classes. Several statistical solutions based on empirical loss functions and probabilistic modeling are provided.  

To install the package in R:
``` 
library(devtools)  
devtools::install_github("ddegras/matchFeat/matchFeat") 
```

## Reference
Degras (2021). Scalable Feature Matching Across Large Data Collections. https://arxiv.org/abs/2101.02035
