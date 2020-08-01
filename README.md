# matchFeats
One-To-One Feature Matching

Statistical methods to match feature vectors between sample units under association ambiguity. Applications include object tracking, video surveillance, remote sensing as well as multilevel modeling. Given a fixed number of classes/distributions, for each unit, exactly one vector of each class is observed without label. The goal is to label the feature vectors using each label exactly once so to produce the best match across units, e.g. by minimizing the variability within classes. Several statistical solutions based on empirical loss functions and probabilistic modeling are provided.  

To install the package in R:
``` 
install.packages("clue") # run this line if you don’t have package 'clue' installed yet
install.packages("foreach") # run this line if you don’t have package 'foreach' installed yet
install.packages("devtools") # run this line if you don’t have package 'devtools' installed yet
library(devtools)  
install_github("ddegras/matchFeats") 
```

