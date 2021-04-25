# STCMTL

R package for STCMTL. 

Motivated by the idea of semisoft clustering of data, we propose a semisoft task clustering approach STCMTL, which can simultaneously reveal the task cluster structure for both pure and mixed tasks as well as select the relevant features. The main assumption behind our approach is that each cluster has some pure tasks, and each mixed task can be represented by a linear combination of pure tasks in different clusters. To solve the resulting non-convex constrained optimization problem, we design an efficient three-step algorithm. The experimental results based on synthetic and real-world datasets validate the effectiveness and efficiency of the proposed approach. Finally, we extend the proposed approach to a robust task-clustering problem.

# Installation 
This package can be installed in R:
```
library("devtools")
install.packages("STCMTL_0.0.0.9000.tar.gz",type = "source",repos = NULL)

```

# Usage

For how to use the package, please download and read Usage.html.
