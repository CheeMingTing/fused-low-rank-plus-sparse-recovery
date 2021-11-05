# Fused low-rank plus sparse decomposition
This toolbox contains the Matlab implementation of fused low-rank plus sparse (fusedLS) decomposition for multi-subject functional connectivity (FC) data. The decomposition is accomplished by solving a fused version of principal component pursuit (PCP) with an additional fusion-type penalty on the differences between the columns of the low-rank matrix. The optimization is solved via a linearized alternating direction method of multipliers (ADMM). The method is useful for separating shared/correlated and subject-specific FC patterns from observed multi-subject FCs, represented by the low-rank and sparse components, respectively. For more details, refer to [Ting, et. al (2021)](https://arxiv.org/abs/2102.10331#) with an application to separating stimulus-induced and background components of dynamic functional connectivity in naturalistic fMRI.

## Contents
The toolbox currently includes the following:
- Implementation of linearized ADMM algorithm for fused low-rank plus sparse decomposition (main function `fused_ls`)
- Demo of fusedLS decomposition applied to simulated multi-subject connectivity data (script `example_simulation.m`)

## Reference
Chee-Ming Ting, Jeremy I Skipper, Steven L Small, Hernando Ombao "Separating Stimulus-Induced and Background Components of Dynamic Functional Connectivity in Naturalistic fMRI"
[arXiv preprint arXiv:2102.10331](https://arxiv.org/abs/2102.10331#).
