# Fused low-rank plus sparse decomposition
This toolbox contains the Matlab implementation of fused low-rank plus sparse (fusedLS) decomposition for multi-subject functional connectivity (FC) data. The decomposition is accomplished by solving a fused version of principal component pursuit (PCP) with an additional fusion-type penalty on the differences between the columns of the low-rank matrix. The optimization is solved via a linearized alternating direction method of multipliers (ADMM). The method is useful for separating shared/correlated and subject-specific FC patterns from observed multi-subject FCs, represented by the low-rank and sparse components, respectively. For more details, refer to [Ting, et. al (2021)](https://arxiv.org/abs/2102.10331#) with an application to separating stimulus-induced and background components of dynamic functional connectivity in naturalistic fMRI.

## Contents
The toolbox currently includes the following:
- Implementation of linearized ADMM algorithm for fused low-rank plus sparse decomposition (main function `fused_ls`)
- Demo of fusedLS decomposition applied to simulated multi-subject connectivity data (script `example_simulation.m`)

## Usage
The function `fused_ls()` implements the linearized ADMM algorithm for fused L+S decomposition (see Section III and IV in [Ting, et. al (2021)](https://arxiv.org/abs/2102.10331#)) of an input data matrix `Z`

`[L, S, rL, err] = fused_ls(Z, lambda1, lambda2, mu, tol, max_iter)`

where the outputs `L` and `S` are respectively the recorved low-rank and sparse components, `rL` is the estimated rank and `err` is the reconstruction error. See documentation in the function header for more details on usage and information on input parameters.

## Demo
A demo of the fusedLS algorithm for separating shared/correlated structure and subject-specific backgroud components in multi-subject FC networks is provided in the script `example_simulation.m`. The multi-subject simulated networks are generated as a sum of correlated FC patterns across subjects and sparse individual background components. See Section IV in [Ting, et. al (2021)](https://arxiv.org/abs/2102.10331#) for details of simulation settings. External function `generateSbm()` from [Dynamic Stochastic Block Models MATLAB Toolbox](https://github.com/IdeasLabUT/Dynamic-Stochastic-Block-Model) is required to simulate the community structure of the underlying correlated FC patterns.

## Citation
Chee-Ming Ting, Jeremy I Skipper, Steven L Small, Hernando Ombao "Separating Stimulus-Induced and Background Components of Dynamic Functional Connectivity in Naturalistic fMRI"
[arXiv preprint arXiv:2102.10331](https://arxiv.org/abs/2102.10331#).
