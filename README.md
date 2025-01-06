# Loewner_vs_Hankel
This MATLAB code compares Loewner Rank Minimization and Hankel Rank Minimization for system identification of LTI systems. It utilizes YALMIP for optimization and was used to conduct the experiments presented in the paper "Identification of Low-Order Systems in a Loewner Framework".

The main file is "Time_Domain_Interpolation.m". It generates a random stable LTI system along with corresponding training data with noise. Then, it parameterizes all LTI systems that are consistent with prior assumptions and could potentially generate the data using Linear Fractional Transformation (LFT). Among these parameterized systems, it identifies a high-order stable system using one of the following two methods:

 - Minimizes the rank of the Hankel matrix formed from the impulse response. This is implemented in "HankelRankMinimization.m".
 - Minimizes the rank of the Loewner matrix formed from the frequency response. This is implemented in "LoewnerRankMinimization.m".

After that, Model Order Reduction (MOR) is applied to the Hankel or Loewner matrices. This is implemented in:

 - "Hankel_Reduction.m"
 - "Loewner_Reduction.m"

Both reduction methods are similar to the Ho-Kalman algorithm. At the end, the results are compared.

If you use this code, please cite the following paper:

**Identification of Low Order Systems in a Loewner Framework**  
Arya Honarpisheh, Rajiv Singh, Jared Miller, Mario Sznaier
[DOI or Conference Link](https://doi.org/10.1016/j.ifacol.2024.08.528)  

BibTeX:
```bibtex
@article{HONARPISHEH2024199,
title = {Identification of Low Order Systems in a Loewner Framework},
journal = {IFAC-PapersOnLine},
volume = {58},
number = {15},
pages = {199-204},
year = {2024},
note = {20th IFAC Symposium on System Identification SYSID 2024},
issn = {2405-8963},
doi = {https://doi.org/10.1016/j.ifacol.2024.08.528},
url = {https://www.sciencedirect.com/science/article/pii/S2405896324013089},
author = {Arya Honarpisheh, Rajiv Singh, Jared Miller, and Mario Sznaier},
}
