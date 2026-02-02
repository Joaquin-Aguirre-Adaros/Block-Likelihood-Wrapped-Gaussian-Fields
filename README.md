# Block-Likelihood-Wrapped-Gaussian-Fields

Code repository for thesis work on block likelihood in the context of Wrapped-Gaussian fields

## Overview

This repository is divided into two principal folders:

1. **simulation_study**: Provides the code and files obtained in the simulation study
2. **application**: Contains the code used for the applicability of this method on a real dataset

## Purpose

The purpose of this work is to provide an alternative for estimating Wrapped-Gaussian fields without requiring Monte Carlo simulations, as needed in:
- Bayesian frameworks
- EM algorithm implementations

We propose **"Block Likelihood"**, a method based on composite likelihood that generalizes the so-called **"pairwise likelihood"**. Pairwise likelihood is considered our benchmark and is currently the only frequentist method that has been tested in previous work (Alegría et al., 2024). Therefore, we compare the performance of both methods throughout this repository.

## Reference

Alegría, A., Caro, S., Bevilacqua, M., Porcu, E., & Clarke, J. (2024). Estimating space and space-time covariance functions for large data sets: a weighted composite likelihood approach. *Journal of the American Statistical Association*, 119(545), 350-360. https://doi.org/10.1080/01621459.2022.2119985

