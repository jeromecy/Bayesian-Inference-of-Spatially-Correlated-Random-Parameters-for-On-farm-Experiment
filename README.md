# Bayesian Inference of Spatially Correlated Random Parameters for On-Farm Experiment

## Overview

This repository contains the code and data for the paper titled "Bayesian inference of spatially correlated random parameters for on-farm experiment" published in Field Crops Research. The paper can be accessed at https://doi.org/10.1016/j.fcr.2022.108477.

## Abstract

**🌾 Context or Problem**  
Accounting for spatial variability is crucial when estimating treatment effects in large on-farm trials. It helps determine the optimal treatment for every part of a paddock, leading to a management strategy that enhances sustainability and profitability.

**🎯 Objective**  
We specify a model with spatially correlated random parameters to account for spatial variability in large on-farm trials. A Bayesian framework is adopted to estimate the posterior distribution of these parameters, allowing for the estimation of spatially-varying treatment effects.

**🔬 Methods**  
Several approaches have been proposed for assessing spatial variability, but they often overlook the issue of model misspecification. We use Bayesian post-sampling tools to diagnose model misspecification, showing that the Gaussian distribution assumption for the response variable is often inadequate. Instead, the Student-t distribution provides more robust inference.

**📊 Results**  
We applied our method to a real on-farm strip trial from Las Rosas, Argentina, aiming to obtain a spatial map of locally-varying optimal nitrogen rates. The analysis revealed that the Student-t distribution is more appropriate than the Gaussian distribution for the response variable.

**🏆 Conclusions**  
Our Bayesian approach is compared with geographically weighted regression, highlighting the differences and advantages of each method.

**🌟 Implications**  
This research provides practical recommendations for designing and analyzing large on-farm trials, contributing valuable insights for precision agriculture and site-specific management.

## Repository Contents

- **Code**: Scripts used for Bayesian inference and analysis.
- **Data**: Datasets used in the study.
- **Manuscript**: The full manuscript of the published paper.

## Citation

If you use this repository in your research, please cite the paper as follows:

```bibtex
@article{cao2024optimal,
  title={Bayesian inference of spatially correlated random parameters for on-farm experiment},
  author={Cao, Zhanglong and Stefanova, Katia and Gibberd, Mark and Rakshit, Suman},
  journal={Field Crops Research},
  volume={318},
  pages={108477},
  year={2022},
  doi={https://doi.org/10.1016/j.fcr.2022.108477},
  publisher={Elsevier}
}
```
