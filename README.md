

---

# sensiverse

[![R-CMD-check](https://github.com/mganslmeier/sensiverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mganslmeier/sensiverse/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/sensiverse)](https://CRAN.R-project.org/package=sensiverse)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## What is sensiverse?

**sensiverse** is an R package to systematically explore *model uncertainty* in quantitative social science by replicating the approach in:

> *Estimating the extent and sources of model uncertainty in political science*  
> Michael Ganslmeier & Tim Vlandas (2025). *Proceedings of the National Academy of Sciences*. DOI: [10.1073/pnas.2414926122](https://www.pnas.org/doi/10.1073/pnas.2414926122)

In that paper, the authors combine Extreme Bounds Analysis with a multiverse perspective to evaluate how empirical conclusions shift when multiple defensible modeling choices are varied in tandem:

- **Dimensions of variation**: choice of control covariates, fixed‐effect structure, standard‐error type, sample subset, dependent variable operationalization  
- They apply this across four political science topics (e.g. welfare generosity, democratization, trust, public goods) to estimate **over 3.6 billion regression coefficients**
- Their empirical conclusion: model uncertainty is substantial, not only in statistical significance but also in effect direction. Among the sources of sensitivity, **sample definition** and **dependent variable specification** tend to dominate over choice of covariates  

The **sensiverse** package operationalizes this framework, allowing a user to:

1. Build a specification grid from a replication‐style dataset defining various model choices
2. Estimate the full model space
3. Aggregate results into *significance shares*  
4. Optionally apply meta‐models (e.g. multinomial logit, nearest neighbor, neural networks) to explain which specification choices drive fragility  
5. Visualize and interpret robustness / fragility similar to the paper’s figures  

---

## Installation

```r
# From CRAN (once published)
install.packages("sensiverse")

# Development version
# install.packages("remotes")
remotes::install_github("mganslmeier/sensiverse")
````

---

## Quick Start Example

Here’s a minimal workflow illustrating core functionality: [https://mganslmeier.github.io/sensiverse-tutorial/](https://mganslmeier.github.io/sensiverse-tutorial/)

---

## Citation

If you use **sensiverse**, please cite:

```bibtex
@article{GanslmeierVlandas2025,
  title   = {Estimating the extent and sources of model uncertainty in political science},
  author  = {Michael Ganslmeier and Tim Vlandas},
  journal = {Proceedings of the National Academy of Sciences},
  year    = {2025},
  doi     = {10.1073/pnas.2414926122}
}
```

---

## License & Authors

This package is licensed under **MIT License** (see `LICENSE`).
Authored by **Michael Ganslmeier** and **Tim Vlandas**.

---

## Acknowledgments & Resources

* Replication materials / code archive: [https://doi.org/10.5281/zenodo.15480536](https://doi.org/10.5281/zenodo.15480536)
* Many thanks to users and reviewers who will help test on a broad set of empirical data

---

## Contact & Contribution

Feel free to open issues or submit pull requests on GitHub. For substantive questions, you can also email **Michael Ganslmeier** (m.ganslmeier@exeter.ac.uk) or **Tim Vlandas** (tim.vlandas@spi.ox.ac.uk).
