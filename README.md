
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsExtremes package tutorial

<!-- badges: start -->
<!-- badges: end -->

The package tsExtremes includes functions for studying extreme features
of heavy-tailed time series. It can be used to compute estimates of the
tail index, the extremal index or cluster lengths of the series. These
are statistics of the series summarizing the magnitude of extreme
records, and the number of concomitant records.

## Installation

<!-- You can install the released version of tsExtremes from [CRAN](https://CRAN.R-project.org) with: -->

To get started, install the package from Github using the command:

``` r
devtools::install_github('GBuritica/tsExtremes')
```

## R Package tutorial

This package tutorial includes two main windows on:

1.  Tail index inference
    - In this vignette you will learn how to compute estimates of the
      tail index $\alpha$ of the series. The implementation is based on
      the Hill-type estimators in (Haan, Mercadier, and Zhou 2016).
2.  Cluster inference
    - In this vignette you will learn how to implement block estimates
      of cluster statistics; see (Buriticá, Mikosch, and Wintenberger
      2022). We consider the example of the extremal index based on the
      estimator proposed in (Buriticá et al. 2021), and also an example
      of cluster lengths (Buriticá and Wintenberger 2024).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-buritica:meyer:mikosch:wintenberger:2021"
class="csl-entry">

Buriticá, G., N. Meyer, T. Mikosch, and O. Wintenberger. 2021. “Some
Variations on the Extremal Index.” *Zap. Nauchn. Semin. POMI.* 30:
52–77.

</div>

<div id="ref-buritica:mikosch:wintenberger:2021" class="csl-entry">

Buriticá, G., T. Mikosch, and O. Wintenberger. 2022. “Large Deviations
of Lp Blocks of Regularly Varyinig Time Series.” *Stochastic Processes
and Their Applications.* 161: 68–101.

</div>

<div id="ref-buritica:wintenberger:2024" class="csl-entry">

Buriticá, G., and O. Wintenberger. 2024. “On the Asymptotics of Extremal
Lp-Blocks Cluster Inference.” *arXiv:2212.13521*.

</div>

<div id="ref-dehaan:mercadier:zhou:2016" class="csl-entry">

Haan, L. de, G. Mercadier, and C. Zhou. 2016. “Adapting Extreme Value
Statistics to Financial Time Series: Dealing Wih Bias and Serial
Dependence.” *Finance and Stochastics* 20: 321–54.

</div>

</div>
