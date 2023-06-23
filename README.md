
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsExtremes

<!-- badges: start -->
<!-- badges: end -->

The package tsExtremes gathers functions for analyzing the extremes of
stationary time series and assessing the time dependencies of extremes.

For example, you can use it compute estimates of:

- the tail index of the series,

- the extremal index of a time series,

- the cluster lengths of the time series.

# Installation

<!-- You can install the released version of tsExtremes from [CRAN](https://CRAN.R-project.org) with: -->

To get started, install the package from Github using the command:

``` r
devtools::install_github('GBuritica/tsExtremes')
```

# R Package tutorial

The outline of this tutorial is as follows:

1.  Tail index inference
    - In this vignette you can learn how to compute unbiased estimates
      of the tail index alpha of the series. The implementation is based
      on the Hill-type estimators in (Haan, Mercadier, and Zhou 2016).
2.  Cluster inference
    - In this window you will learn how to compute the cluster process
      based estimates in (Buriticá et al. 2021) for classical cluster
      statics. We consider the example of inference of the extremal
      index with the alpha cluster process (Buriticá, Mikosch, and
      Wintenberger 2022), and also an example of the cluster lengths.

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
of Lp Blocks of Regularly Varyinig Time Series.”

</div>

<div id="ref-dehaan:mercadier:zhou:2016" class="csl-entry">

Haan, L. de, G. Mercadier, and C. Zhou. 2016. “Adapting Extreme Value
Statistics to Financial Time Series: Dealing Wih Bias and Serial
Dependence.” *Finance and Stochastics* 20: 321–54.

</div>

</div>
