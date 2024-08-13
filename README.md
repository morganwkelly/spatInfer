# spatInfer

The purpose of `spatInfer` is to run spatial regressions that are robust to trends and autocorrelation in the data, and to test the accuracy of the inference method with a spatial noise placebo test. At the same it aims to be extremely simple to use, requiring only a sequence of four commands.

## First Steps

The goal is to run a regression with a spatial basis made up of the first $p$ principal components of a $k \times k$ tensor and with standard errors based on $c$ large clusters. The spatial basis controls for trends and other large scale structure in the variables, and the clusters deal with autocorrelation in residuals.



The parameters $c$, $k$ and $p$ are esimated with two commands.

The first is `optimal_basis`:

``` r
library(spatInfer)

optimal_basis(y ~ x1 + x2, dat)
```

This chooses the combination of $k$ and $p$ that optimally predicts the outcome $y$ based on a Bayes Information Criterion.

Knowing $k$ and $p$, we next run the `placebo` command to choose the optimal number $c$ of large clusters by comparing the p-values of spatial noise placebos against that of the treatment $x1$.

``` r
placebo(y ~ x1 + x2, dat,
            splines=k,
            num_pcs=p
        )
```

For each number of clusters $c$, `placebo` gives the proportion of simulations with p-values below 0.05. The optimal cluster number is where this proportion is in the range 0.05 to 0.07. The command also gives the placebo significance level of the treatment variable: how often the t-statistics of spatial noise placebos are above the real one.

Knowing the spatial basis values $k$ and $p$ and the cluster number $c$ we can now run the spatial basis regression `basis_reg` that gives us standard regression coefficients, confidence intervals and so on.

``` r
basis_reg(y ~ x1 + x2, dat,
           clusters=c,
           splines=k,
           num_pcs=p
          )
```

Besides `placebo`, `spatInfer` gives a second diagnostic `synth`, a synthetic outcome test that tests the hypothesis that the outcome $y$ is spatial noise, and thus independent of the treatment $x1$.

``` r
synth(y ~ x1 + x2, dat,
            splines=k,
            num_pcs=p
)
```
