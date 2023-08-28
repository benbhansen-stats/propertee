Polynomial Simulation like in LRD paper
================
Adam C Sales & Ben B Hansen & Josh Wasserman
August 2022

General dependencies.

``` r
library('knitr')
library('kableExtra')
library('testthat')
source('displaySim.r')
devtools::load_all("../..")
```

The following gives the results in Table 4 of the paper, in addition to
the break-down of RMSE into bias and variance, and analogous results for
normally-distributed errors.

``` r
load("totalPolySim.Rdata")
tab <- prntTab(totalPoly,5,full=TRUE)
#rownames(tab) <- rep(c('level','RMSE','bias','sd'),6)
colnames(tab) <- gsub('(flex|ols)\\.se.','deg=',colnames(tab))#c(rep(paste0('deg=',1:4),2),'')
# colnames(tab)[ncol(tab)] <- 'n/a' # this is for LLR column
kable(tab,format='html',caption='Full results for polynomial simulation',digits=2) %>%
  kable_styling() %>%
  column_spec(6,border_right=TRUE) %>%
  column_spec(11,border_right=TRUE) %>%
  # add_header_above(c(" " = 1, "propertee" = 5, "OLS" = 5, "Loc. Lin." = 1))%>%
  add_header_above(c(" " = 1, "propertee" = 5, "OLS" = 5))%>%
  #group_rows("$t_3$ Error",1,12)%>%group_rows("N(0,1) Error",13,24)%>%
  group_rows("linear",1,5) %>%
  group_rows('antiSym',6,10) %>%
  group_rows('sine',11,15) #%>% group_rows("linear",13,16)%>%group_rows('antiSym',17,20)%>%group_rows('oneSide',21,24)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Full results for polynomial simulation
</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

propertee

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

OLS

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
deg=1
</th>
<th style="text-align:right;">
deg=2
</th>
<th style="text-align:right;">
deg=3
</th>
<th style="text-align:right;">
deg=4
</th>
<th style="text-align:right;">
deg=5
</th>
<th style="text-align:right;">
deg=1
</th>
<th style="text-align:right;">
deg=2
</th>
<th style="text-align:right;">
deg=3
</th>
<th style="text-align:right;">
deg=4
</th>
<th style="text-align:right;">
deg=5
</th>
</tr>
</thead>
<tbody>
<tr grouplength="5">
<td colspan="11" style="border-bottom: 1px solid;">
<strong>linear</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
VarRelBias
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;border-right:1px solid;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;border-right:1px solid;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
RMSE
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;border-right:1px solid;">
0.38
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;border-right:1px solid;">
1.03
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
bias
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;border-right:1px solid;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;border-right:1px solid;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
sd
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;border-right:1px solid;">
0.38
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;border-right:1px solid;">
1.03
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
ciCover
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;border-right:1px solid;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;border-right:1px solid;">
0.94
</td>
</tr>
<tr grouplength="5">
<td colspan="11" style="border-bottom: 1px solid;">
<strong>antiSym</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
VarRelBias
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;border-right:1px solid;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;border-right:1px solid;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
RMSE
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;border-right:1px solid;">
0.40
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;border-right:1px solid;">
0.94
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
bias
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;border-right:1px solid;">
0.13
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;border-right:1px solid;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
sd
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;border-right:1px solid;">
0.38
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;border-right:1px solid;">
0.94
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
ciCover
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;border-right:1px solid;">
0.93
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;border-right:1px solid;">
0.95
</td>
</tr>
<tr grouplength="5">
<td colspan="11" style="border-bottom: 1px solid;">
<strong>sine</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
VarRelBias
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;border-right:1px solid;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;border-right:1px solid;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
RMSE
</td>
<td style="text-align:right;">
1.18
</td>
<td style="text-align:right;">
1.18
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;border-right:1px solid;">
0.37
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;border-right:1px solid;">
0.95
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
bias
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;border-right:1px solid;">
0.01
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;border-right:1px solid;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
sd
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;border-right:1px solid;">
0.37
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;border-right:1px solid;">
0.95
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
ciCover
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;border-right:1px solid;">
0.95
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;border-right:1px solid;">
0.94
</td>
</tr>
</tbody>
</table>

``` r
N_PAPER_SIMS <- 5000
MIN_PAPER_COVERAGE <- 0.93
N_TEST_SIMS <- ncol(totalPoly[[1]])
MIN_COVERAGE_BOUND <- MIN_PAPER_COVERAGE - 2 * sqrt(0.95 * 0.05 / N_PAPER_SIMS) -
  2 * sqrt(0.95 * 0.05 / N_TEST_SIMS)
MIN_COVERAGE_BOUND
```

    ## [1] 0.9176712

``` r
testthat::expect_true(all(rownames(tab)[c(5, 10, 15)] == "ciCover"))
testthat::expect_true(
  all(tab[5, 1:5] > MIN_COVERAGE_BOUND)
)
testthat::expect_true(
  all(tab[10, 3:5] > MIN_COVERAGE_BOUND)
)
testthat::expect_true(
  all(tab[15, 3:5] > MIN_COVERAGE_BOUND)
)
```

Session information

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] propertee_0.0.0.9000 testthat_3.1.4     kableExtra_1.3.4   knitr_1.39        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7        svglite_2.1.0     lattice_0.20-45   prettyunits_1.1.1
    ##  [5] ps_1.6.0          zoo_1.8-11        rprojroot_2.0.2   digest_0.6.28    
    ##  [9] mime_0.12         R6_2.5.1          evaluate_0.15     highr_0.9        
    ## [13] httr_1.4.2        rlang_1.0.5       rstudioapi_0.13   miniUI_0.1.1.1   
    ## [17] callr_3.7.0       urlchecker_1.0.1  rmarkdown_2.17    desc_1.4.1       
    ## [21] devtools_2.4.4    webshot_0.5.2     stringr_1.4.0     htmlwidgets_1.5.4
    ## [25] munsell_0.5.0     shiny_1.7.1       compiler_4.1.2    httpuv_1.6.3     
    ## [29] xfun_0.33         systemfonts_1.0.4 pkgbuild_1.3.1    htmltools_0.5.2  
    ## [33] viridisLite_0.4.0 crayon_1.4.2      withr_2.5.0       later_1.3.0      
    ## [37] waldo_0.4.0       brio_1.1.2        grid_4.1.2        xtable_1.8-4     
    ## [41] lifecycle_1.0.1   magrittr_2.0.1    scales_1.1.1      cli_3.3.0        
    ## [45] stringi_1.7.5     cachem_1.0.6      fs_1.5.2          promises_1.2.0.1 
    ## [49] remotes_2.4.2     xml2_1.3.2        ellipsis_0.3.2    sandwich_3.0-2   
    ## [53] tools_4.1.2       glue_1.6.2        purrr_0.3.4       processx_3.5.2   
    ## [57] pkgload_1.3.0     fastmap_1.1.0     yaml_2.2.1        colorspace_2.0-2 
    ## [61] sessioninfo_1.2.2 rvest_1.0.1       memoise_2.0.1     profvis_0.3.7    
    ## [65] usethis_2.1.6
