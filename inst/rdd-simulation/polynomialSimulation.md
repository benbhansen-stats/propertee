---
title: "Polynomial Simulation like in LRD paper"
author: "Adam C Sales & Ben B Hansen"
date: "05 August, 2022"
output: html_document
---



General dependencies.

```r
library('knitr')
library('kableExtra')
#library(flexida)
source('simulationFunctions.r')
source('ddsandwich.R') ###
source('displaySim.r')
```


```r
devtools::load_all("../..")
```

Initialization. Note that `nreps=0` corresponds to no simulations,
just print results from previously saved simulations.
In order to re-run the simulations, the `nreps`
variable should have been set to a positive integer before initiating this script.

To run the simulations in parallel, using the `parallel` package in
`R`,
register a cluster, called `cl` with the desired number of nodes, with
code similar to the following:

```r
library(parallel)
cl <- makeCluster(5)
```


```r
if (!exists('nreps') ) nreps <- 0
nreps
```

```
## [1] 0
```

```r
if (nreps) {
library('robustbase')
library('rdd')
library('RItools')
library('sandwich')
library('nnet')

clust <- FALSE
if(require('parallel')) if(exists('cl')) if(inherits(cl,"cluster")) clust <- TRUE

if(clust){
  clusterEvalQ(cl,{

    library('robustbase')
    library('rdd')
    library('RItools')
    library('sandwich')
    library('nnet')
    devtools::load_all("../..")
                                        #           library(flexida)
    source('simulationFunctions.r')
    source('ddsandwich.R') ###
    source('displaySim.r')
  })
} else cl <- NULL
}
```
## Run the simulation


```r
if (nreps) {
  set.seed(201609)
  st2 <- system.time(totalPoly <- totalPolySim(nreps,cl))
  save(totalPoly,file=paste0("./totalPolySim",Sys.Date(),".RData"))
  cat(paste0(date(), ', nreps=', nreps, '\n'),
      paste(c(names(st2),'\n', collapse=T)),
      st2,
      file='totalPolySim-runtime.txt', append=TRUE)
} else{
  psims <- sort(grep('totalPolySim',list.files('.'),value=TRUE),decreasing=TRUE)

  load(paste0('./',psims[1]))
}
```

The following gives the results in Table 4 of the paper, in addition
to the break-down of RMSE into bias and variance, and analogous
results for normally-distributed errors.


```r
tab <- prntTab(totalPoly,5,full=TRUE,md=FALSE)
#rownames(tab) <- rep(c('level','RMSE','bias','sd'),6)
colnames(tab) <- gsub('(flex|ols)\\.se.','deg=',colnames(tab))#c(rep(paste0('deg=',1:4),2),'')
colnames(tab)[ncol(tab)] <- 'n/a'
kable(tab,format='html',caption='Full results for polynomial simulation',digits=2)%>%
    kable_styling()%>% column_spec( 6,border_right=TRUE)%>%column_spec(11,border_right=TRUE)%>%
        add_header_above(c(" " = 1, "flexida" = 5, "OLS" = 5, "Loc. Lin." = 1))%>%
            #group_rows("$t_3$ Error",1,12)%>%group_rows("N(0,1) Error",13,24)%>%
            group_rows("linear",1,5)%>%group_rows('antiSym',6,10)%>%group_rows('sine',11,15)#%>%
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Full results for polynomial simulation</caption>
 <thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1"></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">flexida</div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">OLS</div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="1"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">Loc. Lin.</div></th>
</tr>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> deg=1 </th>
   <th style="text-align:right;"> deg=2 </th>
   <th style="text-align:right;"> deg=3 </th>
   <th style="text-align:right;"> deg=4 </th>
   <th style="text-align:right;"> deg=5 </th>
   <th style="text-align:right;"> deg=1 </th>
   <th style="text-align:right;"> deg=2 </th>
   <th style="text-align:right;"> deg=3 </th>
   <th style="text-align:right;"> deg=4 </th>
   <th style="text-align:right;"> deg=5 </th>
   <th style="text-align:right;"> n/a </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="5"><td colspan="12" style="border-bottom: 1px solid;"><strong>linear</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> VarRelBias </td>
   <td style="text-align:right;"> -0.61 </td>
   <td style="text-align:right;"> -0.61 </td>
   <td style="text-align:right;"> -0.76 </td>
   <td style="text-align:right;"> -0.77 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.83 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.01 </td>
   <td style="text-align:right;"> 0.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> RMSE </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.38 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.96 </td>
   <td style="text-align:right;"> 0.49 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> bias </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.01 </td>
   <td style="text-align:right;"> -0.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> sd </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.38 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.96 </td>
   <td style="text-align:right;"> 0.49 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> ciCover </td>
   <td style="text-align:right;"> 0.77 </td>
   <td style="text-align:right;"> 0.76 </td>
   <td style="text-align:right;"> 0.65 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.57 </td>
   <td style="text-align:right;"> 0.95 </td>
   <td style="text-align:right;"> 0.95 </td>
   <td style="text-align:right;"> 0.95 </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.94 </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr grouplength="5"><td colspan="12" style="border-bottom: 1px solid;"><strong>antiSym</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> VarRelBias </td>
   <td style="text-align:right;"> -0.64 </td>
   <td style="text-align:right;"> -0.64 </td>
   <td style="text-align:right;"> -0.78 </td>
   <td style="text-align:right;"> -0.78 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.84 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> -0.03 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> -0.04 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.06 </td>
   <td style="text-align:right;"> 0.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> RMSE </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.41 </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.65 </td>
   <td style="text-align:right;"> 0.81 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.99 </td>
   <td style="text-align:right;"> 0.51 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> bias </td>
   <td style="text-align:right;"> -0.63 </td>
   <td style="text-align:right;"> -0.63 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.13 </td>
   <td style="text-align:right;"> -0.62 </td>
   <td style="text-align:right;"> 0.16 </td>
   <td style="text-align:right;"> 0.17 </td>
   <td style="text-align:right;"> -0.06 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.06 </td>
   <td style="text-align:right;"> 0.02 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> sd </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.38 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.81 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.99 </td>
   <td style="text-align:right;"> 0.51 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> ciCover </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.11 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.53 </td>
   <td style="text-align:right;"> 0.46 </td>
   <td style="text-align:right;"> 0.93 </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.94 </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr grouplength="5"><td colspan="12" style="border-bottom: 1px solid;"><strong>sine</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> VarRelBias </td>
   <td style="text-align:right;"> -0.63 </td>
   <td style="text-align:right;"> -0.63 </td>
   <td style="text-align:right;"> -0.76 </td>
   <td style="text-align:right;"> -0.76 </td>
   <td style="text-align:right;border-right:1px solid;"> -0.83 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.00 </td>
   <td style="text-align:right;"> 0.08 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> RMSE </td>
   <td style="text-align:right;"> 1.19 </td>
   <td style="text-align:right;"> 1.19 </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.38 </td>
   <td style="text-align:right;"> 1.21 </td>
   <td style="text-align:right;"> 0.48 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.97 </td>
   <td style="text-align:right;"> 0.53 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> bias </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> 1.16 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.02 </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> -0.10 </td>
   <td style="text-align:right;"> -0.06 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.01 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> sd </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.38 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.97 </td>
   <td style="text-align:right;"> 0.53 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> ciCover </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.56 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;"> 0.95 </td>
   <td style="text-align:right;"> 0.95 </td>
   <td style="text-align:right;border-right:1px solid;"> 0.94 </td>
   <td style="text-align:right;"> 0.67 </td>
  </tr>
</tbody>
</table>

```r
                #group_rows("linear",13,16)%>%group_rows('antiSym',17,20)%>%group_rows('oneSide',21,24)
```


Session information

```r
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19044)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] flexida_0.0.0.9000 nnet_7.3-16        RItools_0.1-18     SparseM_1.81      
##  [5] rdd_0.57           Formula_1.2-4      AER_1.2-9          survival_3.2-11   
##  [9] car_3.0-11         carData_3.0-4      lmtest_0.9-39      zoo_1.8-9         
## [13] sandwich_3.0-1     robustbase_0.93-9  testthat_3.1.0     kableExtra_1.3.4  
## [17] knitr_1.36        
## 
## loaded via a namespace (and not attached):
##  [1] svd_0.5           httr_1.4.2        pkgload_1.2.3     viridisLite_0.4.0
##  [5] splines_4.1.1     highr_0.9         cellranger_1.1.0  remotes_2.4.1    
##  [9] sessioninfo_1.2.1 pillar_1.6.4      lattice_0.20-44   glue_1.4.2       
## [13] digest_0.6.28     rvest_1.0.1       colorspace_2.0-2  Matrix_1.3-4     
## [17] htmltools_0.5.2   pkgconfig_2.0.3   devtools_2.4.2    haven_2.4.3      
## [21] xtable_1.8-4      purrr_0.3.4       scales_1.1.1      webshot_0.5.2    
## [25] processx_3.5.2    svglite_2.0.0     openxlsx_4.2.4    rio_0.5.27       
## [29] tibble_3.1.5      usethis_2.1.5     ellipsis_0.3.2    cachem_1.0.6     
## [33] withr_2.4.2       cli_3.1.0         magrittr_2.0.1    crayon_1.4.2     
## [37] readxl_1.3.1      memoise_2.0.0     evaluate_0.14     ps_1.6.0         
## [41] fs_1.5.0          fansi_0.5.0       forcats_0.5.1     xml2_1.3.2       
## [45] foreign_0.8-81    pkgbuild_1.2.0    tools_4.1.1       data.table_1.14.2
## [49] prettyunits_1.1.1 hms_1.1.1         lifecycle_1.0.1   stringr_1.4.0    
## [53] munsell_0.5.0     zip_2.2.0         callr_3.7.0       compiler_4.1.1   
## [57] systemfonts_1.0.3 rlang_0.4.11      grid_4.1.1        rstudioapi_0.13  
## [61] rmarkdown_2.11    abind_1.4-5       curl_4.3.2        R6_2.5.1         
## [65] fastmap_1.1.0     utf8_1.2.2        rprojroot_2.0.2   desc_1.4.0       
## [69] stringi_1.7.5     Rcpp_1.0.7        vctrs_0.3.8       DEoptimR_1.0-9   
## [73] xfun_0.26
```
