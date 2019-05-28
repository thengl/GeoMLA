Mapping Prediction Uncertainty using RFsp
================
Hengl, T. and Wright, M.N.

| <a href="https://github.com/thengl"><img src="https://avatars0.githubusercontent.com/u/640722?s=460&v=4" height="100" alt="Tomislav Hengl"></a> | <a href="https://github.com/mnwright"><img src="https://avatars3.githubusercontent.com/u/9598192?s=460&v=4" height="100" alt="Marvin N. Wright"></a> |
| ----------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |

-----

<a href="https://creativecommons.org/licenses/by-sa/4.0/" target="_blank"><img src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" alt=""></a>

-----

There is an increasing interest in using Machine Learning techniques for
the purpose of generating spatial predictions, and for mining
environmental data in general. Machine Learning algorithms, such as
random forests, artificial neural networks and support vector machines
have already shown predictive potential for various environmental
applications (Biau & Scornet, 2016; Nussbaum et al., 2018; Prasad,
Iverson, & Liaw, 2006). This tutorial explains how to use Machine
Learning to efficiently generate spatial predictions and derive
associated uncertainty. Our focus is on using Random Forest as
implemented in the [ranger](https://github.com/imbs-hl/ranger) package
(Wright & Ziegler, 2017), although similar frameworks could be applied
to other tree-based Machine Learning algorithms.

For a complete overview of methods used please refer to: Hengl, T.,
Nussbaum, M., Wright, M. and Heuvelink, G.B.M., 2018. [*“Random Forest
as a Generic Framework for Predictive Modeling of Spatial and
Spatio-temporal Variables”*](https://peerj.com/preprints/26693/), PeerJ
(in review).

## Software installation

Software (required):

  - [R](http://cran.r-project.org/bin/windows/base/) or
    [MRO](https///mran.microsoft.com/download/);

  - [RStudio](http://www.rstudio.com/products/RStudio/);

  - R packages: GSIF, plotKML, ranger, caret, plyr, raster, devtools
    (see: [how to install R
    package](http://www.r-bloggers.com/installing-r-packages/));

  - [QGIS](https://www.qgis.org/en/site/forusers/download.html) to
    visualize predictions (add
    [WMTS](http://gis.sinica.edu.tw/worldmap/wmts) background layers);

R script used in this tutorial you can download from the
**[github](https://github.com/thengl/GeoMLA/tree/master/RF_vs_kriging/R)**.
As a gentle introduction to R programming languange and spatial classes
in R we recommend following [the Geocomputation with R
book](https://geocompr.robinlovelace.net/). To open any spatial layers
in QGIS you will need to use rgdal function `writeGDAL` for rasters or
`writeOGR` for vector layers.

Quantile regression random forest and derivation of standard errors
using Jackknifing is available from ranger version \>0.9.4. To run this
tutorial it is recommended to install
[ranger](https://github.com/imbs-hl/ranger) (Wright & Ziegler, 2017)
directly from github:

``` r
devtools::install_github("imbs-hl/ranger")
```

Other packages that we use here include:

``` r
library(GSIF)
```

    ## GSIF version 0.5-5 (2019-01-04)

    ## URL: http://gsif.r-forge.r-project.org/

``` r
library(rgdal)
```

    ## Loading required package: sp

    ## rgdal: version: 1.4-3, (SVN revision 828)
    ##  Geospatial Data Abstraction Library extensions to R successfully loaded
    ##  Loaded GDAL runtime: GDAL 2.2.3, released 2017/11/20
    ##  Path to GDAL shared files: /usr/share/gdal/2.2
    ##  GDAL binary built with GEOS: TRUE 
    ##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
    ##  Path to PROJ.4 shared files: (autodetected)
    ##  Linking to sp version: 1.3-1

``` r
library(raster)
library(ranger)
quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2)
## color legend:
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", 
        "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", 
        "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
```

## Mapping uncertainty using numeric variables and the ranger package

Consider for example the meuse data set from the
[gstat](https://github.com/edzer/gstat) package:

``` r
demo(meuse, echo=FALSE)
```

We can use number of covariates to help interpolating Zinc variable. We
known from the literature that concentration of metals in soil is
controlled by river flooding and carrying upstream sediments (it is
assumed that the main source of zinc in this case is the river that
occasionally floods the area), hence we can add global surface water
occurrence (Pekel, Cottam, Gorelick, & Belward, 2016), the LiDAR-based
digital elevation model (<http://ahn.nl>) and LGN (Landelijk
Grondgebruiksbestand Nederland) land cover classes, as potential
covariates explaining zinc concentration:

``` r
dir.meuse = "../RF_vs_kriging/data/meuse/"
meuse.grid$SW_occurrence = readGDAL(paste0(dir.meuse, "Meuse_GlobalSurfaceWater_occurrence.tif"))$band1[meuse.grid@grid.index]
## flooding occurrence
meuse.grid$AHN = readGDAL(paste0(dir.meuse, "ahn.asc"))$band1[meuse.grid@grid.index]
## AHN.nl precise elevation
meuse.grid$LGN5 = as.factor(readGDAL(paste0(dir.meuse, "lgn5.asc"))$band1[meuse.grid@grid.index])
## land use class convert to indicators:
meuse.grid@data = cbind(meuse.grid@data, data.frame(model.matrix(~LGN5 - 1, 
    meuse.grid@data)))
```

In addition to the covariate layers, we also derive buffer distances so
we can account for the spatial autocorrelation (see the [RFsp
tutorial](https://github.com/thengl/GeoMLA)):

``` r
grid.dist0 <- GSIF::buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
```

which takes few seconds as it generates 155 gridded maps. Now that we
have prepared all covariates we can fit a spatial prediction model that
includes geographical and process-based
covariates:

``` r
fm1 <- as.formula(paste("zinc ~ ", paste(names(grid.dist0), collapse="+"), 
      " + SW_occurrence + dist + ", paste(paste0("LGN5", levels(meuse.grid$LGN5)), collapse = "+")))
fm1
```

    ## zinc ~ layer.1 + layer.2 + layer.3 + layer.4 + layer.5 + layer.6 + 
    ##     layer.7 + layer.8 + layer.9 + layer.10 + layer.11 + layer.12 + 
    ##     layer.13 + layer.14 + layer.15 + layer.16 + layer.17 + layer.18 + 
    ##     layer.19 + layer.20 + layer.21 + layer.22 + layer.23 + layer.24 + 
    ##     layer.25 + layer.26 + layer.27 + layer.28 + layer.29 + layer.30 + 
    ##     layer.31 + layer.32 + layer.33 + layer.34 + layer.35 + layer.36 + 
    ##     layer.37 + layer.38 + layer.39 + layer.40 + layer.41 + layer.42 + 
    ##     layer.43 + layer.44 + layer.45 + layer.46 + layer.47 + layer.48 + 
    ##     layer.49 + layer.50 + layer.51 + layer.52 + layer.53 + layer.54 + 
    ##     layer.55 + layer.56 + layer.57 + layer.58 + layer.59 + layer.60 + 
    ##     layer.61 + layer.62 + layer.63 + layer.64 + layer.65 + layer.66 + 
    ##     layer.67 + layer.68 + layer.69 + layer.70 + layer.71 + layer.72 + 
    ##     layer.73 + layer.74 + layer.75 + layer.76 + layer.77 + layer.78 + 
    ##     layer.79 + layer.80 + layer.81 + layer.82 + layer.83 + layer.84 + 
    ##     layer.85 + layer.86 + layer.87 + layer.88 + layer.89 + layer.90 + 
    ##     layer.91 + layer.92 + layer.93 + layer.94 + layer.95 + layer.96 + 
    ##     layer.97 + layer.98 + layer.99 + layer.100 + layer.101 + 
    ##     layer.102 + layer.103 + layer.104 + layer.105 + layer.106 + 
    ##     layer.107 + layer.108 + layer.109 + layer.110 + layer.111 + 
    ##     layer.112 + layer.113 + layer.114 + layer.115 + layer.116 + 
    ##     layer.117 + layer.118 + layer.119 + layer.120 + layer.121 + 
    ##     layer.122 + layer.123 + layer.124 + layer.125 + layer.126 + 
    ##     layer.127 + layer.128 + layer.129 + layer.130 + layer.131 + 
    ##     layer.132 + layer.133 + layer.134 + layer.135 + layer.136 + 
    ##     layer.137 + layer.138 + layer.139 + layer.140 + layer.141 + 
    ##     layer.142 + layer.143 + layer.144 + layer.145 + layer.146 + 
    ##     layer.147 + layer.148 + layer.149 + layer.150 + layer.151 + 
    ##     layer.152 + layer.153 + layer.154 + layer.155 + SW_occurrence + 
    ##     dist + LGN51 + LGN52 + LGN53 + LGN54 + LGN55 + LGN56 + LGN59 + 
    ##     LGN511 + LGN516 + LGN518 + LGN519 + LGN520 + LGN523 + LGN526

this is rather long formula since distance to each sampling point is a
separate map. Note also that we use land cover classes `LGN*` as
indicators. We can now produce a regression matrix and fit a ranger
model by using:

``` r
rm.zinc1 <- do.call(cbind, list(meuse@data["zinc"], 
                             over(meuse["zinc"], meuse.grid), 
                             over(meuse["zinc"], grid.dist0)))
```

the model shows that this list of predictors accounts for ca 65% of
variation in Zinc values based on the Out-of-Bag (OOB) samples (Wright &
Ziegler, 2017):

``` r
m1.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, 
                   importance="impurity", seed=1, quantreg= TRUE)
m1.zinc
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(fm1, rm.zinc1, mtry = 22, num.trees = 500, importance = "impurity",      seed = 1, quantreg = TRUE) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  500 
    ## Sample size:                      155 
    ## Number of independent variables:  171 
    ## Mtry:                             22 
    ## Target node size:                 5 
    ## Variable importance mode:         impurity 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       46956.93 
    ## R squared (OOB):                  0.6515079

Note that we set `quantreg=TRUE` which initiates the Quantile Regression
RF approach (Meinshausen, 2006) and helps us estimate also the
prediction error variance i.e. prediction intervals.

Further inspection of the model shows that especially distance to the
river and LGN5 class 16, help with predicting zinc concentrations.
Nevertheless, it seems that buffer distances are most important for
mapping zinc i.e. more important than surface water occurrence, flood
frequency, distance to river and elevation for producing the final
predictions:

``` r
xl <- as.list(ranger::importance(m1.zinc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
```

    ##                [,1]
    ## dist      3829683.2
    ## LGN516    1089382.2
    ## layer.54   979982.3
    ## layer.55   539041.4
    ## layer.53   519892.3
    ## layer.52   432124.7
    ## layer.59   379843.5
    ## layer.82   372548.6
    ## layer.155  344361.8
    ## layer.58   323831.0

We can now predict median value of Zinc and upper and lower 67%
probability prediction intervals (which corresponds to $$1 standard
deviation):

``` r
pred.zinc.rfq = predict(m1.zinc, 
                        cbind(meuse.grid@data, grid.dist0@data), 
                        type="quantiles", quantiles=quantiles)
str(pred.zinc.rfq)
```

    ## List of 5
    ##  $ num.samples              : int 3103
    ##  $ treetype                 : chr "Regression"
    ##  $ num.independent.variables: num 171
    ##  $ num.trees                : num 500
    ##  $ predictions              : num [1:3103, 1:3] 640 612 281 269 472 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:3] "quantile= 0.159" "quantile= 0.5" "quantile= 0.841"
    ##  - attr(*, "class")= chr "ranger.prediction"

where `"quantile= 0.159"` is the lower prediction interval and
`"quantile= 0.841"` is the upper prediction interval. For example for
prediction location 1 the interval is:

``` r
pred.zinc.rfq$predictions[1,]
```

    ## quantile= 0.159   quantile= 0.5 quantile= 0.841 
    ##         640.000        1022.000        1139.295

which shows that the prediction range is relatively wide (note also that
the upper and lower prediction intervals are not necessarily
symetric\!). We can copy the predicted lower and upper intervals to the
spatial object so we can also plot values as maps (maps of predictions
can be found in this [tutorial](https://github.com/thengl/GeoMLA)):

``` r
meuse.grid$zinc_rfq_U = pred.zinc.rfq$predictions[,3]
meuse.grid$zinc_rfq_L = pred.zinc.rfq$predictions[,1]
```

Assuming normal distribution of errors the 67% probability prediction
interval should match 1 s.d. of the prediction error:

![Histogram of s.d. of the prediction error estimated using
QRF.](Workshop_uncertainty_files/figure-gfm/rfq-histogram-1.png)

Compare this numbers with the OOB RMSE and mean s.d. of prediction
error:

``` r
mean(meuse.grid$zinc_rfq_r, na.rm=TRUE); sqrt(m1.zinc$prediction.error)
```

    ## [1] 165.6473

    ## [1] 216.6955

This shows that the mean prediction error is smaller than the OOB RMSE
for about 25%, but in general numbers do match.

An alternative approach to estimating the uncertainty of predictions is
the Jackknifing approach approach (Wager, Hastie, & Efron,
2014):

``` r
m2.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, seed=1, keep.inbag=TRUE)
```

Here the `keep.inbag=TRUE` initaties the Jackknifing approach, which
estimate standard errors of the expected values of predictions, used to
construct confidence intervals. The prediction can be generate by
using:

``` r
pred.zinc.rfj = predict(m2.zinc, cbind(meuse.grid@data, grid.dist0@data), type="se")
str(pred.zinc.rfj)
```

    ## List of 6
    ##  $ predictions              : num [1:3103] 922 925 836 826 888 ...
    ##  $ num.trees                : num 500
    ##  $ num.independent.variables: num 171
    ##  $ num.samples              : int 3103
    ##  $ treetype                 : chr "Regression"
    ##  $ se                       : num [1:3103] 13.6 13.5 13.8 12.9 12.9 ...
    ##  - attr(*, "class")= chr "ranger.prediction"

which adds one extra column called `se` i.e. standard errors. If you
compare OOB RMSE and mean s.d. of prediction error you will notice that
the `se` values are significantly smaller:

``` r
mean(pred.zinc.rfj$se, na.rm=TRUE); sqrt(m2.zinc$prediction.error)
```

    ## [1] 33.41297

    ## [1] 216.6955

We can *scale* the values of `se` so they reflect the mean RMSE by
using:

``` r
meuse.grid$zinc_rfj_r = pred.zinc.rfj$se * 
    sqrt(m2.zinc$prediction.error)/mean(pred.zinc.rfj$se, na.rm=TRUE)
```

If we plot the two maps of errors next to each other we can see
relatively similar patterns:

![Comparison of uncertainty maps based on the QRF vs Jackknife
approaches for the Meuse data
set.](Workshop_uncertainty_files/figure-gfm/jacknife-meuse-maps-1.png)

Note how the single isolated outlier in the lower right corner is
depicted by the RFsp prediction error map. This single isolated high
value of Zinc in that area RFsp has problem explaining as it does not
correlate to any covariates, hence the prediction errors will also be
high.

Another interesting dataset for comparison of RFsp with linear
geostatistical modeling is the Swiss rainfall dataset used in the
Spatial Interpolation Comparison (SIC 1997) exercise, described in
detail in Dubois, Malczewski, & De Cort (2003). This dataset contains
467 measurements of daily rainfall in Switzerland on the 8th of May
1986. Possible covariates include elevation (DEM) and the long term mean
monthly precipitation for May based on the CHELSA climatic images
(Karger et al., 2017) at 1 km:

``` r
sic97.sp = readRDS("../RF_vs_kriging/data/rainfall/sic97.rds")
swiss1km = readRDS("../RF_vs_kriging/data/rainfall/swiss1km.rds")
ov2 = over(y=swiss1km, x=sic97.sp)
```

We can fit a RFsp model for this data set using the same approach from
above. We first derive buffer distances and create the regression
matrix:

``` r
swiss.dist0 <- GSIF::buffer.dist(sic97.sp["rainfall"], 
                                 swiss1km[1], as.factor(1:nrow(sic97.sp))) 
## takes 1+ mins!
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
sw.fm1
```

    ## rainfall ~ layer.1 + layer.2 + layer.3 + layer.4 + layer.5 + 
    ##     layer.6 + layer.7 + layer.8 + layer.9 + layer.10 + layer.11 + 
    ##     layer.12 + layer.13 + layer.14 + layer.15 + layer.16 + layer.17 + 
    ##     layer.18 + layer.19 + layer.20 + layer.21 + layer.22 + layer.23 + 
    ##     layer.24 + layer.25 + layer.26 + layer.27 + layer.28 + layer.29 + 
    ##     layer.30 + layer.31 + layer.32 + layer.33 + layer.34 + layer.35 + 
    ##     layer.36 + layer.37 + layer.38 + layer.39 + layer.40 + layer.41 + 
    ##     layer.42 + layer.43 + layer.44 + layer.45 + layer.46 + layer.47 + 
    ##     layer.48 + layer.49 + layer.50 + layer.51 + layer.52 + layer.53 + 
    ##     layer.54 + layer.55 + layer.56 + layer.57 + layer.58 + layer.59 + 
    ##     layer.60 + layer.61 + layer.62 + layer.63 + layer.64 + layer.65 + 
    ##     layer.66 + layer.67 + layer.68 + layer.69 + layer.70 + layer.71 + 
    ##     layer.72 + layer.73 + layer.74 + layer.75 + layer.76 + layer.77 + 
    ##     layer.78 + layer.79 + layer.80 + layer.81 + layer.82 + layer.83 + 
    ##     layer.84 + layer.85 + layer.86 + layer.87 + layer.88 + layer.89 + 
    ##     layer.90 + layer.91 + layer.92 + layer.93 + layer.94 + layer.95 + 
    ##     layer.96 + layer.97 + layer.98 + layer.99 + layer.100 + layer.101 + 
    ##     layer.102 + layer.103 + layer.104 + layer.105 + layer.106 + 
    ##     layer.107 + layer.108 + layer.109 + layer.110 + layer.111 + 
    ##     layer.112 + layer.113 + layer.114 + layer.115 + layer.116 + 
    ##     layer.117 + layer.118 + layer.119 + layer.120 + layer.121 + 
    ##     layer.122 + layer.123 + layer.124 + layer.125 + layer.126 + 
    ##     layer.127 + layer.128 + layer.129 + layer.130 + layer.131 + 
    ##     layer.132 + layer.133 + layer.134 + layer.135 + layer.136 + 
    ##     layer.137 + layer.138 + layer.139 + layer.140 + layer.141 + 
    ##     layer.142 + layer.143 + layer.144 + layer.145 + layer.146 + 
    ##     layer.147 + layer.148 + layer.149 + layer.150 + layer.151 + 
    ##     layer.152 + layer.153 + layer.154 + layer.155 + layer.156 + 
    ##     layer.157 + layer.158 + layer.159 + layer.160 + layer.161 + 
    ##     layer.162 + layer.163 + layer.164 + layer.165 + layer.166 + 
    ##     layer.167 + layer.168 + layer.169 + layer.170 + layer.171 + 
    ##     layer.172 + layer.173 + layer.174 + layer.175 + layer.176 + 
    ##     layer.177 + layer.178 + layer.179 + layer.180 + layer.181 + 
    ##     layer.182 + layer.183 + layer.184 + layer.185 + layer.186 + 
    ##     layer.187 + layer.188 + layer.189 + layer.190 + layer.191 + 
    ##     layer.192 + layer.193 + layer.194 + layer.195 + layer.196 + 
    ##     layer.197 + layer.198 + layer.199 + layer.200 + layer.201 + 
    ##     layer.202 + layer.203 + layer.204 + layer.205 + layer.206 + 
    ##     layer.207 + layer.208 + layer.209 + layer.210 + layer.211 + 
    ##     layer.212 + layer.213 + layer.214 + layer.215 + layer.216 + 
    ##     layer.217 + layer.218 + layer.219 + layer.220 + layer.221 + 
    ##     layer.222 + layer.223 + layer.224 + layer.225 + layer.226 + 
    ##     layer.227 + layer.228 + layer.229 + layer.230 + layer.231 + 
    ##     layer.232 + layer.233 + layer.234 + layer.235 + layer.236 + 
    ##     layer.237 + layer.238 + layer.239 + layer.240 + layer.241 + 
    ##     layer.242 + layer.243 + layer.244 + layer.245 + layer.246 + 
    ##     layer.247 + layer.248 + layer.249 + layer.250 + layer.251 + 
    ##     layer.252 + layer.253 + layer.254 + layer.255 + layer.256 + 
    ##     layer.257 + layer.258 + layer.259 + layer.260 + layer.261 + 
    ##     layer.262 + layer.263 + layer.264 + layer.265 + layer.266 + 
    ##     layer.267 + layer.268 + layer.269 + layer.270 + layer.271 + 
    ##     layer.272 + layer.273 + layer.274 + layer.275 + layer.276 + 
    ##     layer.277 + layer.278 + layer.279 + layer.280 + layer.281 + 
    ##     layer.282 + layer.283 + layer.284 + layer.285 + layer.286 + 
    ##     layer.287 + layer.288 + layer.289 + layer.290 + layer.291 + 
    ##     layer.292 + layer.293 + layer.294 + layer.295 + layer.296 + 
    ##     layer.297 + layer.298 + layer.299 + layer.300 + layer.301 + 
    ##     layer.302 + layer.303 + layer.304 + layer.305 + layer.306 + 
    ##     layer.307 + layer.308 + layer.309 + layer.310 + layer.311 + 
    ##     layer.312 + layer.313 + layer.314 + layer.315 + layer.316 + 
    ##     layer.317 + layer.318 + layer.319 + layer.320 + layer.321 + 
    ##     layer.322 + layer.323 + layer.324 + layer.325 + layer.326 + 
    ##     layer.327 + layer.328 + layer.329 + layer.330 + layer.331 + 
    ##     layer.332 + layer.333 + layer.334 + layer.335 + layer.336 + 
    ##     layer.337 + layer.338 + layer.339 + layer.340 + layer.341 + 
    ##     layer.342 + layer.343 + layer.344 + layer.345 + layer.346 + 
    ##     layer.347 + layer.348 + layer.349 + layer.350 + layer.351 + 
    ##     layer.352 + layer.353 + layer.354 + layer.355 + layer.356 + 
    ##     layer.357 + layer.358 + layer.359 + layer.360 + layer.361 + 
    ##     layer.362 + layer.363 + layer.364 + layer.365 + layer.366 + 
    ##     layer.367 + layer.368 + layer.369 + layer.370 + layer.371 + 
    ##     layer.372 + layer.373 + layer.374 + layer.375 + layer.376 + 
    ##     layer.377 + layer.378 + layer.379 + layer.380 + layer.381 + 
    ##     layer.382 + layer.383 + layer.384 + layer.385 + layer.386 + 
    ##     layer.387 + layer.388 + layer.389 + layer.390 + layer.391 + 
    ##     layer.392 + layer.393 + layer.394 + layer.395 + layer.396 + 
    ##     layer.397 + layer.398 + layer.399 + layer.400 + layer.401 + 
    ##     layer.402 + layer.403 + layer.404 + layer.405 + layer.406 + 
    ##     layer.407 + layer.408 + layer.409 + layer.410 + layer.411 + 
    ##     layer.412 + layer.413 + layer.414 + layer.415 + layer.416 + 
    ##     layer.417 + layer.418 + layer.419 + layer.420 + layer.421 + 
    ##     layer.422 + layer.423 + layer.424 + layer.425 + layer.426 + 
    ##     layer.427 + layer.428 + layer.429 + layer.430 + layer.431 + 
    ##     layer.432 + layer.433 + layer.434 + layer.435 + layer.436 + 
    ##     layer.437 + layer.438 + layer.439 + layer.440 + layer.441 + 
    ##     layer.442 + layer.443 + layer.444 + layer.445 + layer.446 + 
    ##     layer.447 + layer.448 + layer.449 + layer.450 + layer.451 + 
    ##     layer.452 + layer.453 + layer.454 + layer.455 + layer.456 + 
    ##     layer.457 + layer.458 + layer.459 + layer.460 + layer.461 + 
    ##     layer.462 + layer.463 + layer.464 + layer.465 + layer.466 + 
    ##     layer.467 + CHELSA_rainfall + DEM

``` r
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))
```

We can next fit a RFsp model by using (previously fine-tuned RF
parameters):

``` r
m1.rain <- ranger(sw.fm1, sw.rm[complete.cases(sw.rm),], mtry=27, 
                  min.node.size=2, sample.fraction=0.9930754, 
                  num.trees=150, importance = "impurity", seed=1, quantreg=TRUE)
m1.rain
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(sw.fm1, sw.rm[complete.cases(sw.rm), ], mtry = 27, min.node.size = 2,      sample.fraction = 0.9930754, num.trees = 150, importance = "impurity",      seed = 1, quantreg = TRUE) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  150 
    ## Sample size:                      456 
    ## Number of independent variables:  469 
    ## Mtry:                             27 
    ## Target node size:                 2 
    ## Variable importance mode:         impurity 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       2146.155 
    ## R squared (OOB):                  0.8311812

which shows that the model explains 83% of variation in the daily
precipitation data. Next we predict values and uncertainty by using:

``` r
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), 
                     type="quantiles", quantiles=quantiles)$predictions
## now more computational...
swiss1km$rainfall_rfd1 = rain.rfd1[,2]
## s.d. of the prediction error:
swiss1km$rainfall_rfd1_var = (rain.rfd1[,3]-rain.rfd1[,1])/2
str(swiss1km@data)
```

    ## 'data.frame':    41736 obs. of  5 variables:
    ##  $ CHELSA_rainfall  : int  84 73 70 69 79 99 105 94 92 83 ...
    ##  $ DEM              : num  401 406 395 395 397 ...
    ##  $ border           : chr  "Switzerland" "Switzerland" "Switzerland" "Switzerland" ...
    ##  $ rainfall_rfd1    : num  127 127 127 127 127 127 129 127 127 126 ...
    ##  $ rainfall_rfd1_var: num  12.7 12.2 14 12.7 13.3 ...

this finally gives:

![Predictions and prediction errors for the SIC97 data
set.](Workshop_uncertainty_files/figure-gfm/qrf-sic97-maps-1.png)

the map on the right shows that there are specific zones where the
uncertainty is particulary high. This indicates that the RFsp prediction
error maps are potentially more informative than the geostatistical
error maps (e.g. UK variance map): it can be used to depict local areas
that are significantly more heterogeneous and complex and that require,
either, denser sampling networks or covariates that better represent
local processes in these areas.

![figure](../RF_vs_kriging/img/sic97_rainfall_qgis.jpg) *Figure: SIC97
data set predictions based on RFsp visualized in QGIS.*

So in summary: uncertainty of predictions in RF models can be
efficiently estimated using either the QRF (prediction intervals) or the
Jacknifing approaches (confidence intervals), both are available via the
the ranger package (Wright & Ziegler, 2017). Prediction intervals should
in average match the RMSE estimated using OOB samples or other
Cross-validation approaches. Note however that both approaches are
computationally intensive and could increase the prediction time at the
order of magnitude times. There are additional costs to pay to derive a
reliable and detailed measures of uncertainty.

## Mapping prediction errors for factor/binomial variables

In the next example we look at mapping the uncertainty of predictions of
factor/binomial variables. Consider for example the soil type classes in
the Meuse data set:

``` r
summary(meuse$soil)
```

    ##  1  2  3 
    ## 97 46 12

We can first look at mapping the occurrence of the class
`"1"`:

``` r
meuse@data = cbind(meuse@data, data.frame(model.matrix(~soil-1, meuse@data)))
summary(as.factor(meuse$soil1))
```

    ##  0  1 
    ## 58 97

To produce a map of `soil1` using RFsp we have now two options:

  - *Option 1*: treat binomial variable as numeric variable with 0 / 1
    values (thus a regression problem),
  - *Option 2*: treat binomial variable as factor variable with a single
    class (thus a classification problem),

Both methods are in fact equivalent and should give the same
predictions. There will be however some differences in how are the
uncertainty maps estimated. In the case of *Option 1* we fit a similar
type of model as in the previous example:

``` r
fm.s1 = as.formula(paste("soil1 ~ ", 
                  paste(names(grid.dist0), collapse="+"), 
                  " + SW_occurrence + dist"))
fm.s1
```

    ## soil1 ~ layer.1 + layer.2 + layer.3 + layer.4 + layer.5 + layer.6 + 
    ##     layer.7 + layer.8 + layer.9 + layer.10 + layer.11 + layer.12 + 
    ##     layer.13 + layer.14 + layer.15 + layer.16 + layer.17 + layer.18 + 
    ##     layer.19 + layer.20 + layer.21 + layer.22 + layer.23 + layer.24 + 
    ##     layer.25 + layer.26 + layer.27 + layer.28 + layer.29 + layer.30 + 
    ##     layer.31 + layer.32 + layer.33 + layer.34 + layer.35 + layer.36 + 
    ##     layer.37 + layer.38 + layer.39 + layer.40 + layer.41 + layer.42 + 
    ##     layer.43 + layer.44 + layer.45 + layer.46 + layer.47 + layer.48 + 
    ##     layer.49 + layer.50 + layer.51 + layer.52 + layer.53 + layer.54 + 
    ##     layer.55 + layer.56 + layer.57 + layer.58 + layer.59 + layer.60 + 
    ##     layer.61 + layer.62 + layer.63 + layer.64 + layer.65 + layer.66 + 
    ##     layer.67 + layer.68 + layer.69 + layer.70 + layer.71 + layer.72 + 
    ##     layer.73 + layer.74 + layer.75 + layer.76 + layer.77 + layer.78 + 
    ##     layer.79 + layer.80 + layer.81 + layer.82 + layer.83 + layer.84 + 
    ##     layer.85 + layer.86 + layer.87 + layer.88 + layer.89 + layer.90 + 
    ##     layer.91 + layer.92 + layer.93 + layer.94 + layer.95 + layer.96 + 
    ##     layer.97 + layer.98 + layer.99 + layer.100 + layer.101 + 
    ##     layer.102 + layer.103 + layer.104 + layer.105 + layer.106 + 
    ##     layer.107 + layer.108 + layer.109 + layer.110 + layer.111 + 
    ##     layer.112 + layer.113 + layer.114 + layer.115 + layer.116 + 
    ##     layer.117 + layer.118 + layer.119 + layer.120 + layer.121 + 
    ##     layer.122 + layer.123 + layer.124 + layer.125 + layer.126 + 
    ##     layer.127 + layer.128 + layer.129 + layer.130 + layer.131 + 
    ##     layer.132 + layer.133 + layer.134 + layer.135 + layer.136 + 
    ##     layer.137 + layer.138 + layer.139 + layer.140 + layer.141 + 
    ##     layer.142 + layer.143 + layer.144 + layer.145 + layer.146 + 
    ##     layer.147 + layer.148 + layer.149 + layer.150 + layer.151 + 
    ##     layer.152 + layer.153 + layer.154 + layer.155 + SW_occurrence + 
    ##     dist

``` r
rm.s1 <- do.call(cbind, list(meuse@data["soil1"], 
                             over(meuse["soil1"], meuse.grid), 
                             over(meuse["soil1"], grid.dist0)))
```

which
gives:

``` r
m1.s1 <- ranger(fm.s1, rm.s1, mtry=22, num.trees=500, seed = 1, quantreg=TRUE)
m1.s1
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(fm.s1, rm.s1, mtry = 22, num.trees = 500, seed = 1, quantreg = TRUE) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  500 
    ## Sample size:                      155 
    ## Number of independent variables:  157 
    ## Mtry:                             22 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.05673872 
    ## R squared (OOB):                  0.7592689

In the case of *Option 2* we treat the binomial variable as factor
variable

``` r
rm.s1$soil1c = as.factor(rm.s1$soil1)
summary(rm.s1$soil1c)
```

    ##  0  1 
    ## 58 97

and the model turns into a classification
problem:

``` r
fm.s1c <- as.formula(paste("soil1c ~ ", paste(names(grid.dist0), collapse="+"), 
                           " + SW_occurrence + dist"))
m2.s1 <- ranger(fm.s1c, rm.s1, mtry=22, num.trees=500, 
                seed=1, probability=TRUE, keep.inbag=TRUE)
m2.s1
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(fm.s1c, rm.s1, mtry = 22, num.trees = 500, seed = 1, probability = TRUE,      keep.inbag = TRUE) 
    ## 
    ## Type:                             Probability estimation 
    ## Number of trees:                  500 
    ## Sample size:                      155 
    ## Number of independent variables:  157 
    ## Mtry:                             22 
    ## Target node size:                 10 
    ## Variable importance mode:         none 
    ## Splitrule:                        gini 
    ## OOB prediction error (Brier s.):  0.05711483

which shows that the Out of Bag prediction error (classification error)
is only 0.06 (this number is in the probability scale). Note that, it is
not easy to compare the results of the regression and classification OOB
errors as these are conceptually different. Also note that we turn on
`keep.inbag = TRUE` so that ranger can estimate the classification
errors using the Jackknife-after-Bootstrap method (Wager et al., 2014).
`quantreg=TRUE` obviously would not work here since it is a
classification and not a regression problem.

We next derive prediction errors for the two options:

``` r
pred.soil1_rfb = predict(m1.s1, 
                         cbind(meuse.grid@data, grid.dist0@data), 
                         type="quantiles", quantiles=quantiles)
str(pred.soil1_rfb)
```

    ## List of 5
    ##  $ num.samples              : int 3103
    ##  $ treetype                 : chr "Regression"
    ##  $ num.independent.variables: num 157
    ##  $ num.trees                : num 500
    ##  $ predictions              : num [1:3103, 1:3] 0 0 0 0 0 0 0 0 1 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:3] "quantile= 0.159" "quantile= 0.5" "quantile= 0.841"
    ##  - attr(*, "class")= chr "ranger.prediction"

and copy the values to the spatial object:

``` r
meuse.grid$soil1_rfq_U = pred.soil1_rfb$predictions[,3]
meuse.grid$soil1_rfq = pred.soil1_rfb$predictions[,2]
meuse.grid$soil1_rfq_L = pred.soil1_rfb$predictions[,1]
meuse.grid$soil1_rfq_r = (meuse.grid$soil1_rfq_U - meuse.grid$soil1_rfq_L)/2
mean(meuse.grid$soil1_rfq_r, na.rm=TRUE); sqrt(m1.s1$prediction.error)
```

    ## [1] 0.1197154

    ## [1] 0.2381989

Again, QRF error estimates are somewhat smaller than the OOB RMSE. We
derive the errors also using the Jacknifing
approach:

``` r
pred.soil1_rfc = predict(m2.s1, cbind(meuse.grid@data, grid.dist0@data), type="se")
meuse.grid$soil1_rfc = pred.soil1_rfc$predictions[,2]
```

which can be scaled to RMSE using:

``` r
meuse.grid$soil1_rfc_r = pred.soil1_rfc$se[,2] *
  sqrt(m2.s1$prediction.error)/mean(pred.soil1_rfc$se[,2], na.rm=TRUE)
```

We also derive predictions of the binomial variable using the mean
estimate (which can be often different from the median
estimate\!):

``` r
pred.regr <- predict(m1.s1, cbind(meuse.grid@data, grid.dist0@data), type="response")$predictions
meuse.grid$soil1_rfr <- pred.regr
```

We can finally plot all maps next to each other by using:

![Comparison of uncertainty maps for a binomial variable based on the
QRF vs Jackknife approaches for the Meuse data
set.](Workshop_uncertainty_files/figure-gfm/binomial-meuse-maps-1.png)

Again both approaches to estimating mapping uncertainty give similar
patterns: basically transition zones between class 1 and other soil
types are most uncertain. Note however that the QRF estimate of the
uncertainty (maps on the left) shows much more distinct jumps in values.

## Mapping prediction errors for a factor variable

In the last example we look at how to map uncertainty of predicting
factor-type variable:

``` r
summary(meuse$soil)
```

    ##  1  2  3 
    ## 97 46 12

which has three classes and hence we treat this variable as a
classification problem, which also means that we can now only use the
Jacknifing
approach:

``` r
fm.s = as.formula(paste("soil ~ ", paste(names(grid.dist0), collapse="+"), 
                        " + SW_occurrence + dist"))
rm.s <- do.call(cbind, list(meuse@data["soil"], 
                            over(meuse["soil"], meuse.grid), 
                            over(meuse["soil"], grid.dist0)))
m.s <- ranger(fm.s, rm.s, mtry=22, num.trees=500, seed=1, probability=TRUE, keep.inbag=TRUE)
m.s
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(fm.s, rm.s, mtry = 22, num.trees = 500, seed = 1, probability = TRUE,      keep.inbag = TRUE) 
    ## 
    ## Type:                             Probability estimation 
    ## Number of trees:                  500 
    ## Sample size:                      155 
    ## Number of independent variables:  157 
    ## Mtry:                             22 
    ## Target node size:                 10 
    ## Variable importance mode:         none 
    ## Splitrule:                        gini 
    ## OOB prediction error (Brier s.):  0.08903655

this shows that the model is succesful with the OOB prediction error of
about 0.09. This number is rather abstract so we can also check what is
the actual classification accuracy using hard classes:

``` r
m.s0 <- ranger(fm.s, rm.s, mtry=22, num.trees=150, seed=1)
m.s0
```

    ## Ranger result
    ## 
    ## Call:
    ##  ranger(fm.s, rm.s, mtry = 22, num.trees = 150, seed = 1) 
    ## 
    ## Type:                             Classification 
    ## Number of trees:                  150 
    ## Sample size:                      155 
    ## Number of independent variables:  157 
    ## Mtry:                             22 
    ## Target node size:                 1 
    ## Variable importance mode:         none 
    ## Splitrule:                        gini 
    ## OOB prediction error:             10.32 %

which shows that the classification or mapping accuracy for hard classes
is about 90%. We can produce predictions and uncertainty maps by
using:

``` r
pred.soil_rfc = predict(m.s, cbind(meuse.grid@data, grid.dist0@data), type="se")
```

To plot the prediction we can copy the data to a new object:

``` r
pred.grids = meuse.grid["soil"]
pred.grids@data = do.call(cbind, list(pred.grids@data, 
                                      data.frame(pred.soil_rfc$predictions),
                                      data.frame(pred.soil_rfc$se)))
names(pred.grids) = c("soil", paste0("pred_soil", 1:3), paste0("se_soil", 1:3))
str(pred.grids@data)
```

    ## 'data.frame':    3103 obs. of  7 variables:
    ##  $ soil      : Factor w/ 3 levels "1","2","3": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ pred_soil1: num  0.757 0.759 0.753 0.74 0.765 ...
    ##  $ pred_soil2: num  0.212 0.218 0.22 0.229 0.212 ...
    ##  $ pred_soil3: num  0.0307 0.0232 0.0272 0.0307 0.0232 ...
    ##  $ se_soil1  : num  0.1004 0.1093 0.1061 0.0939 0.072 ...
    ##  $ se_soil2  : num  0.0804 0.0836 0.0821 0.0803 0.0678 ...
    ##  $ se_soil3  : num  0.0331 0.0329 0.0331 0.0331 0.0329 ...

which gives 6 columns in total: 3 columns for predictions and 3 columns
for prediction errors `se`. We can plot the three maps next to each
other by using:

![Predictions of soil types for the meuse data set based on the RFsp:
(above) probability for three soil classes, and (below) derived standard
errors per
class.](Workshop_uncertainty_files/figure-gfm/factor-meuse-maps-1.png)

This shows that some individual points seem to be problematic as they
show high errors for at least 2 classes, and also the lower right corner
of the study area is in average the most difficult to map.

## Spatial cross-validation

When fitting spatial prediction models one should try as much as
possible to use the so-called [spatial
cross-validation](https://mlr-org.github.io/mlr/articles/tutorial/devel/handling_of_spatial_data.html)
to estimate a realistic measure of prediction accuracy. The spatial
cross-validation can be run by using for example the mlr package:

``` r
library(mlr)
```

    ## Loading required package: ParamHelpers

    ## 
    ## Attaching package: 'ParamHelpers'

    ## The following object is masked from 'package:raster':
    ## 
    ##     getValues

    ## 
    ## Attaching package: 'mlr'

    ## The following object is masked from 'package:raster':
    ## 
    ##     resample

``` r
spatial.taskmeuse = makeRegrTask(data = rm.zinc1[,c("zinc","SW_occurrence","dist","AHN")], target = "zinc", coordinates = data.frame(meuse@coords))
spatial.taskmeuse
```

    ## Supervised task: rm.zinc1[, c("zinc", "SW_occurrence", "dist", "AHN")]
    ## Type: regr
    ## Target: zinc
    ## Observations: 155
    ## Features:
    ##    numerics     factors     ordered functionals 
    ##           3           0           0           0 
    ## Missings: FALSE
    ## Has weights: FALSE
    ## Has blocking: FALSE
    ## Has coordinates: TRUE

``` r
learner.rf = makeLearner("regr.ranger")
library("parallelMap")
parallelStartSocket(parallel::detectCores())
```

    ## Starting parallelization in mode=socket with cpus=8.

``` r
resampling = makeResampleDesc("SpRepCV", fold = 5, reps = 5)
cv.meuse = mlr::resample(learner = learner.rf, task = spatial.taskmeuse, resampling = resampling)
```

    ## Exporting objects to slaves for mode socket: .mlr.slave.options

    ## Resampling: repeated spatial cross-validation

    ## Measures:             mse

    ## Mapping in parallel: mode = socket; cpus = 8; elements = 25.

    ## 

    ## Aggregated Result: mse.test.mean=73347.5505395

    ## 

``` r
## compare with non-spatial CV:
nonspatial.taskmeuse = makeRegrTask(data = rm.zinc1[,c("zinc","SW_occurrence","dist","AHN")], target = "zinc")
resampling0 = makeResampleDesc("RepCV", fold = 5, reps = 5)
cv.meuse0 = mlr::resample(learner = learner.rf, task = nonspatial.taskmeuse, resampling = resampling0)
```

    ## Exporting objects to slaves for mode socket: .mlr.slave.options

    ## Resampling: repeated cross-validation

    ## Measures:             mse

    ## Mapping in parallel: mode = socket; cpus = 8; elements = 25.

    ## 

    ## Aggregated Result: mse.test.mean=63021.7024207

    ## 

``` r
parallelStop()
```

    ## Stopped parallelization. All cleaned up.

This shows that the non-spatial CV (as expected) will result in about
10% lower RMSE (253 vs 264). The difference between the spatial vs
non-spatial CV will likely be more significant if the points are
spatially clustered and a strong spatial auto-correlation exists.

## Summary points

The uncertainty of the predictions of random forest for regression-type
problems can be estimated using several approaches:

  - The Jackknife-after-Bootstrap method (see e.g. Wager et al. (2014)).

  - The U-statistics approach of Mentch & Hooker (2016).

  - The Monte Carlo simulations (both target variable and covariates)
    approach of Coulston, Blinn, Thomas, & Wynne (2016).

  - The Quantile Regression Forests (QRF) method (Meinshausen, 2006).

The approaches by Wager et al. (2014) and Mentch & Hooker (2016)
estimate standard errors of the expected values of predictions, used to
construct confidence intervals, while the approaches of Coulston et al.
(2016) and Meinshausen (2006) estimate prediction intervals. The
Quantile Regression Forests (QRF) algorithm estimates the quantiles of
the distribution of the target variable at prediction points. Thus, the
0.025 and 0.975 quantile may be used to derive the lower and upper
limits of a symmetric prediction interval.

In summary: spatial prediction of uncertainty for numeric, binomial and
factor-type variables is straight forward with ranger: buffer distance
and spatial-autocorrelation can be incorporated at once. Compare with
geostatistical packages where GLMs with logit link function and/or
indicator kriging would need to be used, and which requires that
variograms are fitted per class.

The QRF and Jacknifing approaches give about similar predictions of the
uncertainty (similar patterns) but refer to different scales. The QRF
approach seems to be somewhat more attractive as it can be used to
produce whole distributions of prediction errors (by setting the
`quantiles` argument). One should also be careful about the difference
between estimating median and mean values (they can often be
different\!). Also, both the QRF and the Jacknifing approaches can be
computational and increase the computational load by order of mangitude,
hence plan carefully your analysis and prediction with larger data sets.

## References

<div id="refs" class="references">

<div id="ref-Biau2016">

Biau, G., & Scornet, E. (2016). A random forest guided tour. *TEST*,
*25*(2), 197–227.
doi:[10.1007/s11749-016-0481-7](https://doi.org/10.1007/s11749-016-0481-7)

</div>

<div id="ref-COULSTON2016189">

Coulston, J. W., Blinn, C. E., Thomas, V. A., & Wynne, R. H. (2016).
Approximating prediction uncertainty for random forest regression
models. *Photogrammetric Engineering & Remote Sensing*, *82*(3),
189–197.
doi:[10.14358/PERS.82.3.189](https://doi.org/10.14358/PERS.82.3.189)

</div>

<div id="ref-dubois2003mapping">

Dubois, G., Malczewski, J., & De Cort, M. (2003). *Mapping Radioactivity
in the Environment: Spatial Interpolation Comparison 97*. Office for
Official Publications of the European Communities.

</div>

<div id="ref-karger2017climatologies">

Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H.,
Soria-Auza, R. W., … Kessler, M. (2017). Climatologies at high
resolution for the earth’s land surface areas. *Scientific Data*, *4*.

</div>

<div id="ref-meinshausen2006quantile">

Meinshausen, N. (2006). Quantile regression forests. *Journal of Machine
Learning Research*, *7*, 983–999.

</div>

<div id="ref-mentch2016quantifying">

Mentch, L., & Hooker, G. (2016). Quantifying uncertainty in random
forests via confidence intervals and hypothesis tests. *Journal of
Machine Learning Research*, *17*(1), 841–881.

</div>

<div id="ref-nussbaum2018evaluation">

Nussbaum, M., Spiess, K., Baltensweiler, A., Grob, U., Keller, A.,
Greiner, L., … Papritz, A. (2018). Evaluation of digital soil mapping
approaches with large sets of environmental covariates. *Soil*, *4*(1),
1.

</div>

<div id="ref-pekel2016high">

Pekel, J.-F., Cottam, A., Gorelick, N., & Belward, A. S. (2016).
High-resolution mapping of global surface water and its long-term
changes. *Nature*, *504*, 418–422.

</div>

<div id="ref-prasad2006newer">

Prasad, A. M., Iverson, L. R., & Liaw, A. (2006). Newer classification
and regression tree techniques: Bagging and random forests for
ecological prediction. *Ecosystems*, *9*(2), 181–199.

</div>

<div id="ref-wager2014confidence">

Wager, S., Hastie, T., & Efron, B. (2014). Confidence intervals for
random forests: The jackknife and the infinitesimal jackknife. *Journal
of Machine Learning Research*, *15*(1), 1625–1651.

</div>

<div id="ref-wright2017ranger">

Wright, M. N., & Ziegler, A. (2017). ranger: A Fast Implementation of
Random Forests for High Dimensional Data in C++ and R. *Journal of
Statistical Software*, *77*(1), 1–17.

</div>

</div>
