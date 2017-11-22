## Comparison geoR and gstat
## https://stat.ethz.ch/pipermail/r-sig-geo/2017-November/026196.html

library(sp)
library(gstat)
library(geoR)
data(meuse)
dataset= meuse
set.seed(999)

# Split Meuse Dataset into Tzincing and HoldOut Sample datasets
Tzincing_ids <- sample(seq_len(nrow(dataset)), size = (0.7* nrow(dataset)))
Tzincing_sample = dataset[Tzincing_ids,]
Holdout_sample_allvars = dataset[-Tzincing_ids,]

holdoutvars_df <-(dataset[,which(names(dataset) %in% c("x","y","lead","copper","elev","dist"))])
Hold_out_sample = holdoutvars_df[-Tzincing_ids,]
coordinates(Tzincing_sample) <- c('x','y')
coordinates(Hold_out_sample) <- c('x','y')
# Semivariogram modeling
m1  <- variogram(log(zinc)~lead+copper+elev+dist, Tzincing_sample)
m <- vgm("Exp")
m <- fit.variogram(m1, m)

# Apply Univ Krig to Tzincing dataset
prediction_tzincing_data <- krige(log(zinc)~lead+copper+elev+dist, Tzincing_sample, Tzincing_sample, model = m)
prediction_tzincing_data <- expm1(prediction_tzincing_data$var1.pred)

# Apply Univ Krig to Hold Out dataset
prediction_holdout_data <- krige(log(zinc)~lead+copper+elev+dist, Tzincing_sample, Hold_out_sample, model = m)
prediction_holdout_data <- expm1(prediction_holdout_data$var1.pred)

# Computing Predictive errors for Tzincing and Hold Out samples respectively
tzincing_prediction_error_term <- Tzincing_sample$zinc - prediction_tzincing_data
holdout_prediction_error_term <- Holdout_sample_allvars$zinc - prediction_holdout_data

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

# Mean Absolute Error metric :
# UK Predictive errors for Tzincing sample set , and UK Predictive Errors for HoldOut sample set
print(mae(tzincing_prediction_error_term)) #Error for Tzincing sample set
print(mae(holdout_prediction_error_term)) #Error for Hold out sample set

## Compare with how it is done with geoR:
library(geoR)
demo(meuse, echo=FALSE)
set.seed(999)
sel.d = complete.cases(meuse@data[,c("lead","copper","elev", "dist")])
meuse = meuse[sel.d,]
meuse.geo <- as.geodata(meuse["zinc"])
## add covariates:
meuse.geo$covariate = meuse@data[,c("lead","copper","elev", "dist")]
trend = ~ lead+copper+elev+dist
zinc.vgm <- likfit(meuse.geo, lambda=0, trend = trend, ini=c(var(log1p(meuse.geo$data)),800), fix.psiA = FALSE, fix.psiR = FALSE)
zinc.vgm
locs2 = meuse@coords
KC = krige.control(trend.d = trend, trend.l = ~ meuse$lead+meuse$copper+meuse$elev+meuse$dist, obj.model = zinc.vgm)
zinc.uk <- krige.conv(meuse.geo, locations=locs2, krige=KC)

