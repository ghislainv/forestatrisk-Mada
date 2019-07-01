#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# Environmental variables
Sys.unsetenv("DISPLAY") # Remove DISPLAY for Python plot

# Libraries
pkg <- c("coda","reticulate","glue","curl","dplyr","readr","broom",
		 "ggplot2","raster","rasterVis","rgdal","leaflet","kableExtra")
load.pkg <- function(x) {
	if(!require(x, character.only=TRUE)) {
		install.packages(x)
		require(x, character.only=TRUE)
	}
}
loaded <- lapply(pkg,load.pkg)
## Remove useless objects
rm(pkg,load.pkg,loaded)

# Source R functions
source("R/far_plot.R")
source("R/perf.R")

# ================================
# Setup Python virtual environment
# ================================

# Use Python 3
use_python("/usr/bin/python3.5")
# Check python version and virtualenv
py_config()

# Import Python modules
far <- import("forestatrisk")
patsy <- import("patsy")
sm <- import("statsmodels.api")
smf <- import("statsmodels.formula.api")

# ===============================
# Download data
# ===============================

# Make data directory
if (!dir.exists("data")) {
	dir.create("data")
	dir.create("data/model")
	dir.create("data/mada")
	dir.create("data/validation")
}


## ----plot_fcc------------------------------------------------------------
# Make output directory
if (!dir.exists("output")) {
	dir.create("output")
}
# Plot forest cover change 2000-2010
fig <- far$plot$fcc(input_fcc_raster="data/models/2000-2010/fordefor.tif",
					output_file="output/fcc2010.png",
					col=c(255,0,0,255),  # rgba color for deforestation
					figsize=c(5,5),
					dpi=150,
					zoom=c(340000,412000,7420000,7500000))

# ========================================================
# Sample points
# ========================================================

# Training data-set
if (!file.exists("output/sample-2000-2010.txt")) {
  samp <- far$sample(nsamp=20000L, Seed=1234L, csize=10L,
                     var_dir="data/models/2000-2010",
                     input_forest_raster="fordefor.tif",
                     output_file="output/sample-2000-2010.txt",
                     blk_rows=1L)
}
samp <- read.table("output/sample-2000-2010.txt", header=TRUE, sep=",")
set.seed(1234)
train <- sample(1:40000, size=20000, replace=FALSE)
data_train <- samp[train,] %>% dplyr::filter(complete.cases(.))
data_valid <- samp[-train,] %>% dplyr::filter(complete.cases(.))

## ----visu_sampling-------------------------------------------------------
head(data_train)

## ----plot_sample---------------------------------------------------------
# Plot sample points
fig <- far$plot$obs(sample=data_train,
             name_forest_var="fordefor",
             input_fcc_raster="data/models/2000-2010/fordefor.tif",
             output_file="output/obs.png",
             zoom=c(340000,412000,7420000,7500000),
             figsize=c(5,5), #c(11.69,8.27),
             s=5,dpi=300)

## ----correlations--------------------------------------------------------
# Descriptive statistics

# Model formulas
formula_1 <- paste0("fordefor ~ dist_road + dist_town + dist_defor +",
                    "dist_river + dist_edge + altitude + slope - 1")
# Standardized variables (mean=0, std=1)
formula_2 <- paste0("fordefor ~ scale(dist_road) + scale(dist_town) +",
                    "scale(dist_defor) + scale(dist_river) + scale(dist_edge) +",
                    "scale(altitude) + scale(slope) - 1")
formulas <- c(formula_1, formula_2)

# List to store figures
corr_fig <- list()

# Loop on formulas
for (f in 1:length(formulas)) {
    # Output file
    of <- glue("output/correlation_{f}.pdf")
    # Data
    dmat <- patsy$dmatrices(formulas[f], data=data_train, eval_env=-1L,
                                  return_type="dataframe")
    # Plots
    fig <- far$plot$correlation(y=dmat[[1]],data=dmat[[2]],
                         plots_per_page=3L,figsize=c(7,8),
                         dpi=100L,output_file=of)
}

# ========================================================
# hSDM model
# ========================================================

## ----spatial_cells-------------------------------------------------------
# Spatial cells for spatial-autocorrelation
neighborhood <- far$cellneigh_ctry(raster="data/models/2000-2010/fordefor.tif",
                                   vector="data/mada/mada38s.shp",
                                   csize=10L, rank=1L)
nneigh <- neighborhood[[1]]
adj <- neighborhood[[2]]
cell_in <- neighborhood[[3]]
ncell <- neighborhood[[4]]

# Udpate cell number in training dataset
cell_rank <- vector()
for (i in 1:nrow(data_train)) {
  cell_rank[i] <- which(cell_in==data_train$cell[i])-1 # ! cells start at zero
}
data_train$cell <- cell_rank

# Udpate cell number in validation dataset
cell_rank <- vector()
for (i in 1:nrow(data_valid)) {
  cell_rank[i] <- which(cell_in==data_valid$cell[i])-1 # ! cells start at zero
}
data_valid$cell <- cell_rank

# Save data-sets
write.table(data_train, "output/data_train.txt", row.names=FALSE, sep=",")
write.table(data_valid, "output/data_valid.txt", row.names=FALSE, sep=",")

## ----formula-------------------------------------------------------------
# Formula
data_train$trials <- 1  # Set number of trials to one
formula <- paste0("I(1-fordefor) + trials ~ C(sapm) + scale(altitude) +
                  scale(slope) + scale(dist_edge) +
                  scale(dist_defor) +
                  scale(dist_road) + scale(dist_town) + cell")
# formula <- paste0("I(1-fordefor) + trials ~ C(sapm) + scale(altitude) +
#                   scale(slope) + scale(dist_edge) +
#                   scale(dist_defor) + np.power(scale(dist_defor),2) +
#                   scale(dist_road) + scale(dist_town) + cell")

## ----mod_icar------------------------------------------------------------
# Model
mod_icar <- far$model_binomial_iCAR(
  # Observations
  suitability_formula=formula, data=data_train,
  # Spatial structure
  n_neighbors=np_array(nneigh,dtype="int32"), neighbors=np_array(adj,dtype="int32"),
  # Environment
  eval_env=-1L,
  # Chains
  burnin=2000L, mcmc=5000L, thin=5L,
  # Starting values
  beta_start=-99)

## ----summary_icar--------------------------------------------------------
sink(file="output/summary_mod_icar.txt")
print(mod_icar)
sink()
print(mod_icar)

## ----plot_with_py, echo=FALSE--------------------------------------------
traces_fig <- mod_icar$plot(output_file="output/mcmc.pdf",
                            plots_per_page=3L,
                            figsize=c(9,6),
                            dpi=100)

## ----plot_with_r---------------------------------------------------------
require(coda)
mcmc <- as.mcmc(mod_icar$mcmc)
pdf(file="output/mcmc_R.pdf")
plot(mcmc[,c(1:3,9,10)])
dev.off()

## ----rho, fig.height=7, message=FALSE, warning=FALSE---------------------
# Get the spatial random effects
rho <- rep(-9999,ncell)  # -9999 will be considered as nodata
rho[cell_in+1] <- mod_icar$rho

# Resample them
fig <- far$interpolate_rho(rho=r_to_py(rho),
						   input_raster="data/models/2000-2010/fordefor.tif",
                 		   output_file="output/rho.tif",
                 		   csize_orig=10L, csize_new=1L)

# Plot random effects
fig <- far$plot$rho("output/rho_orig.tif",output_file="output/rho_orig.png")
fig <- far$plot$rho("output/rho.tif",output_file="output/rho.png")

# Plot with R
mada <- rgdal::readOGR(dsn="data/mada",layer="mada38s", verbose=FALSE)
r.rho_orig <- raster("output/rho_orig.tif")
r.rho <- raster("output/rho.tif")
rho_plot(r.rho_orig, mada, output_file="output/rho_orig_ggplot.png",
         quantiles_legend=c(0.025,0.975),width=4.5, height=8)
rho_plot(r.rho, mada, output_file="output/rho_ggplot.png",
         quantiles_legend=c(0.025,0.975),width=4.5, height=8)

# ========================================================
# Model comparison
# ========================================================

# Null model
formula_null <- "I(1-fordefor) ~ 1"
dmat_null <- patsy$dmatrices(formula_null, data=r_to_py(data_train), NA_action="drop",
							 return_type="dataframe", eval_env=-1L)
Y <- dmat_null[[1]]
X_null <- dmat_null[[2]]
mod_null <- sm$GLM(Y, X_null, family=sm$families$Binomial())$fit()
print(mod_null$summary())

# Simple glm with no spatial random effects
formula_glm <- paste0("I(1-fordefor) ~ C(sapm) + scale(altitude) + ",
					  "scale(slope) + scale(dist_defor) + ",
					  "scale(dist_edge) + scale(dist_road) + scale(dist_town)")
mod_glm <- smf$glm(formula_glm, r_to_py(data_train),
				   family=sm$families$Binomial(), eval_env=-1L)$fit()

# Summary glm
sink(file="output/summary_mod_binomial_glm.txt")
print(mod_glm$summary())
sink()
print(mod_glm$summary())

# Deviances
deviance_null <- mod_null$deviance
deviance_glm <- mod_glm$deviance
deviance_icar <- mod_icar$deviance
deviance_full <- 0
# Table
dev <- c(deviance_null, deviance_glm, deviance_icar, deviance_full)
mod_dev <- data.frame(model=c("null", "glm", "icar", "full"),
					  deviance=round(dev))
perc <- 100*(1-mod_dev$deviance/deviance_null)
mod_dev$perc <- round(perc)
write.table(mod_dev,file="output/deviance_model_comparison.txt",sep=",",row.names=FALSE)
mod_dev

## ----perf_models---------------------------------------------------------
source("R/perf.R")
set.seed(1234)
performance_null <- performance_index(data_valid,runif(nrow(data_valid)))
performance_nochange <- performance_index(data_valid,rep(0,nrow(data_valid)))
theta_pred <- rep(0,nrow(data_valid))
performance_glm <- performance_index(data_valid,mod_glm$predict(data_valid))
performance_icar <-performance_index(data_valid,mod_icar$predict(data_valid))
# Save
write.table(performance_null, file="output/performance_null.txt",
			sep="\t", row.names=FALSE)
write.table(performance_nochange, file="output/performance_nochange.txt",
			sep="\t", row.names=FALSE)
write.table(performance_glm, file="output/performance_glm.txt",
			sep="\t", row.names=FALSE)
write.table(performance_icar, file="output/performance_icar.txt",
			sep="\t", row.names=FALSE)
# Save for perc=10%
df_mod_valid <- rbind(performance_null[3,],performance_nochange[3,],
					  performance_glm[3,],
					  performance_icar[3,])
df_mod_valid <- df_mod_valid %>%
	dplyr::mutate(mod=c("null","nochange","glm","icar")) %>%
	dplyr::select(10,2:9)
write.table(df_mod_valid, file="output/df_mod_valid.txt", sep="\t", row.names=FALSE)

# Print
df_mod_valid

# ========================================================
# Spatial probability of deforestation in 2010
# ========================================================

## Number of forest pixels in 2000 and 2010
nfor2010_1 <- far$countpix("data/models/2000-2010/fordefor.tif",
						   value=1L, blk_rows=128L)
nfor2010_0 <- far$countpix("data/models/2000-2010/fordefor.tif",
						   value=0L, blk_rows=128L)
nfor2010_npix <- nfor2010_1$npix
nfor2010_ha <- nfor2010_1$area
nfor2000_npix <- nfor2010_1$npix + nfor2010_0$npix
nfor2000_ha <- nfor2010_1$area + nfor2010_0$area

## ----icar probabilities--------------------------------------------------
if (!file.exists("output/prob_icar.tif")) {
  far$predict_raster_binomial_iCAR(
    mod_icar, var_dir="data/models/2010-2017",
    input_cell_raster="output/rho.tif",
    input_forest_raster="data/forest/for2010.tif",
    output_file="output/prob_icar.tif",
    blk_rows=128L)
}

## ----plot icar proba----------------------------------------------------
if (!file.exists("output/prob_icar.png")) {
	fig <- far$plot$prob(
		"output/prob_icar.tif",
		output_file="output/prob_icar.png",
		figsize=c(4,4))
}

## ----glm probabilities------------------------------------------------
if (!file.exists("output/prob_glm.tif")) {
	fig_prob <- far$predict_raster(
		mod_glm, var_dir="data/models/2010-2017",
		input_forest_raster="data/forest/for2010.tif",
		output_file="output/prob_glm.tif",
		blk_rows=128L, transform=TRUE)
}

## ----plot glm proba----------------------------------------------------
if (!file.exists("output/prob_glm.png")) {
	fig <- far$plot$prob(
		"output/prob_glm.tif",
		output_file="output/prob_glm.png",
		figsize=c(4,4))
}

# ========================================================
# Projections in 2017
# ========================================================

## ----defor_area_2010_2017, results="hide"--------------------------------
# Number of deforested pixels
count_d <- far$countpix("data/forest/2010-2017/fordefor.tif", value=0L, blk_rows=128L)

## ----count_d-------------------------------------------------------------
defor_area_2010_2017 <- count_d$area

## ----proj2017, results="hide"--------------------------------------------
mod <- c("icar", "glm")
error_ha <- vector()
for (i in 1:length(mod)) {
	cat(glue("Projections in 2017: {mod[i]}"),"\n")
	if (!file.exists(glue("output/proj2017_{mod[i]}.tif"))) {
		deforest <- far$deforest(input_raster=glue("output/prob_{mod[i]}.tif"),
								 hectares=874211, # defor_area_2010_2017,
								 output_file=glue("output/proj2017_{mod[i]}.tif"),
								 blk_rows=128L)
		error_ha[i] <- deforest[[3]]
	}
}
sink("output/error_ha.txt")
print(error_ha)
sink()

# ================================
# Accuracy of projections in 2017
# ================================

## ----at 30m resolution----------------------------------------------------
mod <- c("icar", "glm")
index_diff2017_30m <- NULL
for (i in 1:length(mod)) {
	cat(glue("Differences in 2017: {mod[i]}_obs"),"\n")
	if (!file.exists(glue("output/diff2017_{mod[i]}_obs.tif"))) {
		far$r_diffproj(inputA=glue("output/proj2017_{mod[i]}.tif"),
					   inputB=glue("data/validation/proj2017_obs.tif"),
					   output_file=glue("output/diff2017_{mod[i]}_obs.tif"),
					   blk_rows=64L)
	}
	confmat <- far$mat_diffproj(glue("output/diff2017_{mod[i]}_obs.tif"))
	acc <- as.data.frame(far$accuracy(confmat)) %>% 
		dplyr::select(OA, EA, FOM, Spe, Sen, TSS, K)
	index_diff2017_30m <- rbind(index_diff2017_30m, acc)
}
acc_mod_30m <- index_diff2017_30m %>%
	dplyr::mutate(res=rep(1,2), mod=mod) %>%
	dplyr::select(8,9,1:7) %>%
	tidyr::gather(OA, EA, FOM, Spe, Sen, TSS, K, key="index",value="value")
write_csv(acc_mod_30m, "output/models_accuracy_30m.csv")

## ----resample_sum--------------------------------------------------------
dirout <- "output/accuracy"
if (!dir.exists(dirout)) {dir.create(dirout)}
resolutions <- c(2, 5, 10, 25, 50, 100, 250, 500, 1000)
models <- c("obs","icar","glm")
# Resample
for (i in 1:length(resolutions)) {
	res <- resolutions[i]
	for (j in 1:length(models)) {
		mod <- models[j]
		r_out <- glue("{dirout}/proj2017_{mod}_sum1_w{res}.tif")
		if (!file.exists(r_out)) {
			far$resample_sum(glue("output/proj2017_{mod}.tif"),
							 glue("{dirout}/proj2017_{mod}_sum1_w{res}.tif"),
							 val=1, window_size=as.integer(res))
			far$resample_sum(glue("output/proj2017_{mod}.tif"),
							 glue("{dirout}/proj2017_{mod}_sum0_w{res}.tif"),
							 val=0, window_size=as.integer(res))
		}
	}
}

## ----confusion_matrix----------------------------------------------------
confmat_icar <- confmat_glm <- list()
# Compute confusion matrix
for (i in 1:length(resolutions)) {
	res <- resolutions[i]
	# icar
	confmat_icar[[i]] <- far$confmat(
		r_obs0=glue("{dirout}/proj2017_obs_sum0_w{res}.tif"),
		r_obs1=glue("{dirout}/proj2017_obs_sum1_w{res}.tif"),
		r_pred0=glue("{dirout}/proj2017_icar_sum0_w{res}.tif"),
		r_pred1=glue("{dirout}/proj2017_icar_sum1_w{res}.tif"))
	# glm
	confmat_glm[[i]] <- far$confmat(
		r_obs0=glue("{dirout}/proj2017_obs_sum0_w{res}.tif"),
		r_obs1=glue("{dirout}/proj2017_obs_sum1_w{res}.tif"),
		r_pred0=glue("{dirout}/proj2017_glm_sum0_w{res}.tif"),
		r_pred1=glue("{dirout}/proj2017_glm_sum1_w{res}.tif"))
}

# Check number of cells
lapply(confmat_icar, sum)
lapply(confmat_glm, sum)

# Save lists
save(confmat_icar, confmat_glm, file="output/confmat.rda")

## ----accuracy------------------------------------------------------------
load("output/confmat.rda")
# Compute indices
acc_mod <- NULL
for (i in 1:length(resolutions)) {
	acc <- as.data.frame(far$accuracy(confmat_glm[[i]])) %>%
		dplyr::select(OA, EA, FOM, Spe, Sen, TSS, K)
	acc_mod <- rbind(acc_mod, acc) 
	acc <- as.data.frame(far$accuracy(confmat_icar[[i]])) %>%
		dplyr::select(OA, EA, FOM, Spe, Sen, TSS, K)
	acc_mod <- rbind(acc_mod, acc) 
}
acc_mod$mod <- rep(c("glm","icar"), 9)
acc_mod$res <- rep(resolutions, each=2)
acc_mod <- acc_mod %>%
	tidyr::gather(OA, EA, FOM, Spe, Sen, TSS, K, key="index",value="value")

# Export
write_csv(acc_mod, "output/models_accuracy.csv")

# Add results at 30m resolution
acc_mod_30m <- read_csv("output/models_accuracy_30m.csv")

# Plot with ggplot2
acc_mod <- read_csv("output/models_accuracy.csv") %>%
	dplyr::bind_rows(acc_mod_30m) %>%
	dplyr::filter(index %in% c("FOM", "OA")) %>%
	dplyr::mutate(res_km=30*res/1000)
p <- ggplot(acc_mod, aes(x=res_km, y=value, color=mod)) +
	geom_line() + xlim(0,15) + facet_wrap(.~index, ncol=4, scales="free")
ggsave("output/models_accuracy.png", p)

# =========================================================
# Model projection difference on the long term (until 2050)
# =========================================================

## ----proj2050, results="hide"--------------------------------------------
mod <- c("icar", "glm")
for (i in 1:length(mod)) {
	cat(glue("Projections in 2050: {mod[i]}"),"\n")
	if (!file.exists(glue("output/proj2050_{mod[i]}.tif"))) {
		deforest <- far$deforest(input_raster=glue("output/prob_{mod[i]}.tif"),
								 hectares=3400000,
								 output_file=glue("output/proj2050_{mod[i]}.tif"),
								 blk_rows=128L)
	}
}

## ----compute_diff2050----------------------------------------------------
comp <- c("icar_glm")
index_diff2050_30m <- NULL
for (i in 1:1) { #length(comp)) {
  cat(glue("Differences in 2050: {comp[i]}"),"\n")
  if (!file.exists(glue("output/diff2050_{comp[i]}.tif"))) {
    m1 <- unlist(strsplit(comp[i],"_"))[1]
    m2 <- unlist(strsplit(comp[i],"_"))[2]
    far$r_diffproj(inputA=glue("output/proj2050_{m1}.tif"),
                   inputB=glue("output/proj2050_{m2}.tif"),
                   output_file=glue("output/diff2050_{comp[i]}.tif"),
                   blk_rows=128L)
  }
  mat <- far$mat_diffproj(glue("output/diff2050_{comp[i]}.tif"))
  index_diff2050_30m <- rbind(index_diff2050_30m, accuracy_indices(mat))
}
diff2050_comp_30m <- index_diff2050_30m %>%
	dplyr::mutate(res=30,comp=comp[1]) %>%
	dplyr::select(comp, res, OA, EA, FOM, Spe, Sen, TSS, K)

## ----diff2050_plot-------------------------------------------------------
for (i in 1:1) { # length(comp)) {
  if (!file.exists(glue("output/diff2050_{comp[i]}.png"))) {
    cat(glue("Plot difference: {comp[i]}"),"\n")
    far$plot$differences(glue("output/diff2050_{comp[i]}.tif"),
                         borders="data/mada/mada38s.shp",
                         output_file=glue("output/diff2050_{comp[i]}.png"))
  }
}

## ----diff2050_icar_glm_ggplot--------------------------------------------
# Plot with ggplot2
r_diff <- resamp2df(input_file="output/diff2050_icar_glm.tif",
                    output_file="output/diff2050_icar_glm_1km.tif",
                    res=1000)
mada <- rgdal::readOGR(dsn="data/mada",layer="mada38s")
rect_df <- data.frame(xmin=c(346000,793000), xmax=c(439000,886000),
                      ymin=c(7387000,7815000), ymax=c(7480000,7908000),
                      id=c(1,2))
p <- diff_plot(input_df=r_diff,
               input_vector=mada,
               output_file="output/diff2050_icar_glm_ggplot.png",
               rect=rect_df)

# =======================================================================
# Using iCAR model to forecast deforestation on 2010--2055 and 2010--2085
# =======================================================================

deforest2055 <- far$deforest(input_raster=glue("output/prob_icar.tif"),
						 hectares=4500000,
						 output_file=glue("output/proj2055_icar.tif"),
						 blk_rows=128L)

deforest2085 <- far$deforest(input_raster=glue("output/prob_icar.tif"),
							 hectares=7500000,
							 output_file=glue("output/proj2085_icar.tif"),
							 blk_rows=128L)

# ========================================================
# End
# ========================================================