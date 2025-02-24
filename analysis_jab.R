########################################################
################## IMMUNOLOGY ANALYSIS 
################### JAB SPECIAL ISSUE 
######################### 2024 
########################################################

# Load packages
library(dplyr)
library(tidyr)
library(car)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(blme)
library(parameters)
library(lattice)
library(performance)
library(emdi)
library(ggfortify)
library(ggplot2)
library(factoextra)
library(openxlsx)
library(ggplot2)
library(readxl)
library(readr)
library(writexl)
library(assignR)
library(MASS)
library(terra)
library(raster)
library(spatstat)
library(rnaturalearth)
library(sf)
library(sp)
library(vcd)

#### PART 1: MIGRATORY STATUS DETERMINATION ####

# Download isoscapes
isoMA <- getIsoscapes(isoType = "GlobalPrecipMA")
isodh.MA <- isoMA[[c("d2h_MA", "d2h_se_MA")]]


# Restrict spatial frame
europe <- ne_countries(continent = "europe", returnclass = "sf")
europe <- st_make_valid(europe)
europe_cropped <- st_crop(europe, xmin = -20, xmax = 45,
                          ymin = 30, ymax = 73)
spdf_europe <- as_Spatial(europe_cropped)


spdf_africa_continent <- ne_countries(continent = "africa")
spdf_africa2 <- subset(spdf_africa_continent, subregion == "Northern Africa" |
                         subregion == "Western Africa" |
                         subregion == "Middle Africa")
spdf_africa2 <- subset(spdf_africa2, !(sovereignt == "Angola"))
spdf_africa2 <- subset(spdf_africa2, !(sovereignt == "Democratic Republic of the Congo"))
spdf_africa3 <- ne_countries(country  = c("Kenya", "Burundi",
                                          "Eritrea", "Ethiopia",
                                          "Somalia", "South Sudan",
                                          "Rwanda", "Uganda",
                                          "Djibouti", "Somaliland"))
spdf_africa_red <- rbind(spdf_africa2, spdf_africa3)

mymask.red <- rbind(spdf_africa_red, spdf_europe)

# Download similar species data
control = subOrigData(marker = "d2H", taxon = c("Perdix perdix",
                                                "Passer domesticus",
                                                "Cyanistes caeruleus",
                                                "Dryocopus maritus",
                                                "Vanellus vanellus",
                                                "Turdus merula",
                                                "Columba palumbus",
                                                "Columba livia",
                                                "Turdus philomelos",
                                                "Locustella luscinioides",
                                                "Acrocephalus scirpaceus",
                                                "Phasianus colchicus",
                                                "Corvus corone",
                                                "Corvus monedula",
                                                "Turtle Dove",
                                                "Fringilla coelebs"),
                      mask = mymask.red)

control.data <- control$data

# Create our own control dataset with chick feathers
data_chicks <- read_delim("own_control_samples.csv", delim = ";", escape_double = FALSE,
                          col_types = cols(d2H = col_number(), 
                          d2H.sd = col_number(), lon = col_number(), 
                          lat = col_number()), trim_ws = TRUE)

atts <- data_chicks[, c(1:7)]
atts$Sample_ID <- as.integer(atts$Sample_ID)
lonlat <- data_chicks[, c(8,9)]

crsref <- "+proj=longlat +datum=WGS84"
pts <- vect(lonlat, crs=crsref)
values(pts) <- atts
control.chicks = pts

control.own <- rbind(control$data, control.chicks)

new_column_order = c("d2H", "d2H.sd", "Site_ID", "Site_name", "State", "Country", "Site_comme", "Dataset_ID", "Sample_ID", "Sample_ID_orig", "Taxon", "Group", "Source_quality", "Age_class", "Material_type", "Matrix", "d18O", "d18O.sd", "Sample_comments", "fieldnum", "nestcode", "Site")
control.own  <- control.own[, new_column_order]
control.own.df <- as.data.frame(control.own)


# Produce a calibrated, sample-type specific isoscape through correlating control samples
d2h_calibration3 = calRaster(known = control.own, isoscape = isodh.MA, mask = mymask.red)


# Import own data
dh2_samples <- read_delim("isotopes.csv", delim = ";", escape_double = FALSE,
                          col_types = cols(id = col_character(), 
                          d2H = col_number(), d2H.sd = col_number()), 
                          trim_ws = TRUE)

dh2_samples <- as.data.frame(dh2_samples)

# Transform to same reference
dh2_samples <- assignR::refTrans(dh2_samples, ref_scale = "OldEC.1_H_1")

# Run bayesian models for probability
dh2_prob = pdRaster(d2h_calibration3, dh2_samples, mask = mymask.red)

# Download probability dataframes
df.prob.red <- values(dh2_prob)
df.coords.red <- xyFromCell(dh2_prob, 1:ncell(dh2_prob))
df.prob.own.MA.red <- data.frame(df.coords.red, df.prob.red)
df.prob.own.MA.red <- na.omit(df.prob.own.MA.red)

write.csv(df.prob.own.MA.red, "df.prob.own.MA.csv")

# Quality analysis
qa = QA(control.own, isodh.MA, valiStation = 8, valiTime = 4, bySite = FALSE, by = 5, mask = mymask.red, name = "own_MA.red")
plot(qa)


# Clearer assignment maps
qtlRaster(dh2_prob3, threshold = 0.05) # (with 95% of security)

qtlRaster(dh2_prob3, threshold = 0.1) # (with 90% of security)





#### PART 2:  PHYSIOLOGICAL PARAMETERS ANALYSIS ####

# Import data
data <- read_excel("data_JAB.xlsx", sheet = "data")

# Calculate body condition index
#reg <- lm(mass ~ tar, data)
#cond <- residuals(reg)
#df.cond <- as.data.frame(cond)
#write.xlsx(df.cond, 'cond.xlsx')


# Prepare data for modelling (without unknown sex, age or behaviour)
data$age <- as.factor(data$age)
datajab <- data %>%
  mutate(sex = na_if(sex, "U"),
         age = na_if(age, "0"),
         behaviour = na_if(behaviour, "un"))


############################# CORRELATION MATRIX ###############################

# Subset dependent variables
physio <- data[, c(14:19, 22)]
physio$elisa_abs <- rowMeans(data[, c("elisa_abs1", "elisa_abs2")], na.rm = TRUE)
physio$h_l <- log(physio$h_l + 1)
physio$elisa_abs <- log(physio$elisa_abs + 1)
physio$agglutination <- log(physio$agglutination)

colnames(physio) <- c("fat", "muscle", "body condition", "H/L ratio",
                      "lysis", "agglutination", "haptoglobin", "IgY levels")
physio <- physio[, c("body condition", "muscle", "fat", "H/L ratio",
                     "IgY levels", "haptoglobin", "lysis", "agglutination")]

# Do correlation matrix
rcor <- rcorr(as.matrix(physio), type = c("pearson"))
coeff <- rcor$r
pvalues <- rcor$P
samples <- rcor$n

# Plot correlation matrix
corrplot(coeff, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.85,  cl.ratio = 0.30)
chart.Correlation(physio, histogram=TRUE, pch=60, number.cex = 2)


############################# BREEDING SEASON MODELS ###########################

# Create data subset
databreed <- datajab %>%
  filter(season == "breeding")

date.breed2 <- databreed$date^2

#### Adjust linear models for hypothesis

## Body condition

mod_cond <- lm(cond ~ date + date.breed2 + locality*behaviour + behaviour:sex + sex, databreed)
summary(mod_cond)
anova(mod_cond)

mod_cond2 <- lm(cond ~ date*locality + date:behaviour + behaviour*sex, databreed)
summary(mod_cond2)
anova(mod_cond2)

mod_cond3 <- lm(cond ~ date*locality + locality:behaviour + behaviour*sex, databreed)
summary(mod_cond3)
anova(mod_cond3)

mod_cond4 <- lm(cond ~ date*locality + locality:behaviour + behaviour+sex, databreed)
summary(mod_cond4)
anova(mod_cond4)

mod_cond5 <- lm(cond ~ date+locality + locality:behaviour + behaviour+sex, databreed)
summary(mod_cond5)
anova(mod_cond5)

mod_cond6 <- lm(cond ~ locality + locality:behaviour + behaviour+sex, databreed)
summary(mod_cond6)
anova(mod_cond6)


AIC(mod_cond)
AIC(mod_cond2)
AIC(mod_cond3)
AIC(mod_cond4)
AIC(mod_cond5)
AIC(mod_cond6)


## Muscle

muscle_log <- log(databreed$muscle)
mod_muscle <- lm(muscle_log ~ date + date.breed2 + locality*behaviour + behaviour:sex + sex + cond, databreed)
summary(mod_muscle)
anova(mod_muscle)

mod_muscle2 <- lm(muscle_log ~ date*behaviour + locality*behaviour + behaviour:sex + sex + cond, databreed)
summary(mod_muscle2)
anova(mod_muscle2)

mod_muscle3 <- lm(muscle_log ~ date*behaviour + locality + behaviour:sex + sex + cond, databreed)
summary(mod_muscle3)
anova(mod_muscle3)

mod_muscle4 <- lm(muscle_log ~ date*behaviour + locality + sex + cond, databreed)
summary(mod_muscle4)
anova(mod_muscle4)

mod_muscle5 <- lm(muscle_log ~ date*behaviour + sex + cond, databreed)
summary(mod_muscle5)
anova(mod_muscle5)

mod_muscle6 <- lm(muscle_log ~ date*behaviour + cond, databreed)
summary(mod_muscle6)
anova(mod_muscle6)


AIC(mod_muscle)
AIC(mod_muscle2)
AIC(mod_muscle3)
AIC(mod_muscle4)
AIC(mod_muscle5)
AIC(mod_muscle6)


## Fat

mod_fat <- glm(fat ~ date + date.breed2 + locality*behaviour + behaviour:sex + sex + cond, family = "poisson", databreed)
summary(mod_fat)
anova(mod_fat, test = "Chisq")

mod_fat2 <- glm(fat ~ date*behaviour + locality*behaviour + behaviour:sex + sex + cond, family = "poisson", databreed)
summary(mod_fat2)
anova(mod_fat2, test = "Chisq")

mod_fat3 <- glm(fat ~ date*behaviour + locality*behaviour + sex + cond, family = "poisson", databreed)
summary(mod_fat3)
anova(mod_fat3, test = "Chisq")

mod_fat4 <- glm(fat ~ date + locality*behaviour + sex + cond, family = "poisson", databreed)
summary(mod_fat4)
anova(mod_fat4, test = "Chisq")

mod_fat5 <- glm(fat ~ date + locality*behaviour + cond, family = "poisson", databreed)
summary(mod_fat5)
anova(mod_fat5, test = "Chisq")

mod_fat6 <- glm(fat ~ date + locality*behaviour, family = "poisson", databreed)
summary(mod_fat6)
anova(mod_fat6, test = "Chisq")

AIC(mod_fat)
AIC(mod_fat2)
AIC(mod_fat3)
AIC(mod_fat4)
AIC(mod_fat5)
AIC(mod_fat6)


## H/L ratio

mod_hl <- lm(log(h_l + 1) ~ date + date.breed2 + date:behaviour + locality*behaviour + behaviour:sex + sex + cond, databreed)
summary(mod_hl)
anova(mod_hl)

mod_hl2 <- lm(log(h_l + 1) ~ date*locality + date:behaviour + locality:behaviour + behaviour*sex + cond, databreed)
summary(mod_hl2)
anova(mod_hl2)

mod_hl3 <- lm(log(h_l + 1) ~ date + locality + date:behaviour + locality:behaviour + behaviour*sex + cond, databreed)
summary(mod_hl3)
anova(mod_hl3)

mod_hl4 <- lm(log(h_l + 1) ~ date + locality + date:behaviour + locality:behaviour + behaviour + sex + cond, databreed)
summary(mod_hl4)
anova(mod_hl4)

mod_hl5 <- lm(log(h_l + 1) ~ date + locality + date:behaviour + locality:behaviour + behaviour + sex, databreed)
summary(mod_hl5)
anova(mod_hl5)

mod_hl6 <- lm(log(h_l + 1) ~ date + locality + locality:behaviour + behaviour + sex, databreed)
summary(mod_hl6)
anova(mod_hl6)

mod_hl7 <- lm(log(h_l + 1) ~ date + locality + behaviour + sex, databreed)
summary(mod_hl7)
anova(mod_hl7)


AIC(mod_hl)
AIC(mod_hl2)
AIC(mod_hl3)
AIC(mod_hl4)
AIC(mod_hl5)
AIC(mod_hl6)
AIC(mod_hl7)


## Agglutination

# Create control correction to avoid models singularity
my.control <- lmerControl(check.conv.grad=.makeCC(action ="ignore", tol=1e-6, relTol=NULL),
                          optimizer="bobyqa", optCtrl=list(maxfun=100000))

# Fit linear mixed-effects models
mod_agg <- lmer(log(agglutination) ~ date + locality*behaviour + behaviour:sex + sex + cond + (1 | batch_hahl/plate_hahl),
                data = databreed,
                control=my.control)
summary(mod_agg)
anova(mod_agg)

mod_agg2 <- lmer(log(agglutination) ~ date + locality*behaviour + sex + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 control=my.control)
summary(mod_agg2)
anova(mod_agg2)

mod_agg3 <- lmer(log(agglutination) ~ date + locality + behaviour + sex + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 control=my.control)
summary(mod_agg3)
anova(mod_agg3)

mod_agg4 <- lmer(log(agglutination) ~ date + locality + sex + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 control=my.control)
summary(mod_agg4)
anova(mod_agg4)

mod_agg5 <- lmer(log(agglutination) ~ date + locality + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 control=my.control)
summary(mod_agg5)
anova(mod_agg5)

mod_agg6 <- lmer(log(agglutination) ~ locality + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 control=my.control)
summary(mod_agg6)
anova(mod_agg6)


AIC(mod_agg)
AIC(mod_agg2)
AIC(mod_agg3)
AIC(mod_agg4)
AIC(mod_agg5)
AIC(mod_agg6)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_agg6))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_agg6)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_agg6)
plotREsim(REsim(mod_agg6))

confint(mod_agg6, method = "boot")

# Check differences with non-mixed best fit model in same analysis
mod_agg.fixed <- lm(log(agglutination) ~ date + locality + cond, databreed)
summary(mod_agg.fixed)
anova(mod_agg.fixed)

AIC(mod_agg6)
AIC(mod_agg.fixed)

r2(mod_agg6)
r2(mod_agg.fixed)



## Lysis

# Fit linear mixed-effects models with bayesian priors to avoid converngence and singularity problems
mod_lys <- blmer(log(lysis + 1) ~ date + locality*behaviour + behaviour:sex + sex + cond + (1 | batch_hahl/plate_hahl),
                 data = databreed,
                 REML = T,
                 control = my.control,
                 cov.prior = invwishart)
summary(mod_lys)
anova(mod_lys)
parameters::p_value(mod_lys)

mod_lys2 <- blmer(log(lysis + 1) ~ date + locality + behaviour + behaviour:sex + sex + cond + (1 | batch_hahl/plate_hahl),
                  data = databreed,
                  REML = T,
                  control = my.control,
                  cov.prior = invwishart)
summary(mod_lys2)
anova(mod_lys2)
parameters::p_value(mod_lys2)

mod_lys3 <- blmer(log(lysis + 1) ~ locality + behaviour + behaviour:sex + sex + cond + (1 | batch_hahl/plate_hahl),
                  data = databreed,
                  REML = T,
                  control = my.control,
                  cov.prior = invwishart)
summary(mod_lys3)
anova(mod_lys3)
parameters::p_value(mod_lys3)

mod_lys4 <- blmer(log(lysis + 1) ~ locality + behaviour + behaviour:sex + sex + (1 | batch_hahl/plate_hahl),
                  data = databreed,
                  REML = T,
                  control = my.control,
                  cov.prior = invwishart)
summary(mod_lys4)
anova(mod_lys4)
parameters::p_value(mod_lys4)
df.residual(mod_lys4)

mod_lys5 <- blmer(log(lysis + 1) ~ behaviour + behaviour:sex + sex + (1 | batch_hahl/plate_hahl),
                  data = databreed,
                  REML = T,
                  control = my.control,
                  cov.prior = invwishart)
summary(mod_lys5)
anova(mod_lys5)
parameters::p_value(mod_lys5)


AIC(mod_lys)
AIC(mod_lys2)
AIC(mod_lys3)
AIC(mod_lys4)
AIC(mod_lys5)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_lys4))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps


# Check differences with non-mixed best fit model
mod_lys.fixed <- lm(lysis ~ behaviour + locality + behaviour:sex + sex, databreed)
summary(mod_lys.fixed)
anova(mod_lys.fixed)

AIC(mod_lys4)
AIC(mod_lys.fixed)

r2(mod_lys4)
r2(mod_lys.fixed)



## Immunoglobulins IgY (ELISA)

# Prepare data for duplicates analysis with samples ID
abs_elisa = c(datajab$elisa_abs1, datajab$elisa_abs2)

dataigy <- datajab %>%
  bind_rows(datajab)
dataigy <- cbind(dataigy, abs_elisa)

databreedigy <- dataigy %>%
  filter(season == "breeding")

# Fit linear mixed-effects models including ID and plate as random effects
mod_elisa <- lmer(formula = log(abs_elisa + 1) ~ date + locality*behaviour + behaviour:sex + sex + cond + (1 | fieldnum) + (1 | plate_igy),
                  data = databreedigy,
                  control=my.control)
summary(mod_elisa)
anova(mod_elisa)

mod_elisa2 <- lmer(formula = log(abs_elisa + 1) ~ date + locality + behaviour + behaviour:sex + sex + cond + (1 | fieldnum) + (1 | plate_igy),
                   data = databreedigy,
                   control=my.control)
summary(mod_elisa2)
anova(mod_elisa2)

mod_elisa3 <- lmer(formula = log(abs_elisa + 1) ~ date + locality + behaviour + sex + cond + (1 | fieldnum) + (1 | plate_igy),
                   data = databreedigy,
                   control=my.control)
summary(mod_elisa3)
anova(mod_elisa3)

mod_elisa4 <- lmer(formula = log(abs_elisa + 1) ~ date + behaviour + sex + cond + (1 | fieldnum) + (1 | plate_igy),
                   data = databreedigy,
                   control=my.control)
summary(mod_elisa4)
anova(mod_elisa4)

mod_elisa5 <- lmer(formula = log(abs_elisa + 1) ~ date + behaviour + cond + (1 | fieldnum) + (1 | plate_igy),
                   data = databreedigy,
                   control=my.control)
summary(mod_elisa5)
anova(mod_elisa5)

mod_elisa6 <- lmer(formula = log(abs_elisa + 1) ~ behaviour + cond + (1 | fieldnum) + (1 | plate_igy),
                   data = databreedigy,
                   control=my.control)
summary(mod_elisa6)
anova(mod_elisa6)


AIC(mod_elisa)
AIC(mod_elisa2)
AIC(mod_elisa3)
AIC(mod_elisa4)
AIC(mod_elisa5)
AIC(mod_elisa6)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_elisa6))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_elisa6)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_elisa6)
plotREsim(REsim(mod_elisa6))

confint(mod_elisa6, method = "boot")

# Check differences with non-mixed best fit model
mod_elisa.fixed <- lm(log(abs_elisa + 1) ~ date + behaviour + cond, databreed)
summary(mod_elisa.fixed)
anova(mod_elisa.fixed)

AIC(mod_elisa6)
AIC(mod_elisa.fixed)

r2(mod_elisa6)
r2(mod_elisa.fixed)



## Haptoglobin

# Fit linear mixed models including plate as random effect
mod_hp <- lmer(hapto_log ~ date + locality + behaviour*sex + cond + (1 | plate_hp),
               data = databreed,
               control = my.control)
summary(mod_hp)
anova(mod_hp)

mod_hp2 <- lmer(hapto_log ~ date + locality + behaviour + sex + cond + (1 | plate_hp),
                data = databreed,
                control = my.control)
summary(mod_hp2)
anova(mod_hp2)

mod_hp3 <- lmer(hapto_log ~ date + behaviour + sex + cond + (1 | plate_hp),
                data = databreed,
                control = my.control)
summary(mod_hp3)
anova(mod_hp3)

mod_hp4 <- lmer(hapto_log ~ date + behaviour + sex + (1 | plate_hp),
                data = databreed,
                control = my.control)
summary(mod_hp4)
anova(mod_hp4)

mod_hp5 <- lmer(hapto_log ~ date + behaviour + (1 | plate_hp),
                data = databreed,
                control = my.control)
summary(mod_hp5)
anova(mod_hp5)


AIC(mod_hp)
AIC(mod_hp2)
AIC(mod_hp3)
AIC(mod_hp4)
AIC(mod_hp5)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_hp5))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_hp5)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_hp5)
plotREsim(REsim(mod_hp5))

confint(mod_hp5, method = "boot")

# Check differences with non-mixed best fit model
mod_hp.fixed <- lm(hapto_log ~ date + behaviour + cond, databreed)
summary(mod_hp.fixed)
anova(mod_hp.fixed)

AIC(mod_hp5)
AIC(mod_hp.fixed)

r2(mod_hp5)
r2(mod_hp.fixed)



################# SEASONAL VARIATION WITHIN RESIDENTS ################

datars <- datajab %>%
  filter(behaviour == "rs")

## Body condition

mod_cond <- lm(cond ~ season + age*sex + season:age + locality:age + season:sex + locality:sex + cycle, datars)
summary(mod_cond)
anova(mod_cond)

mod_cond2 <- lm(cond ~ season + locality + age*sex + season:age + locality:age + season:sex + locality:sex + cycle, datars)
summary(mod_cond2)
anova(mod_cond2)

mod_cond3 <- lm(cond ~ season + locality + age*sex + season:age + locality:age + season:sex + cycle, datars)
summary(mod_cond3)
anova(mod_cond3)

mod_cond4 <- lm(cond ~ season + locality + age*sex + season:age + season:sex + cycle, datars)
summary(mod_cond4)
anova(mod_cond4)

mod_cond5 <- lm(cond ~ season + locality + age + sex + season:age + season:sex + cycle, datars)
summary(mod_cond5)
anova(mod_cond5)

mod_cond6 <- lm(cond ~ season + locality + age + sex + season:age + cycle, datars)
summary(mod_cond6)
anova(mod_cond6)

mod_cond7 <- lm(cond ~ season + locality + age + sex + cycle, datars) # BEST MODEL
summary(mod_cond7)
anova(mod_cond7)

mod_cond8 <- lm(cond ~ locality + age + sex + cycle, datars)
summary(mod_cond8)
anova(mod_cond8)


AIC(mod_cond)
AIC(mod_cond2)
AIC(mod_cond3)
AIC(mod_cond4)
AIC(mod_cond5)
AIC(mod_cond6)
AIC(mod_cond7)
AIC(mod_cond8)


## Muscle

mod_muscle <- lm(log(muscle) ~  season*locality + age*sex + season:age + locality:age + season:sex + locality:sex + cond + cycle, datars)
summary(mod_muscle)
anova(mod_muscle)

mod_muscle2 <- lm(log(muscle) ~  season*locality + age*sex + season:age + locality:age + season:sex + cond + cycle, datars)
summary(mod_muscle2)
anova(mod_muscle2)

mod_muscle3 <- lm(log(muscle) ~  season*locality + age*sex + season:age + locality:age + cond + cycle, datars)
summary(mod_muscle3)
anova(mod_muscle3)

mod_muscle4 <- lm(log(muscle) ~  season*locality + age*sex + season:age + cond + cycle, datars)
summary(mod_muscle4)
anova(mod_muscle4)

mod_muscle5 <- lm(log(muscle) ~  season + locality + age*sex + season:age + cond + cycle, datars)
summary(mod_muscle5)
anova(mod_muscle5)

mod_muscle6 <- lm(log(muscle) ~  season + locality + age*sex + cond + cycle, datars)
summary(mod_muscle6)
anova(mod_muscle6)

mod_muscle7 <- lm(log(muscle) ~  season + age*sex + cond + cycle, datars)
summary(mod_muscle7)
anova(mod_muscle7)

mod_muscle8 <- lm(log(muscle) ~  season + age + sex + cond + cycle, datars)
summary(mod_muscle8)
anova(mod_muscle8)

mod_muscle9 <- lm(log(muscle) ~  season + sex + cond + cycle, datars)
summary(mod_muscle9)
anova(mod_muscle9)

mod_muscle10 <- lm(log(muscle) ~  season + cond + cycle, datars)
summary(mod_muscle10)
anova(mod_muscle10)


AIC(mod_muscle)
AIC(mod_muscle2)
AIC(mod_muscle3)
AIC(mod_muscle4)
AIC(mod_muscle5)
AIC(mod_muscle6)
AIC(mod_muscle7)
AIC(mod_muscle8)
AIC(mod_muscle9)
AIC(mod_muscle10)


## Fat

mod_fat <- glm(fat ~  season*locality + age*sex + season:age + locality:age + season:sex + locality:sex + cycle, family = "poisson", datars)
summary(mod_fat)
anova(mod_fat, test = "Chisq")

mod_fat2 <- glm(fat ~  season*locality + age*sex + locality:age + season:sex + locality:sex + cycle, family = "poisson", datars)
summary(mod_fat2)
anova(mod_fat2, test = "Chisq")

mod_fat3 <- glm(fat ~  season*locality + age + sex + locality:age + season:sex + locality:sex + cycle, family = "poisson", datars)
summary(mod_fat3)
anova(mod_fat3, test = "Chisq")

mod_fat4 <- glm(fat ~  season*locality + age + sex + locality:age + locality:sex + cycle, family = "poisson", datars)
summary(mod_fat4)
anova(mod_fat4, test = "Chisq")

mod_fat5 <- glm(fat ~  season*locality + age + sex + locality:sex + cycle, family = "poisson", datars)
summary(mod_fat5)
anova(mod_fat5, test = "Chisq")

mod_fat6 <- glm(fat ~  season*locality + age + sex + cycle, family = "poisson", datars)
summary(mod_fat6)
anova(mod_fat6, test = "Chisq")

mod_fat7 <- glm(fat ~  season + locality + age + sex + cycle, family = "poisson", datars)
summary(mod_fat7)
anova(mod_fat7, test = "Chisq")

mod_fat8 <- glm(fat ~  season + locality + sex + cycle, family = "poisson", datars)
summary(mod_fat8)
anova(mod_fat8, test = "Chisq")

mod_fat9 <- glm(fat ~  season + locality + cycle, family = "poisson", datars)
summary(mod_fat9)
anova(mod_fat9, test = "Chisq")


AIC(mod_fat)
AIC(mod_fat2)
AIC(mod_fat3)
AIC(mod_fat4)
AIC(mod_fat5)
AIC(mod_fat6)
AIC(mod_fat7)
AIC(mod_fat8)
AIC(mod_fat9)


## H/L ratio

mod_hl <- lm(log(h_l + 1) ~ season*locality + age*sex + season:age + locality:age + season:sex + locality:sex + cycle, datars)
summary(mod_hl)
anova(mod_hl)

mod_hl2 <- lm(log(h_l + 1) ~ season*locality + age*sex + season:age + locality:age + season:sex + cycle, datars)
summary(mod_hl2)
anova(mod_hl2)

mod_hl3 <- lm(log(h_l + 1) ~ season + locality + age*sex + season:age + locality:age + season:sex + cycle, datars)
summary(mod_hl3)
anova(mod_hl3)

mod_hl4 <- lm(log(h_l + 1) ~ season + locality + age*sex + locality:age + season:sex + cycle, datars)
summary(mod_hl4)
anova(mod_hl4)

mod_hl5 <- lm(log(h_l + 1) ~ season + locality + age*sex + season:sex + cycle, datars)
summary(mod_hl5)
anova(mod_hl5)

mod_hl6 <- lm(log(h_l + 1) ~ season + locality + age*sex + cycle, datars)
summary(mod_hl6)
anova(mod_hl6)

mod_hl7 <- lm(log(h_l + 1) ~ season + locality + age + sex + cycle, datars)
summary(mod_hl7)
anova(mod_hl7)

mod_hl8 <- lm(log(h_l + 1) ~ season + age + sex + cycle, datars)
summary(mod_hl8)
anova(mod_hl8)


AIC(mod_hl)
AIC(mod_hl2)
AIC(mod_hl3)
AIC(mod_hl4)
AIC(mod_hl5)
AIC(mod_hl6)
AIC(mod_hl7)
AIC(mod_hl8)



## Agglutination

# Fit linear mixed-effects models including batch and plate as random factors
mod_agg <- lmer(log(agglutination) ~ season*locality + sex + season + locality + season:sex + locality:sex + cycle + (1 | batch_hahl/plate_hahl),
                data = datars,
                control = my.control)
summary(mod_agg)
anova(mod_agg)

mod_agg2 <- lmer(log(agglutination) ~ season*locality + sex + season + locality + season:sex + cycle + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg2)
anova(mod_agg2)

mod_agg3 <- lmer(log(agglutination) ~ season*locality + sex + season + locality + cycle + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg3)
anova(mod_agg3)

mod_agg4 <- lmer(log(agglutination) ~ season + locality + sex + season + locality + cycle + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg4)
anova(mod_agg4)

mod_agg5 <- lmer(log(agglutination) ~ locality + sex + locality + cycle + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg5)
anova(mod_agg5)

mod_agg6 <- lmer(log(agglutination) ~ sex + cycle + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg6)
anova(mod_agg6)

mod_agg7 <- lmer(log(agglutination) ~ sex + (1 | batch_hahl/plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_agg7)
anova(mod_agg7)


AIC(mod_agg)
AIC(mod_agg2)
AIC(mod_agg3)
AIC(mod_agg4)
AIC(mod_agg5)
AIC(mod_agg6)
AIC(mod_agg7)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_agg7))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_agg7)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_agg7)
plotREsim(REsim(mod_agg7))

confint(mod_agg7, method = "boot")


# Check differences with non-mixed best fit model
mod_agg.fixed <- lm(log(agglutination) ~ season*locality + age + cycle, datars)
summary(mod_agg.fixed)
anova(mod_agg.fixed)

AIC(mod_agg7)
AIC(mod_agg.fixed)

r2(mod_agg7)
r2(mod_agg.fixed)



## Lysis

# Fit linear mixed-effects models including plate as random factor
mod_lys <- lmer(log(lysis + 1) ~ season*locality + sex + season:sex + locality:sex + cycle + (1 | plate_hahl),
                data = datars,
                control = my.control)
summary(mod_lys)
anova(mod_lys)

mod_lys2 <- lmer(log(lysis + 1) ~ locality + sex*season + locality:sex + cycle + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_lys2)
anova(mod_lys2)

mod_lys3 <- lmer(log(lysis + 1) ~ locality + sex + season + locality:sex + cycle + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_lys3)
anova(mod_lys3)

mod_lys4 <- lmer(log(lysis + 1) ~ locality + sex + season + cycle + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_lys4)
anova(mod_lys4)

mod_lys5 <- lmer(log(lysis + 1) ~ locality + sex + season + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_lys5)
anova(mod_lys5)

mod_lys6 <- lmer(log(lysis + 1) ~ locality + sex + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_lys6)
anova(mod_lys6)

mod_lys7 <- lmer(log(lysis + 1) ~ locality + (1 | plate_hahl),
                 data = datars,
                 control = my.control)
summary(mod_hp7)
anova(mod_lys7)


AIC(mod_lys)
AIC(mod_lys2)
AIC(mod_lys3)
AIC(mod_lys4)
AIC(mod_lys5)
AIC(mod_lys6)
AIC(mod_lys7)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_lys7))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_lys7)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_lys7)
plotREsim(REsim(mod_lys7))

confint(mod_lys7, method = "boot")

# Check differences with non-mixed best fit model
mod_lys.fixed <- lm(log(lysis + 1) ~ season*locality + age + cycle, datars)
summary(mod_lys.fixed)
anova(mod_lys.fixed)

AIC(mod_lys7)
AIC(mod_lys.fixed)

r2(mod_lys7)
r2(mod_lys.fixed)



## Immunoglobulins IgY (ELISA)

# Fit linear mixed-effects models including ID and plate as random factors
datarsigy <- dataigy %>%
  filter(behaviour == "rs")

mod_elisa <- lmer(log(abs_elisa + 1) ~ season*locality + sex + season:sex + locality:sex + cycle + (1 | fieldnum) + (1 | plate_igy),
                  data = datarsigy,
                  control = my.control)
summary(mod_elisa)
anova(mod_elisa)

mod_elisa2 <- lmer(log(abs_elisa + 1) ~ season*locality + sex + season:sex + cycle + (1 | fieldnum) + (1 | plate_igy),
                   data = datarsigy,
                   control = my.control)
summary(mod_elisa2)
anova(mod_elisa2)

mod_elisa3 <- lmer(log(abs_elisa + 1) ~ season + locality + sex + season:sex + cycle + (1 | fieldnum) + (1 | plate_igy),
                   data = datarsigy,
                   control = my.control)
summary(mod_elisa3)
anova(mod_elisa3)

mod_elisa4 <- lmer(log(abs_elisa + 1) ~ season + locality + sex + cycle + (1 | fieldnum) + (1 | plate_igy),
                   data = datarsigy,
                   control = my.control)
summary(mod_elisa4)
anova(mod_elisa4)

mod_elisa5 <- lmer(log(abs_elisa + 1) ~ season + sex + cycle + (1 | fieldnum) + (1 | plate_igy),
                   data = datarsigy,
                   control = my.control)
summary(mod_elisa5)
anova(mod_elisa5)

mod_elisa6 <- lmer(log(abs_elisa + 1) ~ season + cycle + (1 | fieldnum) + (1 | plate_igy),
                   data = datarsigy,
                   control=my.control)
summary(mod_elisa6)
anova(mod_elisa6)


AIC(mod_elisa)
AIC(mod_elisa2)
AIC(mod_elisa3)
AIC(mod_elisa4)
AIC(mod_elisa5)
AIC(mod_elisa6)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_elisa6))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_elisa6)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_elisa6)
plotREsim(REsim(mod_elisa6))

confint(mod_elisa6, method = "boot")

# Check differences with non-mixed best fit model
mod_elisa.fixed <- lm(log(elisa_abs + 1) ~ season + locality + sex + locality:sex + cycle, datars)
summary(mod_elisa.fixed)
anova(mod_elisa.fixed)

AIC(mod_elisa6)
AIC(mod_elisa.fixed)

r2(mod_elisa6)
r2(mod_elisa.fixed)



## Haptoglobin

# Fit linear mixed-effects models including plate as random factor
mod_hp <- lmer(hapto_log ~ season*locality + sex + season:sex + locality:sex + cycle + (1 | plate_hp),
               data = datars,
               control = my.control)
summary(mod_hp)
anova(mod_hp)

mod_hp2 <- lmer(hapto_log ~ season*locality + sex + season:sex + cycle + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp2)
anova(mod_hp2)

mod_hp3 <- lmer(hapto_log ~ season + locality + sex + season:sex + cycle + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp3)
anova(mod_hp3)

mod_hp4 <- lmer(hapto_log ~ season + locality + sex + cycle + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp4)
anova(mod_hp4)

mod_hp5 <- lmer(hapto_log ~ season + locality + cycle + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp5)
anova(mod_hp5)

mod_hp6 <- lmer(hapto_log ~ season + cycle + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp6)
anova(mod_hp6)

mod_hp7 <- lmer(hapto_log ~ season + (1 | plate_hp),
                data = datars,
                control = my.control)
summary(mod_hp7)
anova(mod_hp7)


AIC(mod_hp)
AIC(mod_hp2)
AIC(mod_hp3)
AIC(mod_hp4)
AIC(mod_hp5)
AIC(mod_hp6)
AIC(mod_hp7)


# Extract random effects proportion of explained variance
var_comps <- as.data.frame(VarCorr(mod_hp7))
var_comps

var_total <- sum(var_comps$vcov)

var_comps$proportion <- var_comps$vcov / var_total
var_comps

ranef(mod_hp7)

# Explore residuals of the model and significance bootstrap intervals
check_model(mod_hp7)
plotREsim(REsim(mod_hp7))

confint(mod_hp7, method = "boot")

# Check differences with non-mixed best fit model
mod_hp.fixed <- lm(hapto_log ~ locality*season + cycle, datars)
summary(mod_hp.fixed)
anova(mod_hp.fixed)

AIC(mod_hp7)
AIC(mod_hp.fixed)

r2(mod_hp7)
r2(mod_hp.fixed)





######################## MEANS COMPARISON TABLE ########################

# Prepare IgY and haptoglobin values for means calculation
databehav_breed$elisa_abs <- rowMeans(databehav_breed[, c("elisa_abs1", "elisa_abs2")], na.rm = TRUE)
datars$elisa_abs <- rowMeans(datars[, c("elisa_abs1", "elisa_abs2")], na.rm = TRUE)

databehav_breed$hapto <- 10^databehav_breed$hapto_log
datars$hapto <- 10^datars$hapto_log

# Means for each variable in different migratory behaviours
behav_cond_means <- aggregate(cond ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_muscle_means <- aggregate(muscle ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_fat_means <- aggregate(fat ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_hl_means <- aggregate(h_l ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_lysis_means <- aggregate(lysis ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_agg_means <- aggregate(agglutination ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_elisa_means <- aggregate(elisa_abs ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)
behav_hp_means <- aggregate(hapto ~ behaviour, data = databehav_breed, FUN = mean, na.rm = TRUE)

behav_cond_sd <- aggregate(cond ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_muscle_sd <- aggregate(muscle ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_fat_sd <- aggregate(fat ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_hl_sd <- aggregate(h_l ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_lysis_sd <- aggregate(lysis ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_agg_sd <- aggregate(agglutination ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_elisa_sd <- aggregate(elisa_abs ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)
behav_hp_sd <- aggregate(hapto ~ behaviour, data = databehav_breed, FUN = sd, na.rm = TRUE)

behav_cond_counts <- aggregate(cond ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_muscle_counts <- aggregate(muscle ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_fat_counts <- aggregate(fat ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_hl_counts <- aggregate(h_l ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_lysis_counts <- aggregate(lysis ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_agg_counts <- aggregate(agglutination ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_elisa_counts <- aggregate(elisa_abs ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))
behav_hp_counts <- aggregate(hapto ~ behaviour, data = databehav_breed, FUN = function(x) length(na.omit(x)))


# Means for each variable in two seasons
season2_cond_means <- aggregate(cond ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_muscle_means <- aggregate(muscle ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_fat_means <- aggregate(fat ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_hl_means <- aggregate(h_l ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_lysis_means <- aggregate(lysis ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_agg_means <- aggregate(agglutination ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_elisa_means <- aggregate(elisa_abs ~ season, data = datars, FUN = mean, na.rm = TRUE)
season2_hp_means <- aggregate(hapto ~ season, data = datars, FUN = mean, na.rm = TRUE)

season2_cond_sd <- aggregate(cond ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_muscle_sd <- aggregate(muscle ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_fat_sd <- aggregate(fat ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_hl_sd <- aggregate(h_l ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_lysis_sd <- aggregate(lysis ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_agg_sd <- aggregate(agglutination ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_elisa_sd <- aggregate(elisa_abs ~ season, data = datars, FUN = sd, na.rm = TRUE)
season2_hp_sd <- aggregate(hapto ~ season, data = datars, FUN = sd, na.rm = TRUE)

season2_cond_counts <- aggregate(cond ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_muscle_counts <- aggregate(muscle ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_fat_counts <- aggregate(fat ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_hl_counts <- aggregate(h_l ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_lysis_counts <- aggregate(lysis ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_agg_counts <- aggregate(agglutination ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_elisa_counts <- aggregate(elisa_abs ~ season, data = datars, FUN = function(x) length(na.omit(x)))
season2_hp_counts <- aggregate(hapto ~ season, data = datars, FUN = function(x) length(na.omit(x)))



###############################  PLOTS  #################################

# Create data subsets for plots
databehav <- datajab %>%
  filter(behaviour != "NA")

databehavsex <- databehav %>%
  filter(sex != "U")

databehav_breed <- databehav %>%
  filter(season == "breeding")

databehavsex_breed <- databehav_breed %>%
  filter(sex != "NA")

databehavsex_breed_lys <- databehavsex_breed %>%
  filter(lysis != "NA")

# Create colors palettes
colors_jazz <- c("#BCEE68", "#ADD8E6")
colors_jazz2 <- c("#274E36","#F28C28")
colors_jazz3 <- c("gold3","#9932CC")


## FIG 1

violin_elisa <- ggplot() +
  geom_violin(data = databehav_breed, mapping = aes(behaviour, log(elisa_abs + 1), fill= behaviour), trim=TRUE, scale="count", alpha = 0.6) +
  scale_fill_manual(name = "Behaviour", values = colors_jazz2, labels = c("Migrants", "Residents")) +
  scale_colour_manual(name = "Behaviour", values = "grey60", labels = c("Migrants", "Residents")) +
  stat_summary(data = databehav_breed, fun.data = "mean_cl_boot", aes(fill = behaviour, x = behaviour, y = log(elisa_abs + 1)), color = "black", size = 0.3,
               position = position_dodge(width = 0.9)) + 
  ylab("IgY levels log10(absorbance)") + xlab("") + 
  scale_x_discrete(labels = c("migrants\nn = 46", "residents\nn = 27")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12))
violin_elisa
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 1A.pdf", violin_elisa, width = 5, height = 5)


hp <- ggplot(databehav_breed, aes(x = date, y = hapto_log)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  labs(x = "Julian date", y = "[Haptoglobin] log10 (mg/ml)")
hp
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/review/fig 1B.pdf", hp, width = 5, height = 5)


### FIG 2 

databehav_long <- databehav %>%
  pivot_longer(cols = c(muscle, fat), names_to = "measurement", values_to = "value")

stats <- databehav_long %>%
  group_by(behaviour, measurement) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(
    error = qt(0.975, df = n - 1) * sd / sqrt(n),  # Error estándar del intervalo de confianza del 95%
    ymin = mean - error,
    ymax = mean + error
  )

stats <- stats %>%
  mutate(measurement = factor(measurement, levels = c("muscle", "fat")))

combi1 <- ggplot(stats, aes(x = behaviour, y = mean, fill = behaviour, color = behaviour)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), size = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = colors_jazz2) + 
  scale_color_manual(values = colors_jazz2) + 
  facet_wrap(~measurement) +
  ylab("Scores") + 
  xlab("") + 
  scale_x_discrete(labels = c("migrants\nn = 50", "residents\nn = 33", "migrants\nn = 50", "residents\nn = 32")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16, face = "bold")
  )
combi1
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 2A.pdf", combi1) # Saving 6.32 x 5.96 in image


muscledate <- ggplot(databehav_breed, aes(x = date, y = muscle, color = behaviour)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = colors_jazz2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  labs(x = "Julian date", y = "Muscle scores")
muscledate
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 2B.pdf", muscledate, width = 5, height = 5)


fatdate <- ggplot(databehav_breed, aes(x = date, y = fat, color = behaviour)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = colors_jazz2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  labs(x = "Julian date", y = "Fat scores")
fatdate
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/review/fig 2C.pdf", fatdate, width = 5, height = 5)


### FIG 3

violin_hl <- ggplot() +
  geom_violin(data = datars, mapping = aes(season, log(h_l), fill= season), trim=FALSE, scale="count", alpha = 0.8) +
  geom_boxplot(data = datars, mapping = aes(season, log(h_l)), width = 0.1, outlier.shape = 16, outlier.size = 2, outlier.colour = "red", alpha = 0.7,  
               fill = NA, color = NA) +
  scale_fill_manual(name = "Season", values = colors_jazz, labels = c("Breeding", "Winter")) +
  scale_colour_manual(name = "Season", values = "grey60", labels = c("Breeding", "Winter")) +
  stat_summary(data = datars, fun.data = "mean_cl_boot", aes(fill = season, x = season, y = log(h_l)), color = "black", size = 0.3,
               position = position_dodge(width = 0.9)) + 
  ylab("log10 H/L ratio") + xlab("") + 
  scale_x_discrete(labels = c("breeding\nn = 32", "winter\nn = 45")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12))
  violin_hl
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 3A.pdf", violin_hl, width = 5, height = 5)


violin_m <- ggplot() +
  geom_violin(data = datars, mapping = aes(season, muscle, fill= season), trim=FALSE, scale="count", alpha = 0.8) +
  geom_boxplot(data = datars, mapping = aes(season, muscle), width = 0.1, outlier.shape = 16, outlier.size = 2, outlier.colour = "red", alpha = 0.7,  
               fill = NA, color = NA) +
  scale_fill_manual(name = "Season", values = colors_jazz, labels = c("Breeding", "Winter")) +
  scale_colour_manual(name = "Season", values = "grey60", labels = c("Breeding", "Winter")) +
  stat_summary(data = datars, fun.data = "mean_cl_boot", aes(fill = season, x = season, y = muscle), color = "black", size = 0.3,
               position = position_dodge(width = 0.9)) + 
  ylab("Muscle score") + xlab("") + 
  scale_x_discrete(labels = c("breeding\nn = 33", "winter\nn = 49")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 12))
violin_m
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 3B.pdf", violin_m, width = 5, height = 5)





### FIG 4

stats4 <- databehavsex %>%
  group_by(behaviour, sex) %>%
  summarise(
    mean = mean(log(lysis+1), na.rm = TRUE),
    sd = sd(log(lysis+1), na.rm = TRUE),
    n = n()
  ) %>%
  mutate(
    error = qt(0.975, df = n - 1) * sd / sqrt(n),
    ymin = mean - error,
    ymax = mean + error
  )

stats4 <- stats4 %>%
  mutate(behaviour = factor(behaviour, 
                            levels = c("mg", "rs"),
                            labels = c("Migrants", "Residents")))

stats4 <- stats4 %>%
  mutate(sex = factor(sex, 
                      levels = c("M", "F"),
                      labels = c("males", "females")))

combi4 <- ggplot(stats4, aes(x = sex, y = mean, fill = sex, color = sex)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), size = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = colors_jazz3) + 
  scale_color_manual(values = colors_jazz3) + 
  facet_wrap(~behaviour) +
  ylab("Complement activity log10(lysis score)") + 
  xlab("") + 
  scale_x_discrete(labels = c("males\nn = 18", "females\nn = 7")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16, face = "bold")
  )
combi4
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 4A.pdf", combi4)


stats2 <- databehavsex %>%
  group_by(behaviour, sex) %>%
  summarise(
    mean = mean(log(h_l), na.rm = TRUE),
    sd = sd(log(h_l), na.rm = TRUE),
    n = n()
  ) %>%
  mutate(
    error = qt(0.975, df = n - 1) * sd / sqrt(n),
    ymin = mean - error,
    ymax = mean + error
  )

stats2 <- stats2 %>%
  mutate(behaviour = factor(behaviour, 
                            levels = c("mg", "rs"),
                            labels = c("Migrants", "Residents")))

stats2 <- stats2 %>%
  mutate(sex = factor(sex, 
                      levels = c("M", "F"),
                      labels = c("males", "females")))

combi2 <- ggplot(stats2, aes(x = sex, y = mean, fill = sex, color = sex)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), size = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = colors_jazz3) + 
  scale_color_manual(values = colors_jazz3) + 
  facet_wrap(~behaviour) +
  ylab("log10 H/L ratio") + 
  xlab("") + 
  scale_x_discrete(labels = c("males\nn = 32", "females\nn = 38")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16, face = "bold")
  )
combi2
ggsave("~/Desktop/Papers/hoopoe_physiologyJAB/figures/fig 4B.pdf", combi2)


