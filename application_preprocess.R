#### This script is for processing chemical and biomarker data at visit 1 ####
## load necessary packages
library(haven)
library(sas7bdat)
#### ------------------ Part 1: chemical data imputation ---------------------
demo <- read_sas("~/DUA/projects/FGS/data/formsingle.sas7bdat",
                  col_select = c("SUBJECT_ID","BMIGRP_fm002"))
names(demo) = c("ID", "BMIGRP")
demo$ID = as.integer(demo$ID)
chem <- read.sas7bdat("../data/chem_gdm_base_12022021.sas7bdat")
chemical_info <- read.csv("../data/chem_LOD.csv")
PFAS <- chem[,c(chemical_info$Chemicals[chemical_info$chemClass_short == "PFA"],
                "ID", "GDM")]
apply(PFAS,2,function(x) sum(is.na(x)))

# Remove PFBS for too many missing values
PFAS = PFAS[,-12]
nPFAS = ncol(PFAS) - 2
PFAS.missing = apply(PFAS,1,function(x) all(is.na(x[1:nPFAS])))
sum(PFAS.missing)
# Remove 41 samples who have no observations for any PFAS
PFAS = PFAS[!PFAS.missing,]

# Merge obesity information
PFAS = dplyr::inner_join(PFAS,demo, by = "ID")

# Match in LOD information

LODmat = matrix(NA,nrow(PFAS), nPFAS)
for(k in 1:nPFAS){
  infoK = chemical_info[chemical_info$Chemicals == names(PFAS)[k],]
  LODmat[PFAS$BMIGRP == 1, k] = infoK$LOD_nonobese
  LODmat[PFAS$BMIGRP == 2, k] = infoK$LOD_obese
}
LOD_percent = (colSums(PFAS[,1:nPFAS] < LODmat,na.rm = T) + apply(PFAS[,1:nPFAS],2, function(x) sum(is.na(x)))) / 
  nrow(PFAS)

# check why some of NMeFOSAA is NA for some samples
PFAS[is.na(PFAS$NMeFOSAA),]
# obese women have all NMeFOSAA NA. Consider remove this variable, as it is hard to impute. 
# NA is not caused by non-detectability. Probably the lab for obese women could not assay it.
PFAS_remove = names(LOD_percent)[LOD_percent >0.7]
PFAS_remove = c(PFAS_remove, "NMeFOSAA")
LODmat = LODmat[,-which(names(PFAS)%in%PFAS_remove)]
PFAS = PFAS[,-which(names(PFAS)%in%PFAS_remove)]
PFAS_comp = PFAS
for(i in 1:(ncol(PFAS)-3)){
  PFASi = PFAS[,i]
  PFASi[PFASi < LODmat[,i]] <- LODmat[PFASi < LODmat[,i], i]/sqrt(2)
  PFAS_comp[,i] <- PFASi
}
write.csv(PFAS_comp,file = "../data/working data/PFAS_baseline.csv",row.names = F)
#### ------------------ Part 2: thyroid function biomarkers at baseline ---------------------
thyroid = read_sas("../data/visit0.sas7bdat")

thyroid = thyroid[thyroid$ID %in% PFAS_comp$ID,]
thyroid = thyroid[,c("ID", "Age", "ENRRACE_fm002", "PreBMI", "Education", "nulliparity",
                     "GestationalAge","TSH","fT4","fT3")]
thyroid$R_fT4_fT3 = thyroid$fT4/thyroid$fT3

#### ------------------ Part 3: other chemical covariate ------------ ####
PFAS_comp = read.csv("../data/working data/PFAS_baseline.csv")
cov = chem[,c("ID","TotalLipid")]
data <- dplyr::inner_join(PFAS_comp,thyroid, by = "ID")
data <- dplyr::inner_join(data,cov,by = "ID")
data$logfT3 = log(data$fT3)
data$logfT4 = log(data$fT4)
data$logTSH = log(data$TSH)
data$logR_fT4_fT3 = log(data$R_fT4_fT3)
data_control = data[data$GDM == 0,]
apply(data_control,2, function(x) sum(is.na(x)))
data_GDM = data[data$GDM==1,]

write.csv(data,file = "../data/working data/data.csv", row.names = F)

#### ------------------ Part 4: adjust for covariates -----------------####
# Logarithm transformation on the chemicals (and biomarkers), 
# adjust for variates including Age,
# GDM, obesity, race, education, nulliparity and GA.
data <- read.csv("../data/working data/data.csv")
## fT3
lm.logfT3 = lm(logfT3~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                 BMIGRP + nulliparity + GestationalAge, data = data)
lm.fT3 = lm(fT3~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
              BMIGRP + nulliparity + GestationalAge, data = data)
summary(lm.logfT3)
summary(lm.fT3)
logfT3.m = data$logfT3
logfT3.m[!is.na(logfT3.m)] <-lm.logfT3$residuals

data.logfT3 <- data.frame(log(data[,1:8]),
                          Y = logfT3.m)

fT3.m = data$fT3
fT3.m[!is.na(fT3.m)] <-lm.fT3$residuals

data.fT3 <- data.frame(log(data[,1:8]),
                       Y = fT3.m)
head(data.fT3)
## fT4
lm.logfT4 = lm(logfT4~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                 BMIGRP + nulliparity + GestationalAge, data = data)
lm.fT4 = lm(fT4~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                 BMIGRP + nulliparity + GestationalAge, data = data)

summary(lm.logfT4)
summary(lm.fT4)
logfT4.m = data$logfT4
logfT4.m[!is.na(logfT4.m)] <-lm.logfT4$residuals

data.logfT4 <- data.frame(log(data[,1:8]),
                          Y = logfT4.m)

fT4.m = data$fT4
fT4.m[!is.na(fT4.m)] <-lm.fT4$residuals

data.fT4 <- data.frame(log(data[,1:8]),
                          Y = fT4.m)

## TSH
lm.logTSH = lm(logTSH~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                 BMIGRP + nulliparity + GestationalAge, data = data)
summary(lm.logTSH)
logTSH.m = data$logTSH
logTSH.m[!is.na(logTSH.m)] <-lm.logTSH$residuals

data.logTSH <- data.frame(log(data[,1:8]),
                          Y = logTSH.m)

lm.TSH = lm(TSH~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                 BMIGRP + nulliparity + GestationalAge, data = data)
summary(lm.TSH)
TSH.m = data$TSH
TSH.m[!is.na(TSH.m)] <-lm.TSH$residuals

data.TSH <- data.frame(log(data[,1:8]),
                          Y = TSH.m)

## R_fT4/fT3
lm.logR_fT4_fT3 = lm(logR_fT4_fT3~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
              BMIGRP + nulliparity + GestationalAge, data = data)
summary(lm.logR_fT4_fT3)
logR_fT4_fT3.m = data$logR_fT4_fT3
logR_fT4_fT3.m[!is.na(logR_fT4_fT3.m)] <-lm.logR_fT4_fT3$residuals

data.logR_fT4_fT3 <- data.frame(log(data[,1:8]),
                          Y = logR_fT4_fT3.m)
lm.R_fT4_fT3 = lm(R_fT4_fT3~ Age + GDM + as.factor(ENRRACE_fm002) + as.factor(Education)+
                       BMIGRP + nulliparity + GestationalAge, data = data)
summary(lm.R_fT4_fT3)
R_fT4_fT3.m = data$R_fT4_fT3
R_fT4_fT3.m[!is.na(R_fT4_fT3.m)] <-lm.R_fT4_fT3$residuals

data.R_fT4_fT3 <- data.frame(log(data[,1:8]),
                             Y = R_fT4_fT3.m)

#### ------------------- Part 5: rescale X --------------------------- ####
rescaleX = function(x){
  minX = apply(x,2,min)
  maxX = apply(x,2,max)
  N = nrow(x)
  
  return(list(x.scale = (x - matrix(rep(minX, each = N), nrow = N))/
              matrix(rep(maxX-minX, each = N), nrow = N),
         x.min = minX,
         x.max = maxX))
}
scale.info = vector(length = 8, "list")
i = 1
for(m in c("logfT3", "logfT4", "logTSH", "logR_fT4_fT3",
           "fT3", "fT4", "TSH", "R_fT4_fT3")){
  data = get(paste0("data.",m))
  x.scale = rescaleX(data[,-ncol(data)])
  data[,-ncol(data)] = x.scale$x.scale
  assign(paste0("data.",m),data)
  
  scale.info[[i]] <- rbind(x.scale$x.min, x.scale$x.max)
  i = i + 1
}
names(scale.info) <- c("logfT3", "logfT4", "logTSH", "logR_fT4_fT3",
                       "fT3", "fT4", "TSH", "R_fT4_fT3")
scale.info

save("scale.info", file = "scale_info.Rdata")
write.csv(data.logfT3, "../data/working data/logfT3.csv",row.names = F)
write.csv(data.fT3, "../data/working data/fT3.csv",row.names = F)
write.csv(data.logfT4, "../data/working data/logfT4.csv",row.names = F)
write.csv(data.fT4, "../data/working data/fT4.csv",row.names = F)
write.csv(data.logTSH, "../data/working data/logTSH.csv",row.names = F)
write.csv(data.TSH, "../data/working data/TSH.csv",row.names = F)
write.csv(data.logR_fT4_fT3, "../data/working data/logR_fT4_fT3.csv",row.names = F)
write.csv(data.R_fT4_fT3, "../data/working data/R_fT4_fT3.csv",row.names = F)


