### Load libraries
library()
library()


### Load data
wd <- getwd()
patients <- read.csv(paste(wd, "/Patients_HCPAnywhere.csv", sep = ""), header = TRUE, sep = ",")
diagnostics <- read.csv(paste(wd, "/HCPAnywhere_Diagnostics.csv", sep = ""), header = TRUE, sep = ",")
plombage <- read.csv(paste(wd, "/HCPAnywhere_Plomb_Data.csv", sep = ""), header = TRUE, sep = ",")
