
#####################################
#              ALL DATASETS        ##
#####################################

setwd("C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/DBDA2Eprograms")
setwd("C:/Users/muntazir mehdi/OneDrive - RMIT University/PhD Topic Clustering/DBDA2Eprograms")

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.

# Load necessary libraries
library(rjags)
library(coda)
library(parallel)
library(runjags)
library(tictoc)
library(ggplot2)
library(GGally)
library(DAAG)
library(mlbench)
library(irr)
library(psych)
library(rrcov)
library(mclust)
library(MASS)
library(titanic)
library(rattle)
library(factoextra)
library(corrplot)
library(gridExtra)
library(clustertend)
library(hopkins)
library(clValid)
library(clusterSim)
library(dirichletprocess)
library(fossil)
library(kamila)
library(moments)
library(ade4)

############################################################
# Load the Elisa dataset
data("elisa",package ="isdals")
elisa 
elisa_data <- as.matrix(elisa[, -c(1)])
true_lab <- as.numeric(elisa[,1])

# Normalize the CO2 data
elisa_data <- scale(elisa_data)

# Plotting the true label clustering
p1 <- fviz_cluster(list(data = elisa_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "ELISA")+ 
  theme(legend.position = "none") 

############################################################
# Load the Dune dataset
data("dunedata", package ="ade4")
dunedata$envir 
Dune_data <- as.matrix(dunedata$envir[, -c(4,5)])
true_lab_Dune=as.numeric(dunedata$envir[,5])

# Normalize the Dune dataset
Dune_scaled_data <- scale(Dune_data)
Dune_scaled_data

# Plotting the true label clustering
p2 <- fviz_cluster(list(data = Dune_scaled_data, cluster = true_lab_Dune), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "DUNE")+ 
  theme(legend.position = "none") 

############################################################
# Load the Puromycin dataset
data("Puromycin")
Puromycin
Pu_data <- as.matrix(Puromycin[,c(-3)]) # Convert all columns to numeric
true_lab <- as.numeric(Puromycin$state)    # Use the 'class' column as the true labels

# Normalize the Puromycin dataset
Pu_data <- scale(Pu_data)

# Plotting the true label clustering
p3 <- fviz_cluster(list(data = Pu_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "PUROMYCIN")+ 
  theme(legend.position = "none") 
############################################################
# Load the Chazeb dataset
data("chazeb",package ="ade4")
meau1 <- data.frame(chazeb$tab,chazeb$cla)
names(meau1)[ncol(meau1)] <- "cla"
Chaz_data <- as.matrix(meau1[, -c(4,7)])
true_lab_Chaz=as.numeric(meau1[,7])

# Normalize the corvus dataset
Chaz_scaled_data <- scale(Chaz_data)

# Plotting the true label clustering
p4 <- fviz_cluster(list(data = Chaz_scaled_data, cluster = true_lab_Chaz), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "CHAZEB")+ 
  theme(legend.position = "none") 
############################################################
# Load the Meau dataset
data("meau", package ="ade4")
meau
meau1 <- data.frame(meau$env,meau$design$season)
names(meau1)[ncol(meau1)] <- "Season"
Meau_data <- as.matrix(meau1[, -c(7,8,11)])
true_lab_Meau=as.numeric(meau1[,11])

# Normalize the corvus dataset
Meau_scaled_data <- scale(Meau_data)

# Plotting the true label clustering
p5 <- fviz_cluster(list(data = Meau_scaled_data, cluster = true_lab_Meau), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "MEAU")+ 
  theme(legend.position = "none") 
############################################################
# Load the Corvus dataset
data("corvus", package ="ade4")
corvus
Cor_data <- as.matrix(corvus[, -c(3,4)])
true_lab=as.numeric(corvus[,3])

# Normalize the Corvus dataset
Cor_data <- scale(Cor_data)

# Plotting the true label clustering
p6 <- fviz_cluster(list(data = Cor_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "CORVUS")+ 
  theme(legend.position = "none") 
############################################################
# Load the Macroloire dataset
data("macroloire",package ="ade4")
macroloire 
Cor_data <- as.matrix(macroloire$envir[, c(2,3)])
true_lab=as.numeric(macroloire$envir[,5])

# Normalize the Corvus dataset
Cor_data <- scale(Cor_data)

# Plotting the true label clustering
p7 <- fviz_cluster(list(data = Cor_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "MACROLOIRE")+ 
  theme(legend.position = "none") 
############################################################
# Load the pap dataset
data("pap", package ="ade4")
pap
pap1 <- data.frame(pap$tab,pap$taxo$superfamille,pap$taxo$famille,pap$taxo$genre)
names(pap1)[5:7] <- c("superfamille", "famille", "genre")
Cor_data <- as.matrix(pap1[, -c(5:7)])
true_lab=as.numeric(pap1[,5])

# Normalize the Corvus dataset
Cor_data <- scale(Cor_data)

# Plotting the true label clustering
p8 <- fviz_cluster(list(data = Cor_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "PAP")+ 
  theme(legend.position = "none") 
############################################################
# Load the carniherbi dataset
data("carniherbi49",package ="ade4")
carniherbi49                            
help(carniherbi49)
str(carniherbi49)
uniques <- lapply(carniherbi49, unique)
uniques
Carni_data <- as.matrix(carniherbi49$tab2[, -c(1)])
true_lab_Carni=as.numeric(carniherbi49$tab2[,1])

# Normalize the Corvus dataset
Carni_scaled_data <- scale(Carni_data)

# Plotting the true label clustering
p9 <- fviz_cluster(list(data = Carni_scaled_data, cluster = true_lab_Carni), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "CARNIHERBI49")+ 
  theme(legend.position = "none") 
############################################################
# Load the oribatid dataset
data("oribatid", package ="ade4")
oribatid1 <- oribatid$env
Cor_data <- as.matrix(oribatid1[, -c(1:3)])
true_lab=as.numeric(oribatid1[,3])

# Normalize the Corvus dataset
Cor_data <- scale(Cor_data)

# Plotting the true label clustering
p10 <- fviz_cluster(list(data = Cor_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "ORBATID")+ 
  theme(legend.position = "none") 
############################################################
# Load the Aravo dataset
data(aravo,package = "ade4")
Son_data <- as.matrix(aravo$env[,-c(3,5)]) # Convert all columns to numeric
true_lab_Son <- as.numeric(aravo$env$ZoogD)    # Use the 'class' column as the true labels

# Normalize the Sonar data
Son_scaled_data <- scale(Son_data)

# Plotting the true label clustering
p11 <- fviz_cluster(list(data = Son_scaled_data, cluster = true_lab_Son), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "ARAVO")+ 
  theme(legend.position = "none") 

############################################################
# Load the Ecomor dataset
data(ecomor, package = "ade4")
ecomor                            
help(ecomor)
str(ecomor)
uniques <- lapply(ecomor, unique)
uniques
Ecomo_data <- as.matrix(ecomor$morpho)
true_lab_Ecomo=as.numeric(ecomor$taxo$Ordre)

# Normalize the Corvus dataset
Ecomo_scaled_data <- scale(Ecomo_data)

# Plotting the true label clustering
p12 <- fviz_cluster(list(data = Ecomo_scaled_data, cluster = true_lab_Ecomo), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "ECOMOR")+ 
  theme(legend.position = "none") 
#################################################################
# Load the Cats dataset
data("cats")
cats
cat_data <- as.matrix(cats[,-1]) # Convert all columns to numeric
true_lab <- as.numeric(cats$Sex)    # Use the 'class' column as the true labels

# Normalize the Cats dataset
cat_data <- scale(cat_data)

# Plotting the true label clustering
p13 <- fviz_cluster(list(data = cat_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "CAT")+ 
  theme(legend.position = "none") 
#################################################################
# Load the Iris dataset
data(iris)
iris
iris_data <- as.matrix(iris[, c(1,2,4)])
true_lab=iris[,5]
true_lab=as.numeric(true_lab)

# Normalize the iris data
iris_data <- scale(iris_data)

# Plotting the true label clustering
p14 <- fviz_cluster(list(data = iris_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "IRIS")+ 
  theme(legend.position = "none") 
#################################################################
# Load the Fish Catch dataset
data("fish")

#after dropping the highly correlated variables, the variables Length2, Height, and Width were
#used for the analysis 
Fish_data <- as.matrix(fish[,c(3,5,6)])
true_lab_Fish=fish[,7]

# Normalize the Fish Catch data
Fish_scaled_data <- scale(Fish_data)

# Plotting the true label clustering
p15 <- fviz_cluster(list(data = Fish_scaled_data, cluster = true_lab_Fish), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "FISH")+ 
  theme(legend.position = "none") 

#################################################################
# Load the Wine dataset
data("wine", package = "rattle")
W_data <- as.matrix(wine[,-1]) # Convert all columns to numeric
true_lab <- as.numeric(wine$Type)    # Use the 'class' column as the true labels

# Normalize the Wine dataset
W_data <- scale(W_data)

# Plotting the true label clustering
p16 <- fviz_cluster(list(data = W_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "WINE")+ 
  theme(legend.position = "none") 

############################################################
# Load the AIS dataset
data("ais")
ais
AIS_data <- as.matrix(ais[,c(1:2,4:6,8:10)])
true_lab <- as.numeric(ais$sex)

# Normalize the AIS data
AIS_data <- scale(AIS_data)

# Plotting the true label clustering
p17 <- fviz_cluster(list(data = AIS_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "AIS")+ theme(legend.position = "none") 

############################################################
# Load the dataset
data(lascaux,package = "ade4")
help(lascaux)
str(lascaux)
uniques <- lapply(lascaux, unique)
uniques
Lascaux_data <- as.matrix(lascaux$tap)
true_lab_Lascaux= as.numeric(lascaux$sex)

# Normalize the Corvus dataset
Lascaux_scaled_data <- scale(Lascaux_data)

# Plotting the true label clustering
p18 <- fviz_cluster(list(data = Lascaux_scaled_data, cluster = true_lab_Lascaux), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "LASCAUX")+ 
  theme(legend.position = "none") 
############################################################
# Load the Boston dataset
data("Boston")
Boston
Boston$chas[Boston$chas == 0] = "non-tract B.R"
Boston$chas[Boston$chas == 1] = "tract B.R"
Bos_data <- as.matrix(Boston[,-4])  # Exclude the 'class' column for clustering
true_lab <- as.numeric(factor(Boston[,4], levels = c("non-tract B.R", "tract B.R")))     # Use the 'class' column as the true labels

# Normalize the Boston data
Bos_data <- scale(Bos_data)

# Plotting the true label clustering
p19 <- fviz_cluster(list(data = Bos_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "BOSTON")+ 
  theme(legend.position = "none") 
############################################################
# Load the Chick Weight dataset
data("ChickWeight")
ChickWeight
CW_data1 <- sapply(ChickWeight, as.numeric) # Convert all columns to numeric
CW_data <- CW_data1[,-c(3,4)]
true_lab <- as.numeric(ChickWeight$Diet)    # Use the 'class' column as the true labels

# Normalize the Chick Weight dataset
CW_data <- scale(CW_data)

# Plotting the true label clustering
p20 <- fviz_cluster(list(data = CW_data, cluster = true_lab), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "CHICK WEIGHT")+ 
  theme(legend.position = "none") 
############################################################

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11,p12, ncol = 3)
grid.arrange(p13, p14, p15, p16, p17,p18,p19,p20, ncol = 3)

################################################################
