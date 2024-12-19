#################################################################
#                       Lascaux DATASET                      ##
#################################################################

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
library(mclust)
library(factoextra)
library(gridExtra)
library(clustertend)
library(hopkins)
library(clValid)
library(clusterSim)
library(dirichletprocess)
library(fossil)
library(corrplot)
library(moments)
library(MASS)
library(ade4)
library(aricode)

# Load the Corvus dataset
data(lascaux,package = "ade4")
help(lascaux)
str(lascaux)
uniques <- lapply(lascaux, unique)
uniques
Lascaux_data <- as.matrix(lascaux$tap)
true_lab_Lascaux= as.numeric(lascaux$sex)

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Lascaux_data)

# 2. Check for missing values
sum(is.na(Lascaux_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(Lascaux_data)) {
  hist(Lascaux_data[, i], main=colnames(Lascaux_data)[i], xlab=colnames(Lascaux_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Lascaux_data, main="Lascaux Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_Lascaux])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(Lascaux_data)
Cor_normalized$Season <- lascaux$sex

par(mfrow=c(1,2))  # Reset the plotting area for boxplots
for (i in 1:ncol(Lascaux_data)) {
  boxplot(Lascaux_data[, i] ~ Cor_normalized$Season, main=colnames(Lascaux_data)[i], xlab="Region", ylab=colnames(Lascaux_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Lascaux_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Corvus dataset
Lascaux_scaled_data <- scale(Lascaux_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Lascaux_scaled_data, 2, skewness)
kurtosis_values <- apply(Lascaux_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

par(mfrow=c(1, 2))
# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=1)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=1)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(Lascaux_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Lascaux_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Lascaux_scaled_data[,i])$out
  print(paste("Feature:", colnames(Lascaux_scaled_data)[i]))
  print(paste("Outliers:", length(outliers1[[i]])))
  total_outliers1 <- total_outliers1 + length(outliers1[[i]])
  outlier_rows1 <- unique(c(outlier_rows1, outliers1[[i]]))  # Collect unique rows with outliers
}

# Print the total number of outliers after the loop
print(paste("Total outliers:", total_outliers1))

# Print the number of rows with outliers
print(paste("Number of rows with outliers:", length(outlier_rows1)))

# Create the detect_outlier function
detect_outlier <- function(x) {
  
  # Calculate the first quantile
  Quantile1 <- quantile(x, probs = 0.25)
  
  # Calculate the third quantile
  Quantile3 <- quantile(x, probs = 0.75)
  
  # Calculate interquartile range (IQR)
  IQR <- Quantile3 - Quantile1
  
  # Return TRUE if an element is an outlier, FALSE otherwise
  outliers <- x > Quantile3 + (IQR * 1.5) | x < Quantile1 - (IQR * 1.5)
  return(outliers)
}

# Apply detect_outlier function to each column
outliers_matrix <- apply(Lascaux_scaled_data, 2, detect_outlier)
outliers_matrix

# Count the number of rows that have at least one TRUE (outlier)
rows_with_outliers <- rowSums(outliers_matrix) > 0
num_rows_with_outliers <- sum(rows_with_outliers)

# Print the result
print(paste("Number of rows with at least one outlier:", num_rows_with_outliers))

###############################################################################################################
#                                          DPGMVLG
###############################################################################################################

#Plotting true labels
fviz_cluster(list(data = Lascaux_scaled_data, cluster = true_lab_Lascaux), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Lascaux_scaled_data)
D <- ncol(Lascaux_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Lascaux_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Lascaux_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Lascaux_scaled_data))^(1 / (D - 1))

model_string <- "
Data {
  C <- 1000000000
  for (i in 1:P) {
     zeros[i] <- 0
  }
  v_g <- 1.4262       
  delta_g <- deltaData
}
model {
  for (i in 1:P) {
   for (g in 1:G){
      for ( j in 1:D){
        summ2[i,j,g] <- exp(mu_g[j, g] * x[i, j])/lambda_g[j, g]
      }
      logLik[i, g] <- z[i, g] * (log(pi_g[g])+ v_g * log(delta_g) + sum(log(mu_g[ ,g])) - v_g * sum(log(lambda_g[, g])) - loggam(v_g) ) + 
                      v_g * z[i, g] * sum(mu_g[1:D, g] * x[i, 1:D]) - z[i, g] * sum(summ2[i, , g])  +  
                      z[i, g] * (1 - delta_g) * prod(1/lambda_g[ , g]) * exp(-2.07728 * (D-1) + sum(mu_g[ , g] * x[i, ]) ) 
    
      z[i, g] <- ifelse(cluster[i] == g, 1, 0)
   }
   cluster[i] ~ dcat(pi_g) 
   zeros[i] ~ dpois(-sum(logLik[i,]) + C)
  }
   for (g in 1:G) {
    for (j in 1:D) {
       mu_g[j, g] ~ dgamma(4,2)          
       lambda_g[j, g] ~ dgamma(1,4)
    }
    alpha[g] <- 5
   }
   pi_g[1:G] ~ ddirch(alpha/sum(alpha))
}
"
# Write the model to a file
writeLines(model_string, con = "TEMPmodel.txt")

# Parameters to monitor
params <- c("delta_g","mu_g", "lambda_g", "cluster", "z")

inits <- list(list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 22021),
              list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 32019))

# Track the runtime for running the JAGS model
tic("JAGS Model Runtime")
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=params ,
                        data=data_list ,
                        n.chains=2 ,
                        adapt = 3000,
                        burnin=2000,
                        sample=1000, 
                        inits = inits,
                        thin = 2)
run_time_Lascaux <- toc()
codaSamples_Lascaux = as.mcmc.list( runJagsOut )
summaryChains_Lascaux <- summary(codaSamples_Lascaux)

diagMCMC( codaObject=codaSamples_Lascaux , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Lascaux , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Lascaux$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Lascaux$statistics),1], P, G)
z_mode_Lascaux <- apply(matrix(summaryChains_Lascaux$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Lascaux$statistics),1], P, G),1, which.max)
z_mode_Lascaux

plot(Lascaux_scaled_data, col= z_mode_Lascaux, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Lascaux), col=unique(z_mode_Lascaux), pch=16, title="Cluster")
table( true_lab_Lascaux , z_mode_Lascaux)

# To switch to the same labels in the true clustering
new <- c(1,2)
old <- c(2,1)
z_mode_Lascaux[z_mode_Lascaux %in% old] <- new[match(z_mode_Lascaux,old,nomatch = 0)]

plot(Lascaux_scaled_data, col= z_mode_Lascaux, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Lascaux), col=unique(z_mode_Lascaux), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Lascaux , z_mode_Lascaux)
kappa2(data.frame(rater1 = true_lab_Lascaux, rater2 = z_mode_Lascaux))

calculate_dp_Gmvlg_clustering_metrics <- function(dataLascaux_1, true_clusters_Lascaux) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Lascaux <- table(true_clusters_Lascaux, z_mode_Lascaux)
  kappa_result_Lascaux <- kappa2(data.frame(rater1 = true_clusters_Lascaux, rater2 = z_mode_Lascaux))
  ari_result_Lascaux <- adjustedRandIndex(true_clusters_Lascaux, z_mode_Lascaux)
  nmi_result_Lascaux <- NMI(true_clusters_Lascaux, z_mode_Lascaux)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Lascaux,
    Kappa = kappa_result_Lascaux$value,
    ARI = ari_result_Lascaux,
    NMI = nmi_result_Lascaux,
    CPU_RUN_TIME = run_time_Lascaux$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsLascaux_1 <- calculate_dp_Gmvlg_clustering_metrics(dataLascaux_1 = Lascaux_scaled_data, true_clusters_Lascaux = true_lab_Lascaux)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsLascaux_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsLascaux_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsLascaux_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsLascaux_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsLascaux_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsLascaux_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Lascaux <- eclust(Lascaux_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Lascaux <- toc()
kmeans_clusters_Lascaux <- kmeans_result_Lascaux$cluster
table(true_lab_Lascaux, kmeans_clusters_Lascaux)
new <- 1:2
old <- c(2,1)
kmeans_clusters_Lascaux[kmeans_clusters_Lascaux %in% old] <- new[match(kmeans_clusters_Lascaux, old, nomatch = 0)]
table(true_lab_Lascaux, kmeans_clusters_Lascaux)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Lascaux <- eclust(Lascaux_scaled_data, "clara", G, graph = FALSE)
clara_time_Lascaux <- toc()
clara_clusters_Lascaux <- clara_result_Lascaux$cluster
table(true_lab_Lascaux, clara_clusters_Lascaux)
#new <- 1:2
#old <- c(2,1)
#clara_clusters_Lascaux[clara_clusters_Lascaux %in% old] <- new[match(clara_clusters_Lascaux, old, nomatch = 0)]
table(true_lab_Lascaux, clara_clusters_Lascaux)

# PAM clustering
tic("PAM Runtime")
pam_result_Lascaux <- eclust(Lascaux_scaled_data, "pam", G, graph = FALSE)
pam_time_Lascaux <- toc()
pam_clusters_Lascaux <- pam_result_Lascaux$cluster
table(true_lab_Lascaux, pam_clusters_Lascaux)
#new <- 1:2
#old <- c(2,1)
#pam_clusters_Lascaux[pam_clusters_Lascaux %in% old] <- new[match(pam_clusters_Lascaux, old, nomatch = 0)]
table(true_lab_Lascaux, pam_clusters_Lascaux)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Lascaux <- hclust(dist(Lascaux_scaled_data), method = "ward.D2")
hclust_time_Lascaux <- toc()
hclust_clusters_Lascaux <- cutree(hclust_result_Lascaux, k = G)
table(true_lab_Lascaux, hclust_clusters_Lascaux)
#new <- 1:2
#old <- c(2,1)
#hclust_clusters_Lascaux[hclust_clusters_Lascaux %in% old] <- new[match(hclust_clusters_Lascaux, old, nomatch = 0)]
table(true_lab_Lascaux, hclust_clusters_Lascaux)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Lascaux <- Mclust(Lascaux_scaled_data, G = G)
mclust_time_Lascaux <- toc()
summary(mclust_result_Lascaux)
mclust_clusters_Lascaux <- mclust_result_Lascaux$classification
table(true_lab_Lascaux, mclust_clusters_Lascaux)
#new <- 1:2
#old <- c(2,1)
#mclust_clusters_Lascaux[mclust_clusters_Lascaux %in% old] <- new[match(mclust_clusters_Lascaux, old, nomatch = 0)]
table(true_lab_Lascaux, mclust_clusters_Lascaux)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Lascaux_scaled_data, alphaPriors = c(4,2))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Lascaux <- toc()
# Extract clusters 
dpMVN_clusters_Lascaux <- as.numeric(dp$clusterLabels)
new <- c(1,2)
old <- c(2,1)
dpMVN_clusters_Lascaux[dpMVN_clusters_Lascaux %in% old] <- new[match(dpMVN_clusters_Lascaux, old, nomatch = 0)]
print(dpMVN_clusters_Lascaux)
table(true_lab_Lascaux, dpMVN_clusters_Lascaux)

calculate_clustering_metricsLascaux_2 <- function(dataLascaux_2, true_clusters_Lascaux, estimated_clusters_Lascaux_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Lascaux_list)) {
    clusters <- estimated_clusters_Lascaux_list[[method_name]]
    
    # Calculate metrics
    table_result_Lascaux <- table(true_clusters_Lascaux, clusters)
    kappa_result_Lascaux <- kappa2(data.frame(rater1 = true_clusters_Lascaux, rater2 = clusters))
    ari_result_Lascaux <- adjustedRandIndex(true_clusters_Lascaux, clusters)
    nmi_result_Lascaux <- NMI(true_clusters_Lascaux, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Lascaux,
      Kappa = kappa_result_Lascaux$value,
      ARI = ari_result_Lascaux,
      NMI = nmi_result_Lascaux,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Lascaux <- list(
  KMeans = kmeans_clusters_Lascaux,
  CLARA = clara_clusters_Lascaux,
  PAM = pam_clusters_Lascaux,
  Hierarchical = hclust_clusters_Lascaux,
  Mclust = mclust_clusters_Lascaux,
  DPMVN = dpMVN_clusters_Lascaux,
  True = true_lab_Lascaux
)

times_list_Pios <- list(
  KMeans = kmeans_time_Lascaux,
  CLARA = clara_time_Lascaux,
  PAM = pam_time_Lascaux,
  Hierarchical = hclust_time_Lascaux,
  Mclust = mclust_time_Lascaux,
  DPMVN = DPMVN_time_Lascaux
)

# Call the function
clustering_metricsLascaux_2 <- calculate_clustering_metricsLascaux_2(dataLascaux_2 = Lascaux_scaled_data, true_clusters_Lascaux = true_lab_Lascaux, estimated_clusters_Lascaux_list = cluster_result_Lascaux, times_list = times_list_Pios)

# Print the results for each method
for (method_name in names(clustering_metricsLascaux_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsLascaux_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsLascaux_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsLascaux_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsLascaux_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsLascaux_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Lascaux_K <- data.frame(clustering_metricsLascaux_2$True$Kappa, dp_Gmvlg_metricsLascaux_1$dp_Gmvlg$Kappa, clustering_metricsLascaux_2$KMeans$Kappa, clustering_metricsLascaux_2$CLARA$Kappa, clustering_metricsLascaux_2$PAM$Kappa, clustering_metricsLascaux_2$Hierarchical$Kappa, clustering_metricsLascaux_2$Mclust$Kappa, clustering_metricsLascaux_2$DPMVN$Kappa)
row_Lascaux_ARI <- data.frame(clustering_metricsLascaux_2$True$ARI, dp_Gmvlg_metricsLascaux_1$dp_Gmvlg$ARI, clustering_metricsLascaux_2$KMeans$ARI, clustering_metricsLascaux_2$CLARA$ARI, clustering_metricsLascaux_2$PAM$ARI, clustering_metricsLascaux_2$Hierarchical$ARI, clustering_metricsLascaux_2$Mclust$ARI, clustering_metricsLascaux_2$DPMVN$ARI)
row_Lascaux_NMI <- data.frame(clustering_metricsLascaux_2$True$NMI, dp_Gmvlg_metricsLascaux_1$dp_Gmvlg$NMI, clustering_metricsLascaux_2$KMeans$NMI, clustering_metricsLascaux_2$CLARA$NMI, clustering_metricsLascaux_2$PAM$NMI, clustering_metricsLascaux_2$Hierarchical$NMI, clustering_metricsLascaux_2$Mclust$NMI, clustering_metricsLascaux_2$DPMVN$NMI)
row_Lascaux_Cpu <- data.frame(dp_Gmvlg_metricsLascaux_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsLascaux_2$KMeans$CPU_RUN_TIME, clustering_metricsLascaux_2$CLARA$CPU_RUN_TIME, clustering_metricsLascaux_2$PAM$CPU_RUN_TIME, clustering_metricsLascaux_2$Hierarchical$CPU_RUN_TIME, clustering_metricsLascaux_2$Mclust$CPU_RUN_TIME, clustering_metricsLascaux_2$DPMVN$CPU_RUN_TIME)

tableLascaux_1 <- clustering_metricsLascaux_2$True$Table
tableLascaux_2 <- dp_Gmvlg_metricsLascaux_1$dp_Gmvlg$Table
tableLascaux_3 <- clustering_metricsLascaux_2$KMeans$Table
tableLascaux_4 <- clustering_metricsLascaux_2$CLARA$Table
tableLascaux_5 <- clustering_metricsLascaux_2$PAM$Table
tableLascaux_6 <- clustering_metricsLascaux_2$Hierarchical$Table
tableLascaux_7 <- clustering_metricsLascaux_2$Mclust$Table
tableLascaux_8 <- clustering_metricsLascaux_2$DPMVN$Table

colnames(row_Lascaux_K) <- NULL
colnames(row_Lascaux_ARI) <- NULL
colnames(row_Lascaux_NMI) <- NULL
colnames(row_Lascaux_Cpu) <- NULL

row_Lascaux_K <- as.matrix(row_Lascaux_K)
row_Lascaux_ARI <- as.matrix(row_Lascaux_ARI)
row_Lascaux_NMI <- as.matrix(row_Lascaux_NMI)
row_Lascaux_Cpu <- as.matrix(row_Lascaux_Cpu)

kappa_table_Lascaux <- rbind(row_Lascaux_K)
ARI_table_Lascaux <- rbind(row_Lascaux_ARI)
NMI_table_Lascaux <- rbind(row_Lascaux_NMI)
cpu_runtime_table_Lascaux <- rbind(row_Lascaux_Cpu)
clustering_table_Lascaux <- rbind(tableLascaux_1,tableLascaux_2,tableLascaux_3,tableLascaux_4,tableLascaux_5,tableLascaux_6,tableLascaux_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Lascaux)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Lascaux)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Lascaux)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Lascaux)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Lascaux)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/18_LASCAUX/Cluster_metricsLascaux.xlsx", overwrite = TRUE)

