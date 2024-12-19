#################################################################
#         AUSTRALIAN INSTITUTE OF SPORTS (AIS)  DATASET        ##
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

# Load the AIS dataset
data("ais")
help(ais)
Ais_data <- as.matrix(ais[,c(1:2,4:6,8:10)])
true_lab_Ais <- as.numeric(ais$sex)

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Ais_data)

# 2. Check for missing values
sum(is.na(Ais_data))

# 3. Distribution of each feature
par(mfrow=c(3, 3))  # Set up the plotting area for 11 histograms
for (i in 1:ncol(Ais_data)) {
  hist(Ais_data[, i], main=colnames(Ais_data)[i], xlab=colnames(Ais_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Ais_data, main="AIS Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Ais])

# 5. Boxplots for each feature by sex
# Add normalized data to the original data frame
ais_normalized <- as.data.frame(Ais_data)
ais_normalized$sex <- ais$sex

par(mfrow=c(3, 3))  # Reset the plotting area for 11 boxplots
for (i in 1:ncol(Ais_data)) {
  boxplot(Ais_data[, i] ~ ais$sex, main=colnames(Ais_data)[i], xlab="Sex", ylab=colnames(Ais_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Ais_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the AIS data
Ais_scaled_data <- scale(Ais_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Ais_scaled_data, 2, skewness)
kurtosis_values <- apply(Ais_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Ais_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Ais_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Ais_scaled_data[,i])$out
  print(paste("Feature:", colnames(Ais_scaled_data)[i]))
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
outliers_matrix <- apply(Ais_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Ais_scaled_data, cluster = true_lab_Ais), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Ais_scaled_data)
D <- ncol(Ais_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Ais_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Ais_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Ais_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(12,13)         
       lambda_g[j, g] ~ dgamma(1,3)    
    }
    alpha[g] ~ dbeta(2,6)
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
run_time_Ais<- toc()
codaSamples_Ais = as.mcmc.list( runJagsOut )
summaryChains_Ais <- summary(codaSamples_Ais)

diagMCMC( codaObject=codaSamples_Ais , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Ais , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Ais , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Ais , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Ais , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Ais , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Ais , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Ais , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Ais$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_Ais$statistics),1], P, G)
z_mode_Ais <- apply(matrix(summaryChains_Ais$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_Ais$statistics),1], P, G),1, which.max)
z_mode_Ais

par(mfrow=c(1, 1))
plot(Ais_scaled_data, col= z_mode_Ais, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ais), col=unique(z_mode_Ais), pch=16, title="Cluster")
table( true_lab_Ais , z_mode_Ais)

# To switch to the same labels in the true clustering
new <- 1:2
old <- c(2,1)
z_mode_Ais[z_mode_Ais %in% old] <- new[match(z_mode_Ais,old,nomatch = 0)]

plot(Ais_scaled_data, col= z_mode_Ais, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ais), col=unique(z_mode_Ais), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Ais , z_mode_Ais)
kappa2(data.frame(rater1 = true_lab_Ais, rater2 = z_mode_Ais))

calculate_dp_Gmvlg_clustering_metrics <- function(dataAis_1, true_clusters_Ais) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Ais <- table(true_clusters_Ais, z_mode_Ais)
  kappa_result_Ais <- kappa2(data.frame(rater1 = true_clusters_Ais, rater2 = z_mode_Ais))
  ari_result_Ais <- adjustedRandIndex(true_clusters_Ais, z_mode_Ais)
  nmi_result_Ais <- NMI(true_clusters_Ais, z_mode_Ais)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Ais,
    Kappa = kappa_result_Ais$value,
    ARI = ari_result_Ais,
    NMI = nmi_result_Ais,
    CPU_RUN_TIME = run_time_Ais$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsAis_1 <- calculate_dp_Gmvlg_clustering_metrics(dataAis_1 = Ais_scaled_data, true_clusters_Ais = true_lab_Ais)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsAis_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsAis_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsAis_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsAis_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsAis_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsAis_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Ais <- eclust(Ais_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Ais <- toc()
kmeans_clusters_Ais <- kmeans_result_Ais$cluster
table(true_lab_Ais, kmeans_clusters_Ais)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Ais <- eclust(Ais_scaled_data, "clara", G, graph = FALSE)
clara_time_Ais <- toc()
clara_clusters_Ais <- clara_result_Ais$cluster
table(true_lab_Ais, clara_clusters_Ais)

# PAM clustering
tic("PAM Runtime")
pam_result_Ais <- eclust(Ais_scaled_data, "pam", G, graph = FALSE)
pam_time_Ais <- toc()
pam_clusters_Ais <- pam_result_Ais$cluster
table(true_lab_Ais, pam_clusters_Ais)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Ais <- hclust(dist(Ais_scaled_data), method = "ward.D2")
hclust_time_Ais <- toc()
hclust_clusters_Ais <- cutree(hclust_result_Ais, k = G)
table(true_lab_Ais, hclust_clusters_Ais)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Ais <- Mclust(Ais_scaled_data, G = G)
mclust_time_Ais <- toc()
summary(mclust_result_Ais)
mclust_clusters_Ais <- mclust_result_Ais$classification
table(true_lab_Ais, mclust_clusters_Ais)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Ais_scaled_data, alphaPriors = c(1000,2))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Ais <- toc()
# Extract clusters 
dpMVN_clusters_Ais <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Ais)
table(true_lab_Ais, dpMVN_clusters_Ais)

calculate_clustering_metricsAis_2 <- function(dataAis_2, true_clusters_Ais, estimated_clusters_Ais_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Ais_list)) {
    clusters <- estimated_clusters_Ais_list[[method_name]]
    
    # Calculate metrics
    table_result_Ais <- table(true_clusters_Ais, clusters)
    kappa_result_Ais <- kappa2(data.frame(rater1 = true_clusters_Ais, rater2 = clusters))
    ari_result_Ais <- adjustedRandIndex(true_clusters_Ais, clusters)
    nmi_result_Ais <- NMI(true_clusters_Ais, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Ais,
      Kappa = kappa_result_Ais$value,
      ARI = ari_result_Ais,
      NMI = nmi_result_Ais,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Ais <- list(
  KMeans = kmeans_clusters_Ais,
  CLARA = clara_clusters_Ais,
  PAM = pam_clusters_Ais,
  Hierarchical = hclust_clusters_Ais,
  Mclust = mclust_clusters_Ais,
  DPMVN = dpMVN_clusters_Ais,
  True = true_lab_Ais
)

times_list_Ais <- list(
  KMeans = kmeans_time_Ais,
  CLARA = clara_time_Ais,
  PAM = pam_time_Ais,
  Hierarchical = hclust_time_Ais,
  Mclust = mclust_time_Ais,
  DPMVN = DPMVN_time_Ais
)

# Call the function
clustering_metricsAis_2 <- calculate_clustering_metricsAis_2(dataAis_2 = Ais_scaled_data, true_clusters_Ais = true_lab_Ais, estimated_clusters_Ais_list = cluster_result_Ais, times_list = times_list_Ais)

# Print the results for each method
for (method_name in names(clustering_metricsAis_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsAis_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsAis_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsAis_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsAis_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsAis_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Ais_K <- data.frame(clustering_metricsAis_2$True$Kappa, dp_Gmvlg_metricsAis_1$dp_Gmvlg$Kappa, clustering_metricsAis_2$KMeans$Kappa, clustering_metricsAis_2$CLARA$Kappa, clustering_metricsAis_2$PAM$Kappa, clustering_metricsAis_2$Hierarchical$Kappa, clustering_metricsAis_2$Mclust$Kappa, clustering_metricsAis_2$DPMVN$Kappa)
row_Ais_ARI <- data.frame(clustering_metricsAis_2$True$ARI, dp_Gmvlg_metricsAis_1$dp_Gmvlg$ARI, clustering_metricsAis_2$KMeans$ARI, clustering_metricsAis_2$CLARA$ARI, clustering_metricsAis_2$PAM$ARI, clustering_metricsAis_2$Hierarchical$ARI, clustering_metricsAis_2$Mclust$ARI, clustering_metricsAis_2$DPMVN$ARI)
row_Ais_NMI <- data.frame(clustering_metricsAis_2$True$NMI, dp_Gmvlg_metricsAis_1$dp_Gmvlg$NMI, clustering_metricsAis_2$KMeans$NMI, clustering_metricsAis_2$CLARA$NMI, clustering_metricsAis_2$PAM$NMI, clustering_metricsAis_2$Hierarchical$NMI, clustering_metricsAis_2$Mclust$NMI, clustering_metricsAis_2$DPMVN$NMI)
row_Ais_Cpu <- data.frame(dp_Gmvlg_metricsAis_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsAis_2$KMeans$CPU_RUN_TIME, clustering_metricsAis_2$CLARA$CPU_RUN_TIME, clustering_metricsAis_2$PAM$CPU_RUN_TIME, clustering_metricsAis_2$Hierarchical$CPU_RUN_TIME, clustering_metricsAis_2$Mclust$CPU_RUN_TIME, clustering_metricsAis_2$DPMVN$CPU_RUN_TIME)

tableAis_1 <- clustering_metricsAis_2$True$Table
tableAis_2 <- dp_Gmvlg_metricsAis_1$dp_Gmvlg$Table
tableAis_3 <- clustering_metricsAis_2$KMeans$Table
tableAis_4 <- clustering_metricsAis_2$CLARA$Table
tableAis_5 <- clustering_metricsAis_2$PAM$Table
tableAis_6 <- clustering_metricsAis_2$Hierarchical$Table
tableAis_7 <- clustering_metricsAis_2$Mclust$Table
tableAis_8 <- clustering_metricsAis_2$DPMVN$Table

colnames(row_Ais_K) <- NULL
colnames(row_Ais_ARI) <- NULL
colnames(row_Ais_NMI) <- NULL
colnames(row_Ais_Cpu) <- NULL

row_Ais_K <- as.matrix(row_Ais_K)
row_Ais_ARI <- as.matrix(row_Ais_ARI)
row_Ais_NMI <- as.matrix(row_Ais_NMI)
row_Ais_Cpu <- as.matrix(row_Ais_Cpu)

kappa_table_Ais <- rbind(row_Ais_K)
ARI_table_Ais <- rbind(row_Ais_ARI)
NMI_table_Ais <- rbind(row_Ais_NMI)
cpu_runtime_table_Ais <- rbind(row_Ais_Cpu)
clustering_table_Ais <- rbind(tableAis_1,tableAis_2,tableAis_3,tableAis_4,tableAis_5,tableAis_6,tableAis_7, tableAis_8)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Ais)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Ais)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Ais)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Ais)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Ais)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/17_AIS/Cluster_metricsAis.xlsx", overwrite = TRUE)
