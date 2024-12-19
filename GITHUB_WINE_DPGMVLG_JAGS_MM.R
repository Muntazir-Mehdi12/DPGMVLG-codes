#################################################################
#                       Wine DATASET                           ##
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
library(rattle)
library(pgmm)
library(moments)
library(aricode)

# Load the Wine dataset
data("wine",package = "rattle")
help("wine",package = "rattle")
wine_data <- as.matrix(wine[,-1]) # Convert all columns to numeric
true_lab_wine <- as.numeric(wine$Type)    # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(wine_data)

# 2. Check for missing values
sum(is.na(wine_data))

# 3. Distribution of each feature
par(mfrow=c(4, 4))  # Adjusted for the number of features
for (i in 1:ncol(wine_data)) {
  hist(wine_data[, i], main=colnames(wine_data)[i], xlab=colnames(wine_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(wine_data, main="Wine Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[true_lab_wine])

# 5. Boxplots for each feature by 'Type'
# Add normalized data to the original data frame
wine_normalized <- as.data.frame(wine_data)
wine_normalized$Type <- as.factor(wine$Type)

par(mfrow=c(4, 4))  # Reset the plotting area for boxplots
for (i in 1:ncol(wine_data)) {
  boxplot(wine_data[, i] ~ wine_normalized$Type, main=colnames(wine_data)[i], xlab="Type", ylab=colnames(wine_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(wine_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Wine dataset
wine_scaled_data <- scale(wine_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(wine_scaled_data, 2, skewness)
kurtosis_values <- apply(wine_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=2)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=2)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(wine_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(wine_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(wine_scaled_data[,i])$out
  print(paste("Feature:", colnames(wine_scaled_data)[i]))
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
outliers_matrix <- apply(wine_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = wine_scaled_data, cluster = true_lab_wine), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(wine_scaled_data)
D <- ncol(wine_scaled_data)

#Try with different number of clusters
G <- 3  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = wine_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(wine_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(wine_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(1,2)        
       lambda_g[j, g] ~ dgamma(2,2)
    }
    alpha[g] ~ dbeta(10,5)
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
run_time_wine <- toc()
codaSamples_wine = as.mcmc.list( runJagsOut )
summaryChains_wine <- summary(codaSamples_wine)

diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_wine$statistics[(1+13*G+13*G+P+1):nrow(summaryChains_wine$statistics),1], P, G)
z_mode_wine <- apply(matrix(summaryChains_wine$statistics[(1+13*G+13*G+P+1):nrow(summaryChains_wine$statistics),1], P, G),1, which.max)
z_mode_wine

plot(wine_scaled_data, col= z_mode_wine, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_wine), col=unique(z_mode_wine), pch=16, title="Cluster")
table( true_lab_wine , z_mode_wine)

# To switch to the same labels in the true clustering
new <- c(1,3)
old <- c(3,1)
z_mode_wine[z_mode_wine %in% old] <- new[match(z_mode_wine,old,nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
z_mode_wine[z_mode_wine %in% old] <- new[match(z_mode_wine,old,nomatch = 0)]

plot(wine_scaled_data, col= z_mode_wine, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_wine), col=unique(z_mode_wine), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_wine , z_mode_wine)
kappa2(data.frame(rater1 = true_lab_wine, rater2 = z_mode_wine))

calculate_dp_Gmvlg_clustering_metrics <- function(datawine_1, true_clusters_wine) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_wine <- table(true_clusters_wine, z_mode_wine)
  kappa_result_wine <- kappa2(data.frame(rater1 = true_clusters_wine, rater2 = z_mode_wine))
  ari_result_wine <- adjustedRandIndex(true_clusters_wine, z_mode_wine)
  nmi_result_wine <- NMI(true_clusters_wine, z_mode_wine)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_wine,
    Kappa = kappa_result_wine$value,
    ARI = ari_result_wine,
    NMI = nmi_result_wine,
    CPU_RUN_TIME = run_time_wine$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricswine_1 <- calculate_dp_Gmvlg_clustering_metrics(datawine_1 = wine_scaled_data, true_clusters_wine = true_lab_wine)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricswine_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_wine <- eclust(wine_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_wine <- toc()
kmeans_clusters_wine <- kmeans_result_wine$cluster
table(true_lab_wine, kmeans_clusters_wine)
new <- 1:2
old <- c(2,1)
kmeans_clusters_wine[kmeans_clusters_wine %in% old] <- new[match(kmeans_clusters_wine, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_wine[kmeans_clusters_wine %in% old] <- new[match(kmeans_clusters_wine, old, nomatch = 0)]
table(true_lab_wine, kmeans_clusters_wine)

# CLARA clustering
tic("CLARA Runtime")
clara_result_wine <- eclust(wine_scaled_data, "clara", G, graph = FALSE)
clara_time_wine <- toc()
clara_clusters_wine <- clara_result_wine$cluster
table(true_lab_wine, clara_clusters_wine)

# PAM clustering
tic("PAM Runtime")
pam_result_wine <- eclust(wine_scaled_data, "pam", G, graph = FALSE)
pam_time_wine <- toc()
pam_clusters_wine <- pam_result_wine$cluster
table(true_lab_wine, pam_clusters_wine)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_wine <- hclust(dist(wine_scaled_data), method = "ward.D2")
hclust_time_wine <- toc()
hclust_clusters_wine <- cutree(hclust_result_wine, k = G)
table(true_lab_wine, hclust_clusters_wine)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_wine <- Mclust(wine_scaled_data, G = G)
mclust_time_wine <- toc()
summary(mclust_result_wine)
mclust_clusters_wine <- mclust_result_wine$classification
table(true_lab_wine, mclust_clusters_wine)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(wine_scaled_data, alphaPriors = c(50,30))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_wine <- toc()
# Extract clusters 
dpMVN_clusters_wine <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_wine)
table(true_lab_wine, dpMVN_clusters_wine)

calculate_clustering_metricswine_2 <- function(datawine_2, true_clusters_wine, estimated_clusters_wine_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_wine_list)) {
    clusters <- estimated_clusters_wine_list[[method_name]]
    
    # Calculate metrics
    table_result_wine <- table(true_clusters_wine, clusters)
    kappa_result_wine <- kappa2(data.frame(rater1 = true_clusters_wine, rater2 = clusters))
    ari_result_wine <- adjustedRandIndex(true_clusters_wine, clusters)
    nmi_result_wine <- NMI(true_clusters_wine, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_wine,
      Kappa = kappa_result_wine$value,
      ARI = ari_result_wine,
      NMI = nmi_result_wine,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_wine <- list(
  KMeans = kmeans_clusters_wine,
  CLARA = clara_clusters_wine,
  PAM = pam_clusters_wine,
  Hierarchical = hclust_clusters_wine,
  Mclust = mclust_clusters_wine,
  DPMVN = dpMVN_clusters_wine,
  True = true_lab_wine
)

times_list_wine <- list(
  KMeans = kmeans_time_wine,
  CLARA = clara_time_wine,
  PAM = pam_time_wine,
  Hierarchical = hclust_time_wine,
  Mclust = mclust_time_wine,
  DPMVN = DPMVN_time_wine
)

# Call the function
clustering_metricswine_2 <- calculate_clustering_metricswine_2(datawine_2 = wine_scaled_data, true_clusters_wine = true_lab_wine, estimated_clusters_wine_list = cluster_result_wine, times_list = times_list_wine)

# Print the results for each method
for (method_name in names(clustering_metricswine_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricswine_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricswine_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricswine_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricswine_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricswine_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_wine_K <- data.frame(clustering_metricswine_2$True$Kappa, dp_Gmvlg_metricswine_1$dp_Gmvlg$Kappa, clustering_metricswine_2$KMeans$Kappa, clustering_metricswine_2$CLARA$Kappa, clustering_metricswine_2$PAM$Kappa, clustering_metricswine_2$Hierarchical$Kappa, clustering_metricswine_2$Mclust$Kappa, clustering_metricswine_2$DPMVN$Kappa)
row_wine_ARI <- data.frame(clustering_metricswine_2$True$ARI, dp_Gmvlg_metricswine_1$dp_Gmvlg$ARI, clustering_metricswine_2$KMeans$ARI, clustering_metricswine_2$CLARA$ARI, clustering_metricswine_2$PAM$ARI, clustering_metricswine_2$Hierarchical$ARI, clustering_metricswine_2$Mclust$ARI, clustering_metricswine_2$DPMVN$ARI)
row_wine_NMI <- data.frame(clustering_metricswine_2$True$NMI, dp_Gmvlg_metricswine_1$dp_Gmvlg$NMI, clustering_metricswine_2$KMeans$NMI, clustering_metricswine_2$CLARA$NMI, clustering_metricswine_2$PAM$NMI, clustering_metricswine_2$Hierarchical$NMI, clustering_metricswine_2$Mclust$NMI, clustering_metricswine_2$DPMVN$NMI)
row_wine_Cpu <- data.frame(dp_Gmvlg_metricswine_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricswine_2$KMeans$CPU_RUN_TIME, clustering_metricswine_2$CLARA$CPU_RUN_TIME, clustering_metricswine_2$PAM$CPU_RUN_TIME, clustering_metricswine_2$Hierarchical$CPU_RUN_TIME, clustering_metricswine_2$Mclust$CPU_RUN_TIME, clustering_metricswine_2$DPMVN$CPU_RUN_TIME)

tablewine_1 <- clustering_metricswine_2$True$Table
tablewine_2 <- dp_Gmvlg_metricswine_1$dp_Gmvlg$Table
tablewine_3 <- clustering_metricswine_2$KMeans$Table
tablewine_4 <- clustering_metricswine_2$CLARA$Table
tablewine_5 <- clustering_metricswine_2$PAM$Table
tablewine_6 <- clustering_metricswine_2$Hierarchical$Table
tablewine_7 <- clustering_metricswine_2$Mclust$Table
tablewine_8 <- clustering_metricswine_2$DPMVN$Table

colnames(row_wine_K) <- NULL
colnames(row_wine_ARI) <- NULL
colnames(row_wine_NMI) <- NULL
colnames(row_wine_Cpu) <- NULL

row_wine_K <- as.matrix(row_wine_K)
row_wine_ARI <- as.matrix(row_wine_ARI)
row_wine_NMI <- as.matrix(row_wine_NMI)
row_wine_Cpu <- as.matrix(row_wine_Cpu)

kappa_table_wine <- rbind(row_wine_K)
ARI_table_wine <- rbind(row_wine_ARI)
NMI_table_wine <- rbind(row_wine_NMI)
cpu_runtime_table_wine <- rbind(row_wine_Cpu)
clustering_table_wine <- rbind(tablewine_1,tablewine_2,tablewine_3,tablewine_4,tablewine_5,tablewine_6,tablewine_7, tablewine_8)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_wine)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_wine)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_wine)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_wine)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_wine)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/16_WINE/Cluster_metricswine.xlsx", overwrite = TRUE)
