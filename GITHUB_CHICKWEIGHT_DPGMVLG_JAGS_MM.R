#################################################################
#                       Chick Weight DATASET                   ##
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
library(MASS)
library(aricode)
library(moments)

# Load the Chick Weight dataset
data("ChickWeight")
data1 <- sapply(ChickWeight, as.numeric) # Convert all columns to numeric
CW_data <- data1[,-c(3,4)]
true_lab_CW <- as.numeric(ChickWeight$Diet)    # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(CW_data)

# 2. Check for missing values
sum(is.na(CW_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(CW_data)) {
  hist(CW_data[, i], main=colnames(CW_data)[i], xlab=colnames(CW_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(CW_data, main="Chick Weight Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_CW])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
cw_normalized <- as.data.frame(CW_data)
cw_normalized$Diet <- ChickWeight$Diet

par(mfrow=c(1, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(CW_data)) {
  boxplot(CW_data[, i] ~ ChickWeight$Diet, main=colnames(CW_data)[i], xlab="Diet", ylab=colnames(CW_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(CW_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Chick Weight dataset
CW_scaled_data <- scale(CW_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(CW_scaled_data, 2, skewness)
kurtosis_values <- apply(CW_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(CW_scaled_data))
combined_data
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(CW_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(CW_scaled_data[,i])$out
  print(paste("Feature:", colnames(CW_scaled_data)[i]))
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
outliers_matrix <- apply(CW_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = CW_scaled_data, cluster = true_lab_CW), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(CW_scaled_data)
D <- ncol(CW_scaled_data)

#Try with different number of clusters
G <- 4  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = CW_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(CW_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(CW_scaled_data))^(1 / (D - 1))

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
      logLik[i, g] <- z[i, g] * (log(pi_g[g])+ v_g * log(delta_g) + sum(log(mu_g[ ,g])) - v_g * sum(log(lambda_g[, g])) - loggam(v_g) ) + # Unnecessary: - loggam(v_g)
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
       lambda_g[j, g] ~ dgamma(1,5)
    }
     alpha[g] ~ dgamma(1,5)
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
run_time_CW <- toc()
codaSamples_CW = as.mcmc.list( runJagsOut )
summaryChains_CW <- summary(codaSamples_CW)

diagMCMC( codaObject=codaSamples_CW , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_CW , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_CW , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_CW , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_CW , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_CW , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_CW , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_CW , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_CW$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_CW$statistics),1], P, G)
z_mode_CW <- apply(matrix(summaryChains_CW$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_CW$statistics),1], P, G),1, which.max)
z_mode_CW

plot(CW_scaled_data, col= z_mode_CW, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_CW), col=unique(z_mode_CW), pch=16, title="Cluster")
table( true_lab_CW , z_mode_CW)

# To switch to the same labels in the true clustering
new <- c(1,4)
old <- c(4,1)
z_mode_CW[z_mode_CW %in% old] <- new[match(z_mode_CW,old,nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
z_mode_CW[z_mode_CW %in% old] <- new[match(z_mode_CW,old,nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
z_mode_CW[z_mode_CW %in% old] <- new[match(z_mode_CW,old,nomatch = 0)]

plot(CW_scaled_data, col= z_mode_CW, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_CW), col=unique(z_mode_CW), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_CW , z_mode_CW)
kappa2(data.frame(rater1 = true_lab_CW, rater2 = z_mode_CW))

calculate_dp_Gmvlg_clustering_metrics <- function(dataCW_1, true_clusters_CW) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_CW <- table(true_clusters_CW, z_mode_CW)
  kappa_result_CW <- kappa2(data.frame(rater1 = true_clusters_CW, rater2 = z_mode_CW))
  ari_result_CW <- adjustedRandIndex(true_clusters_CW, z_mode_CW)
  nmi_result_CW <- NMI(true_clusters_CW, z_mode_CW)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_CW,
    Kappa = kappa_result_CW$value,
    ARI = ari_result_CW,
    NMI = nmi_result_CW,
    CPU_RUN_TIME = run_time_CW$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsCW_1 <- calculate_dp_Gmvlg_clustering_metrics(dataCW_1 = CW_scaled_data, true_clusters_CW = true_lab_CW)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsCW_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsCW_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsCW_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsCW_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsCW_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsCW_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_CW <- eclust(CW_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_CW <- toc()
kmeans_clusters_CW <- kmeans_result_CW$cluster
table(true_lab_CW, kmeans_clusters_CW)
new <- c(1,2)
old <- c(2,1)
kmeans_clusters_CW[kmeans_clusters_CW %in% old] <- new[match(kmeans_clusters_CW, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_CW[kmeans_clusters_CW %in% old] <- new[match(kmeans_clusters_CW, old, nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
kmeans_clusters_CW[kmeans_clusters_CW %in% old] <- new[match(kmeans_clusters_CW, old, nomatch = 0)]
table(true_lab_CW, kmeans_clusters_CW)

# CLARA clustering
tic("CLARA Runtime")
clara_result_CW <- eclust(CW_scaled_data, "clara", G, graph = FALSE)
clara_time_CW <- toc()
clara_clusters_CW <- clara_result_CW$cluster
table(true_lab_CW, clara_clusters_CW)
new <- c(3,4)
old <- c(4,3)
clara_clusters_CW[clara_clusters_CW %in% old] <- new[match(clara_clusters_CW, old, nomatch = 0)]
table(true_lab_CW, clara_clusters_CW)

# PAM clustering
tic("PAM Runtime")
pam_result_CW <- eclust(CW_scaled_data, "pam", G, graph = FALSE)
pam_time_CW <- toc()
pam_clusters_CW <- pam_result_CW$cluster
new <- c(1,2)
old <- c(2,1)
pam_clusters_CW[pam_clusters_CW %in% old] <- new[match(pam_clusters_CW, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
pam_clusters_CW[pam_clusters_CW %in% old] <- new[match(pam_clusters_CW, old, nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
pam_clusters_CW[pam_clusters_CW %in% old] <- new[match(pam_clusters_CW, old, nomatch = 0)]
table(true_lab_CW, pam_clusters_CW)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_CW <- hclust(dist(CW_scaled_data), method = "ward.D2")
hclust_time_CW <- toc()
hclust_clusters_CW <- cutree(hclust_result_CW, k = G)
new <- c(3,4)
old <- c(4,3)
hclust_clusters_CW[hclust_clusters_CW %in% old] <- new[match(hclust_clusters_CW, old, nomatch = 0)]
table(true_lab_CW, hclust_clusters_CW)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_CW <- Mclust(CW_scaled_data, G = G)
mclust_time_CW <- toc()
summary(mclust_result_CW)
mclust_clusters_CW <- mclust_result_CW$classification
table(true_lab_CW, mclust_clusters_CW)
new <- c(1,3)
old <- c(3,1)
mclust_clusters_CW[mclust_clusters_CW %in% old] <- new[match(mclust_clusters_CW, old, nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
mclust_clusters_CW[mclust_clusters_CW %in% old] <- new[match(mclust_clusters_CW, old, nomatch = 0)]
table(true_lab_CW, mclust_clusters_CW)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(CW_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_CW <- toc()
# Extract clusters 
dpMVN_clusters_CW <- as.numeric(dp$clusterLabels)
new <- c(1,2)
old <- c(2,1)
dpMVN_clusters_CW[dpMVN_clusters_CW %in% old] <- new[match(dpMVN_clusters_CW, old, nomatch = 0)]
#new <- c(3,4)
#old <- c(4,3)
#dpMVN_clusters_CW[dpMVN_clusters_CW %in% old] <- new[match(dpMVN_clusters_CW, old, nomatch = 0)]
print(dpMVN_clusters_CW)
table(true_lab_CW, dpMVN_clusters_CW)
kappa2(data.frame(rater1 = true_lab_CW, rater2 = dpMVN_clusters_CW))

calculate_clustering_metricsCW_2 <- function(dataCW_2, true_clusters_CW, estimated_clusters_CW_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_CW_list)) {
    clusters <- estimated_clusters_CW_list[[method_name]]
    
    # Calculate metrics
    table_result_CW <- table(true_clusters_CW, clusters)
    kappa_result_CW <- kappa2(data.frame(rater1 = true_clusters_CW, rater2 = clusters))
    ari_result_CW <- adjustedRandIndex(true_clusters_CW, clusters)
    nmi_result_CW <- NMI(true_clusters_CW, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_CW,
      Kappa = kappa_result_CW$value,
      ARI = ari_result_CW,
      NMI = nmi_result_CW,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_CW <- list(
  KMeans = kmeans_clusters_CW,
  CLARA = clara_clusters_CW,
  PAM = pam_clusters_CW,
  Hierarchical = hclust_clusters_CW,
  Mclust = mclust_clusters_CW,
  DPMVN = dpMVN_clusters_CW,
  True = true_lab_CW
)

times_list_CW <- list(
  KMeans = kmeans_time_CW,
  CLARA = clara_time_CW,
  PAM = pam_time_CW,
  Hierarchical = hclust_time_CW,
  Mclust = mclust_time_CW,
  DPMVN = DPMVN_time_CW
)

# Call the function
clustering_metricsCW_2 <- calculate_clustering_metricsCW_2(dataCW_2 = CW_scaled_data, true_clusters_CW = true_lab_CW, estimated_clusters_CW_list = cluster_result_CW, times_list = times_list_CW)

# Print the results for each method
for (method_name in names(clustering_metricsCW_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsCW_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsCW_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsCW_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsCW_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsCW_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_CW_K <- data.frame(clustering_metricsCW_2$True$Kappa, dp_Gmvlg_metricsCW_1$dp_Gmvlg$Kappa, clustering_metricsCW_2$KMeans$Kappa, clustering_metricsCW_2$CLARA$Kappa, clustering_metricsCW_2$PAM$Kappa, clustering_metricsCW_2$Hierarchical$Kappa, clustering_metricsCW_2$Mclust$Kappa, clustering_metricsCW_2$DPMVN$Kappa)
row_CW_ARI <- data.frame(clustering_metricsCW_2$True$ARI, dp_Gmvlg_metricsCW_1$dp_Gmvlg$ARI, clustering_metricsCW_2$KMeans$ARI, clustering_metricsCW_2$CLARA$ARI, clustering_metricsCW_2$PAM$ARI, clustering_metricsCW_2$Hierarchical$ARI, clustering_metricsCW_2$Mclust$ARI, clustering_metricsCW_2$DPMVN$ARI)
row_CW_NMI <- data.frame(clustering_metricsCW_2$True$NMI, dp_Gmvlg_metricsCW_1$dp_Gmvlg$NMI, clustering_metricsCW_2$KMeans$NMI, clustering_metricsCW_2$CLARA$NMI, clustering_metricsCW_2$PAM$NMI, clustering_metricsCW_2$Hierarchical$NMI, clustering_metricsCW_2$Mclust$NMI, clustering_metricsCW_2$DPMVN$NMI)
row_CW_CPU <- data.frame(dp_Gmvlg_metricsCW_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsCW_2$KMeans$CPU_RUN_TIME, clustering_metricsCW_2$CLARA$CPU_RUN_TIME, clustering_metricsCW_2$PAM$CPU_RUN_TIME, clustering_metricsCW_2$Hierarchical$CPU_RUN_TIME, clustering_metricsCW_2$Mclust$CPU_RUN_TIME, clustering_metricsCW_2$DPMVN$CPU_RUN_TIME)

tableCW_1 <- clustering_metricsCW_2$True$Table
tableCW_2 <- dp_Gmvlg_metricsCW_1$dp_Gmvlg$Table
tableCW_3 <- clustering_metricsCW_2$KMeans$Table
tableCW_4 <- clustering_metricsCW_2$CLARA$Table
tableCW_5 <- clustering_metricsCW_2$PAM$Table
tableCW_6 <- clustering_metricsCW_2$Hierarchical$Table
tableCW_7 <- clustering_metricsCW_2$Mclust$Table
tableCW_8 <- clustering_metricsCW_2$DPMVN$Table

colnames(row_CW_K) <- NULL
colnames(row_CW_ARI) <- NULL
colnames(row_CW_NMI) <- NULL
colnames(row_CW_CPU) <- NULL

row_CW_K <- as.matrix(row_CW_K)
row_CW_ARI <- as.matrix(row_CW_ARI)
row_CW_NMI <- as.matrix(row_CW_NMI)
row_CW_CPU <- as.matrix(row_CW_CPU)

kappa_table_CW <- rbind(row_CW_K)
ARI_table_CW <- rbind(row_CW_ARI)
NMI_table_CW <- rbind(row_CW_NMI)
CPU_runtime_table_CW <- rbind(row_CW_CPU)
clustering_table_CW <- rbind(tableCW_1,tableCW_2,tableCW_3,tableCW_4,tableCW_5,tableCW_6,tableCW_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_CW)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_CW)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_CW)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", CPU_runtime_table_CW)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_CW)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/20_CHICK WEIGHT/Cluster_metricsCW.xlsx", overwrite = TRUE)
