#################################################################
#                       PAP DATASET                            ##
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

# Load the pap dataset
data("pap",package ="ade4")
pap1 <- data.frame(pap$tab,pap$taxo$superfamille)
names(pap1)[ncol(pap1)] <- "superfamille"
pap_data <- as.matrix(pap1[, -5])
true_lab_pap=as.numeric(pap1[,5])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(pap_data)

# 2. Check for missing values
sum(is.na(pap_data))

# 3. Distribution of each feature
par(mfrow=c(2, 2))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(pap_data)) {
  hist(pap_data[, i], main=colnames(pap_data)[i], xlab=colnames(pap_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(pap_data, main="PAP Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_pap])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(pap_data)
Cor_normalized$superfamille <- pap$taxo$superfamille

par(mfrow=c(2,2))  # Reset the plotting area for boxplots
for (i in 1:ncol(pap_data)) {
  boxplot(pap_data[, i] ~ Cor_normalized$superfamille, main=colnames(pap_data)[i], xlab="superfamille", ylab=colnames(pap_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(pap_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the pap dataset
pap_scaled_data <- scale(pap_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(pap_scaled_data, 2, skewness)
kurtosis_values <- apply(pap_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=1)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=1)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(pap_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(pap_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(pap_scaled_data[,i])$out
  print(paste("Feature:", colnames(pap_scaled_data)[i]))
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
outliers_matrix <- apply(pap_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = pap_scaled_data, cluster = true_lab_pap), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(pap_scaled_data)
D <- ncol(pap_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = pap_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(pap_scaled_data))^(1 / (D - 1)),
  G = G
)

deltaData <- det(cor(pap_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(15,18)          
       lambda_g[j, g] ~ dgamma(15,15) 
     }
   }
   alpha[1] ~ dgamma(10,20)
   alpha[2] ~ dgamma(10,15)
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
run_time_pap <- toc()
codaSamples_pap = as.mcmc.list( runJagsOut )
summaryChains_pap <- summary(codaSamples_pap)

diagMCMC( codaObject=codaSamples_pap , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_pap , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_pap , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_pap , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_pap , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_pap , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_pap , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_pap , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_pap$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_pap$statistics),1], P, G)
z_mode_pap <- apply(matrix(summaryChains_pap$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_pap$statistics),1], P, G),1, which.max)
z_mode_pap

plot(pap_scaled_data, col= z_mode_pap, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_pap), col=unique(z_mode_pap), pch=16, title="Cluster")
table( true_lab_pap , z_mode_pap)

# To switch to the same labels in the true clustering
#new <- 1:2
#old <- c(2,1)
#z_mode_pap[z_mode_pap %in% old] <- new[match(z_mode_pap,old,nomatch = 0)]

plot(pap_scaled_data, col= z_mode_pap, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_pap), col=unique(z_mode_pap), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_pap , z_mode_pap)
kappa2(data.frame(rater1 = true_lab_pap, rater2 = z_mode_pap))

calculate_dp_Gmvlg_clustering_metrics <- function(datapap_1, true_clusters_pap) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics

  # Calculate clustering metrics and the contingency table
  table_result_pap <- table(true_clusters_pap, z_mode_pap)
  kappa_result_pap <- kappa2(data.frame(rater1 = true_clusters_pap, rater2 = z_mode_pap))
  ari_result_pap <- adjustedRandIndex(true_clusters_pap, z_mode_pap)
  nmi_result_pap <- NMI(true_clusters_pap, z_mode_pap)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_pap,
    Kappa = kappa_result_pap$value,
    ARI = ari_result_pap,
    NMI = nmi_result_pap,
    CPU_RUN_TIME = run_time_pap$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricspap_1 <- calculate_dp_Gmvlg_clustering_metrics(datapap_1 = pap_scaled_data, true_clusters_pap = true_lab_pap)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricspap_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricspap_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricspap_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricspap_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricspap_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricspap_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_pap <- eclust(pap_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_pap <- toc()
kmeans_clusters_pap <- kmeans_result_pap$cluster
table(true_lab_pap, kmeans_clusters_pap)

# CLARA clustering
tic("CLARA Runtime")
clara_result_pap <- eclust(pap_scaled_data, "clara", G, graph = FALSE)
clara_time_pap <- toc()
clara_clusters_pap <- clara_result_pap$cluster
table(true_lab_pap, clara_clusters_pap)
new <- 1:2
old <- c(2,1)
clara_clusters_pap[clara_clusters_pap %in% old] <- new[match(clara_clusters_pap, old, nomatch = 0)]
table(true_lab_pap, clara_clusters_pap)

# PAM clustering
tic("PAM Runtime")
pam_result_pap <- eclust(pap_scaled_data, "pam", G, graph = FALSE)
pam_time_pap <- toc()
pam_clusters_pap <- pam_result_pap$cluster
table(true_lab_pap, pam_clusters_pap)
new <- 1:2
old <- c(2,1)
pam_clusters_pap[pam_clusters_pap %in% old] <- new[match(pam_clusters_pap, old, nomatch = 0)]
table(true_lab_pap, pam_clusters_pap)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_pap <- hclust(dist(pap_scaled_data), method = "ward.D2")
hclust_time_pap <- toc()
hclust_clusters_pap <- cutree(hclust_result_pap, k = G)
table(true_lab_pap, hclust_clusters_pap)
new <- 1:2
old <- c(2,1)
hclust_clusters_pap[hclust_clusters_pap %in% old] <- new[match(hclust_clusters_pap, old, nomatch = 0)]
table(true_lab_pap, hclust_clusters_pap)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_pap <- Mclust(pap_scaled_data, G = G)
mclust_time_pap <- toc()
summary(mclust_result_pap)
mclust_clusters_pap <- mclust_result_pap$classification
table(true_lab_pap, mclust_clusters_pap)
new <- 1:2
old <- c(2,1)
mclust_clusters_pap[mclust_clusters_pap %in% old] <- new[match(mclust_clusters_pap, old, nomatch = 0)]
table(true_lab_pap, mclust_clusters_pap)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(pap_scaled_data, alphaPriors = c(6,22))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_pap <- toc()
# Extract clusters 
dpMVN_clusters_pap <- as.numeric(dp$clusterLabels)
#new <- c(2,3)
#old <- c(3,2)
#dpMVN_clusters_pap[dpMVN_clusters_pap %in% old] <- new[match(dpMVN_clusters_pap, old, nomatch = 0)]
print(dpMVN_clusters_pap)
table(true_lab_pap, dpMVN_clusters_pap)

calculate_clustering_metricspap_2 <- function(datapap_2, true_clusters_pap, estimated_clusters_pap_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_pap_list)) {
    clusters <- estimated_clusters_pap_list[[method_name]]
    
    # Calculate metrics
    table_result_pap <- table(true_clusters_pap, clusters)
    kappa_result_pap <- kappa2(data.frame(rater1 = true_clusters_pap, rater2 = clusters))
    ari_result_pap <- adjustedRandIndex(true_clusters_pap, clusters)
    nmi_result_pap <- NMI(true_clusters_pap, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_pap,
      Kappa = kappa_result_pap$value,
      ARI = ari_result_pap,
      NMI = nmi_result_pap,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_pap <- list(
  KMeans = kmeans_clusters_pap,
  CLARA = clara_clusters_pap,
  PAM = pam_clusters_pap,
  Hierarchical = hclust_clusters_pap,
  Mclust = mclust_clusters_pap,
  DPMVN = dpMVN_clusters_pap,
  True = true_lab_pap
)

times_list_pap <- list(
  KMeans = kmeans_time_pap,
  CLARA = clara_time_pap,
  PAM = pam_time_pap,
  Hierarchical = hclust_time_pap,
  Mclust = mclust_time_pap,
  DPMVN = DPMVN_time_pap
)

# Call the function
clustering_metricspap_2 <- calculate_clustering_metricspap_2(datapap_2 = pap_scaled_data, true_clusters_pap = true_lab_pap, estimated_clusters_pap_list = cluster_result_pap, times_list = times_list_pap)

# Print the results for each method
for (method_name in names(clustering_metricspap_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricspap_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricspap_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricspap_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricspap_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricspap_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_pap_K <- data.frame(clustering_metricspap_2$True$Kappa, dp_Gmvlg_metricspap_1$dp_Gmvlg$Kappa, clustering_metricspap_2$KMeans$Kappa, clustering_metricspap_2$CLARA$Kappa, clustering_metricspap_2$PAM$Kappa, clustering_metricspap_2$Hierarchical$Kappa, clustering_metricspap_2$Mclust$Kappa, clustering_metricspap_2$DPMVN$Kappa)
row_pap_ARI <- data.frame(clustering_metricspap_2$True$ARI, dp_Gmvlg_metricspap_1$dp_Gmvlg$ARI, clustering_metricspap_2$KMeans$ARI, clustering_metricspap_2$CLARA$ARI, clustering_metricspap_2$PAM$ARI, clustering_metricspap_2$Hierarchical$ARI, clustering_metricspap_2$Mclust$ARI, clustering_metricspap_2$DPMVN$ARI)
row_pap_NMI <- data.frame(clustering_metricspap_2$True$NMI, dp_Gmvlg_metricspap_1$dp_Gmvlg$NMI, clustering_metricspap_2$KMeans$NMI, clustering_metricspap_2$CLARA$NMI, clustering_metricspap_2$PAM$NMI, clustering_metricspap_2$Hierarchical$NMI, clustering_metricspap_2$Mclust$NMI, clustering_metricspap_2$DPMVN$NMI)
row_pap_Cpu <- data.frame(dp_Gmvlg_metricspap_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricspap_2$KMeans$CPU_RUN_TIME, clustering_metricspap_2$CLARA$CPU_RUN_TIME, clustering_metricspap_2$PAM$CPU_RUN_TIME, clustering_metricspap_2$Hierarchical$CPU_RUN_TIME, clustering_metricspap_2$Mclust$CPU_RUN_TIME, clustering_metricspap_2$DPMVN$CPU_RUN_TIME)

tablepap_1 <- clustering_metricspap_2$True$Table
tablepap_2 <- dp_Gmvlg_metricspap_1$dp_Gmvlg$Table
tablepap_3 <- clustering_metricspap_2$KMeans$Table
tablepap_4 <- clustering_metricspap_2$CLARA$Table
tablepap_5 <- clustering_metricspap_2$PAM$Table
tablepap_6 <- clustering_metricspap_2$Hierarchical$Table
tablepap_7 <- clustering_metricspap_2$Mclust$Table
tablepap_8 <- clustering_metricspap_2$DPMVN$Table

colnames(row_pap_K) <- NULL
colnames(row_pap_ARI) <- NULL
colnames(row_pap_NMI) <- NULL
colnames(row_pap_Cpu) <- NULL

row_pap_K <- as.matrix(row_pap_K)
row_pap_ARI <- as.matrix(row_pap_ARI)
row_pap_NMI <- as.matrix(row_pap_NMI)
row_pap_Cpu <- as.matrix(row_pap_Cpu)

kappa_table_pap <- rbind(row_pap_K)
ARI_table_pap <- rbind(row_pap_ARI)
NMI_table_pap <- rbind(row_pap_NMI)
cpu_runtime_table_pap <- rbind(row_pap_Cpu)
clustering_table_pap <- rbind(tablepap_1,tablepap_2,tablepap_3,tablepap_4,tablepap_5,tablepap_6,tablepap_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_pap)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_pap)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_pap)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_pap)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_pap)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/08_PAP/Cluster_metricspap.xlsx", overwrite = TRUE)
