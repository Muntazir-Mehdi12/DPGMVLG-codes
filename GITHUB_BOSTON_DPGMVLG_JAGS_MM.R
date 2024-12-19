#################################################################
#                        Boston Housing DATASET                ##
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
library(factoextra)
library(aricode)
library(moments)

# Load the Boston dataset
data("Boston")
Boston$chas[Boston$chas == 0] = "non-tract B.R"
Boston$chas[Boston$chas == 1] = "tract B.R"
Bos_data <- as.matrix(Boston[,-4])  # Exclude the 'class' column for clustering
true_lab_Bos <- as.numeric(factor(Boston[,4], levels = c("non-tract B.R", "tract B.R")))     # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Bos_data)

# 2. Check for missing values
sum(is.na(Bos_data))

# 3. Distribution of each feature
par(mfrow=c(4, 4))  # Set up the plotting area for 13 histograms
for (i in 1:ncol(Bos_data)) {
  hist(Bos_data[, i], main=colnames(Bos_data)[i], xlab=colnames(Bos_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Bos_data, main="Boston Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[true_lab_Bos %% 3 + 1])

# 5. Boxplots for each feature by 'rad' category
# Add normalized data to the original data frame
bos_normalized <- as.data.frame(Bos_data)
bos_normalized$chas <- Boston$chas


par(mfrow=c(4, 4))  # Reset the plotting area for boxplots
for (i in 1:ncol(Bos_data)) {
  boxplot(Bos_data[, i] ~ Boston$chas, main=colnames(Bos_data)[i], xlab="Radial Highway", ylab=colnames(Bos_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Bos_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Boston data
Bos_scaled_data <- scale(Bos_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Bos_scaled_data, 2, skewness)
kurtosis_values <- apply(Bos_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

par(mfrow=c(1, 2))
# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=2)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=2)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(Bos_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Bos_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Bos_scaled_data[,i])$out
  print(paste("Feature:", colnames(Bos_scaled_data)[i]))
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
outliers_matrix <- apply(Bos_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Bos_scaled_data, cluster = true_lab_Bos), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Bos_scaled_data)
D <- ncol(Bos_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Bos_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Bos_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Bos_scaled_data))^(1 / (D - 1))

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
     alpha[g] ~ dgamma(2,6)
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
run_time_Bos <- toc()
codaSamples_Bos = as.mcmc.list( runJagsOut )
summaryChains_Bos <- summary(codaSamples_Bos)

diagMCMC( codaObject=codaSamples_Bos , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Bos , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Bos , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Bos , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Bos , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Bos , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Bos , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Bos , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Bos$statistics[(1+13*G+13*G+P+1):nrow(summaryChains_Bos$statistics),1], P, G)
z_mode_Bos <- apply(matrix(summaryChains_Bos$statistics[(1+13*G+13*G+P+1):nrow(summaryChains_Bos$statistics),1], P, G),1, which.max)
z_mode_Bos

plot(Bos_scaled_data, col= z_mode_Bos, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Bos), col=unique(z_mode_Bos), pch=16, title="Cluster")
table( true_lab_Bos , z_mode_Bos)

# To switch to the same labels in the true clustering
new <- 1:2
old <- c(2,1)
z_mode_Bos[z_mode_Bos %in% old] <- new[match(z_mode_Bos,old,nomatch = 0)]

plot(Bos_scaled_data, col= z_mode_Bos, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Bos), col=unique(z_mode_Bos), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Bos , z_mode_Bos)
kappa2(data.frame(rater1 = true_lab_Bos, rater2 = z_mode_Bos))

calculate_dp_Gmvlg_clustering_metrics <- function(dataBos_1, true_clusters_Bos) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Bos <- table(true_clusters_Bos, z_mode_Bos)
  kappa_result_Bos <- kappa2(data.frame(rater1 = true_clusters_Bos, rater2 = z_mode_Bos))
  ari_result_Bos <- adjustedRandIndex(true_clusters_Bos, z_mode_Bos)
  nmi_result_Bos <- NMI(true_clusters_Bos, z_mode_Bos)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Bos,
    Kappa = kappa_result_Bos$value,
    ARI = ari_result_Bos,
    NMI = nmi_result_Bos,
    CPU_RUN_TIME = run_time_Bos$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsBos_1 <- calculate_dp_Gmvlg_clustering_metrics(dataBos_1 = Bos_scaled_data, true_clusters_Bos = true_lab_Bos)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsBos_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsBos_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsBos_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsBos_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsBos_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsBos_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Bos <- eclust(Bos_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Bos <- toc()
kmeans_clusters_Bos <- kmeans_result_Bos$cluster
table(true_lab_Bos, kmeans_clusters_Bos)
#new <- 1:2
#old <- c(2,1)
#kmeans_clusters_Bos[kmeans_clusters_Bos %in% old] <- new[match(kmeans_clusters_Bos, old, nomatch = 0)]
#table(true_lab_Bos, kmeans_clusters_Bos)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Bos <- eclust(Bos_scaled_data, "clara", G, graph = FALSE)
clara_time_Bos <- toc()
clara_clusters_Bos <- clara_result_Bos$cluster
table(true_lab_Bos, clara_clusters_Bos)
#new <- 1:2
#old <- c(2,1)
#clara_clusters_Bos[clara_clusters_Bos %in% old] <- new[match(clara_clusters_Bos, old, nomatch = 0)]
#table(true_lab_Bos, clara_clusters_Bos)

# PAM clustering
tic("PAM Runtime")
pam_result_Bos <- eclust(Bos_scaled_data, "pam", G, graph = FALSE)
pam_time_Bos <- toc()
pam_clusters_Bos <- pam_result_Bos$cluster
table(true_lab_Bos, pam_clusters_Bos)
#new <- 1:2
#old <- c(2,1)
#pam_clusters_Bos[pam_clusters_Bos %in% old] <- new[match(pam_clusters_Bos, old, nomatch = 0)]
#table(true_lab_Bos, pam_clusters_Bos)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Bos <- hclust(dist(Bos_scaled_data), method = "ward.D2")
hclust_time_Bos <- toc()
hclust_clusters_Bos <- cutree(hclust_result_Bos, k = G)
table(true_lab_Bos, hclust_clusters_Bos)
#new <- 1:2
#old <- c(2,1)
#hclust_clusters_Bos[hclust_clusters_Bos %in% old] <- new[match(hclust_clusters_Bos, old, nomatch = 0)]
#table(true_lab_Bos, hclust_clusters_Bos)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Bos <- Mclust(Bos_scaled_data, G = G)
mclust_time_Bos <- toc()
summary(mclust_result_Bos)
mclust_clusters_Bos <- mclust_result_Bos$classification
table(true_lab_Bos, mclust_clusters_Bos)
new <- 1:2
old <- c(2,1)
mclust_clusters_Bos[mclust_clusters_Bos %in% old] <- new[match(mclust_clusters_Bos, old, nomatch = 0)]
table(true_lab_Bos, mclust_clusters_Bos)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Bos_scaled_data, alphaPriors =  c(1,500))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 100)
DPMVN_time_Bos <- toc()
# Extract clusters 
dpMVN_clusters_Bos <- as.numeric(dp$clusterLabels)
new <- 1:2
old <- c(2,1)
dpMVN_clusters_Bos[dpMVN_clusters_Bos %in% old] <- new[match(dpMVN_clusters_Bos, old, nomatch = 0)]
print(dpMVN_clusters_Bos)
table(true_lab_Bos, dpMVN_clusters_Bos)

calculate_clustering_metricsBos_2 <- function(dataBos_2, true_clusters_Bos, estimated_clusters_Bos_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Bos_list)) {
    clusters <- estimated_clusters_Bos_list[[method_name]]
    
    # Calculate metrics
    table_result_Bos <- table(true_clusters_Bos, clusters)
    kappa_result_Bos <- kappa2(data.frame(rater1 = true_clusters_Bos, rater2 = clusters))
    ari_result_Bos <- adjustedRandIndex(true_clusters_Bos, clusters)
    nmi_result_Bos <- NMI(true_clusters_Bos, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Bos,
      Kappa = kappa_result_Bos$value,
      ARI = ari_result_Bos,
      NMI = nmi_result_Bos,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Bos <- list(
  KMeans = kmeans_clusters_Bos,
  CLARA = clara_clusters_Bos,
  PAM = pam_clusters_Bos,
  Hierarchical = hclust_clusters_Bos,
  Mclust = mclust_clusters_Bos,
  DPMVN = dpMVN_clusters_Bos,
  True = true_lab_Bos
)

times_list_Bos <- list(
  KMeans = kmeans_time_Bos,
  CLARA = clara_time_Bos,
  PAM = pam_time_Bos,
  Hierarchical = hclust_time_Bos,
  Mclust = mclust_time_Bos,
  DPMVN = DPMVN_time_Bos
)

# Call the function
clustering_metricsBos_2 <- calculate_clustering_metricsBos_2(dataBos_2 = Bos_scaled_data, true_clusters_Bos = true_lab_Bos, estimated_clusters_Bos_list = cluster_result_Bos, times_list = times_list_Bos)

# Print the results for each method
for (method_name in names(clustering_metricsBos_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsBos_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsBos_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsBos_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsBos_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsBos_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Bos_K <- data.frame(clustering_metricsBos_2$True$Kappa, dp_Gmvlg_metricsBos_1$dp_Gmvlg$Kappa, clustering_metricsBos_2$KMeans$Kappa, clustering_metricsBos_2$CLARA$Kappa, clustering_metricsBos_2$PAM$Kappa, clustering_metricsBos_2$Hierarchical$Kappa, clustering_metricsBos_2$Mclust$Kappa, clustering_metricsBos_2$DPMVN$Kappa)
row_Bos_ARI <- data.frame(clustering_metricsBos_2$True$ARI, dp_Gmvlg_metricsBos_1$dp_Gmvlg$ARI, clustering_metricsBos_2$KMeans$ARI, clustering_metricsBos_2$CLARA$ARI, clustering_metricsBos_2$PAM$ARI, clustering_metricsBos_2$Hierarchical$ARI, clustering_metricsBos_2$Mclust$ARI, clustering_metricsBos_2$DPMVN$ARI)
row_Bos_NMI <- data.frame(clustering_metricsBos_2$True$NMI, dp_Gmvlg_metricsBos_1$dp_Gmvlg$NMI, clustering_metricsBos_2$KMeans$NMI, clustering_metricsBos_2$CLARA$NMI, clustering_metricsBos_2$PAM$NMI, clustering_metricsBos_2$Hierarchical$NMI, clustering_metricsBos_2$Mclust$NMI, clustering_metricsBos_2$DPMVN$NMI)
row_Bos_CPU <- data.frame(dp_Gmvlg_metricsBos_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsBos_2$KMeans$CPU_RUN_TIME, clustering_metricsBos_2$CLARA$CPU_RUN_TIME, clustering_metricsBos_2$PAM$CPU_RUN_TIME, clustering_metricsBos_2$Hierarchical$CPU_RUN_TIME, clustering_metricsBos_2$Mclust$CPU_RUN_TIME, clustering_metricsBos_2$DPMVN$CPU_RUN_TIME)

tableBos_1 <- clustering_metricsBos_2$True$Table
tableBos_2 <- dp_Gmvlg_metricsBos_1$dp_Gmvlg$Table
tableBos_3 <- clustering_metricsBos_2$KMeans$Table
tableBos_4 <- clustering_metricsBos_2$CLARA$Table
tableBos_5 <- clustering_metricsBos_2$PAM$Table
tableBos_6 <- clustering_metricsBos_2$Hierarchical$Table
tableBos_7 <- clustering_metricsBos_2$Mclust$Table
tableBos_8 <- clustering_metricsBos_2$DPMVN$Table

colnames(row_Bos_K) <- NULL
colnames(row_Bos_ARI) <- NULL
colnames(row_Bos_NMI) <- NULL
colnames(row_Bos_CPU) <- NULL

row_Bos_K <- as.matrix(row_Bos_K)
row_Bos_ARI <- as.matrix(row_Bos_ARI)
row_Bos_NMI <- as.matrix(row_Bos_NMI)
row_Bos_CPU <- as.matrix(row_Bos_CPU)

kappa_table_Bos <- rbind(row_Bos_K)
ARI_table_Bos <- rbind(row_Bos_ARI)
NMI_table_Bos <- rbind(row_Bos_NMI)
CPU_runtime_table_Bos <- rbind(row_Bos_CPU)
clustering_table_Bos <- rbind(tableBos_1,tableBos_2,tableBos_3,tableBos_4,tableBos_5,tableBos_6,tableBos_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Bos)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Bos)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Bos)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", CPU_runtime_table_Bos)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Bos)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/19_BOSTON HOUSING/Cluster_metricsBos.xlsx", overwrite = TRUE)
