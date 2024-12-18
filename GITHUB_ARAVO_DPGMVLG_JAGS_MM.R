#################################################################
#                       ARAVO DATASET                          ##
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
library(moments)
library(aricode)

# Load the Aravo dataset
data(aravo,package = "ade4")
Aravo_data <- as.matrix(aravo$env[,-c(3,5)]) # Convert all columns to numeric
true_lab_Aravo <- as.numeric(aravo$env$ZoogD)    # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Aravo_data)

# 2. Check for missing values
sum(is.na(Aravo_data))

# 3. Distribution of each feature
par(mfrow=c(2, 3))  # Set up the plotting area for 60 histograms
for (i in 1:ncol(Aravo_data)) {
  hist(Aravo_data[, i], main=colnames(Aravo_data)[i], xlab=colnames(Aravo_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Aravo_data, main="ARAVO Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Aravo])

# 5. Boxplots for each feature by 'Class'
# Add normalized data to the original data frame
sonar_normalized <- as.data.frame(Aravo_data)
sonar_normalized$Class <- aravo$env$ZoogD

par(mfrow=c(2, 3))  # Reset the plotting area for boxplots
for (i in 1:ncol(Aravo_data)) {
  boxplot(Aravo_data[, i] ~ aravo$env$ZoogD, main=colnames(Aravo_data)[i], xlab="Class", ylab=colnames(Aravo_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Aravo_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix[,], method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Aravo data
Aravo_scaled_data <- scale(Aravo_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Aravo_scaled_data, 2, skewness)
kurtosis_values <- apply(Aravo_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Aravo_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Aravo_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Aravo_scaled_data[,i])$out
  print(paste("Feature:", colnames(Aravo_scaled_data)[i]))
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
outliers_matrix <- apply(Aravo_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Aravo_scaled_data, cluster = true_lab_Aravo), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Aravo_scaled_data)
D <- ncol(Aravo_scaled_data)

#Try with different number of clusters
G <- 3  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Aravo_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Aravo_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Aravo_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(4,3)          
       lambda_g[j, g] ~ dgamma(1,6) 
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
run_time_Aravo <- toc()
codaSamples_Aravo = as.mcmc.list( runJagsOut )
summaryChains_Aravo <- summary(codaSamples_Aravo)

diagMCMC( codaObject=codaSamples_Aravo , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Aravo , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Aravo$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_Aravo$statistics),1], P, G)
z_mode_Aravo <- apply(matrix(summaryChains_Aravo$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_Aravo$statistics),1], P, G),1, which.max)
z_mode_Aravo

plot(Aravo_scaled_data, col= z_mode_Aravo, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Aravo), col=unique(z_mode_Aravo), pch=16, title="Cluster")
table( true_lab_Aravo , z_mode_Aravo)
# To switch to the same labels in the true clustering
new <- c(1,2)
old <- c(2,1)
z_mode_Aravo[z_mode_Aravo %in% old] <- new[match(z_mode_Aravo,old,nomatch = 0)]

plot(Aravo_scaled_data, col= z_mode_Aravo, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Aravo), col=unique(z_mode_Aravo), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Aravo , z_mode_Aravo)
kappa2(data.frame(rater1 = true_lab_Aravo, rater2 = z_mode_Aravo))

calculate_dp_Gmvlg_clustering_metrics <- function(dataAravo_1, true_clusters_Aravo) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Aravo <- table(true_clusters_Aravo, z_mode_Aravo)
  kappa_result_Aravo <- kappa2(data.frame(rater1 = true_clusters_Aravo, rater2 = z_mode_Aravo))
  ari_result_Aravo <- adjustedRandIndex(true_clusters_Aravo, z_mode_Aravo)
  nmi_result_Aravo <- NMI(true_clusters_Aravo, z_mode_Aravo)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Aravo,
    Kappa = kappa_result_Aravo$value,
    ARI = ari_result_Aravo,
    NMI = nmi_result_Aravo,
    CPU_RUN_TIME = run_time_Aravo$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsAravo_1 <- calculate_dp_Gmvlg_clustering_metrics(dataAravo_1 = Aravo_scaled_data, true_clusters_Aravo = true_lab_Aravo)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsAravo_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsAravo_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsAravo_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsAravo_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsAravo_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsAravo_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Aravo <- eclust(Aravo_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Aravo <- toc()
kmeans_clusters_Aravo <- kmeans_result_Aravo$cluster
table(true_lab_Aravo, kmeans_clusters_Aravo)
new <- c(1,3)
old <- c(3,1)
kmeans_clusters_Aravo[kmeans_clusters_Aravo %in% old] <- new[match(kmeans_clusters_Aravo, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_Aravo[kmeans_clusters_Aravo %in% old] <- new[match(kmeans_clusters_Aravo, old, nomatch = 0)]
table(true_lab_Aravo, kmeans_clusters_Aravo)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Aravo <- eclust(Aravo_scaled_data, "clara", G, graph = FALSE)
clara_time_Aravo <- toc()
clara_clusters_Aravo <- clara_result_Aravo$cluster
table(true_lab_Aravo, clara_clusters_Aravo)
new <- c(1,3)
old <- c(3,1)
clara_clusters_Aravo[clara_clusters_Aravo %in% old] <- new[match(clara_clusters_Aravo, old, nomatch = 0)]
table(true_lab_Aravo, clara_clusters_Aravo)

# PAM clustering
tic("PAM Runtime")
pam_result_Aravo <- eclust(Aravo_scaled_data, "pam", G, graph = FALSE)
pam_time_Aravo <- toc()
pam_clusters_Aravo <- pam_result_Aravo$cluster
table(true_lab_Aravo, pam_clusters_Aravo)
new <- c(1,3)
old <- c(3,1)
pam_clusters_Aravo[pam_clusters_Aravo %in% old] <- new[match(pam_clusters_Aravo, old, nomatch = 0)]
table(true_lab_Aravo, pam_clusters_Aravo)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Aravo <- hclust(dist(Aravo_scaled_data), method = "ward.D2")
hclust_time_Aravo <- toc()
hclust_clusters_Aravo <- cutree(hclust_result_Aravo, k = G)
table(true_lab_Aravo, hclust_clusters_Aravo)
new <- c(1,3)
old <- c(3,1)
hclust_clusters_Aravo[hclust_clusters_Aravo %in% old] <- new[match(hclust_clusters_Aravo, old, nomatch = 0)]
table(true_lab_Aravo, hclust_clusters_Aravo)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Aravo <- Mclust(Aravo_scaled_data, G = G)
mclust_time_Aravo <- toc()
summary(mclust_result_Aravo)
mclust_clusters_Aravo <- mclust_result_Aravo$classification
table(true_lab_Aravo, mclust_clusters_Aravo)
new <- 1:2
old <- c(2,1)
mclust_clusters_Aravo[mclust_clusters_Aravo %in% old] <- new[match(mclust_clusters_Aravo, old, nomatch = 0)]
table(true_lab_Aravo, mclust_clusters_Aravo)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Aravo_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Aravo <- toc()
# Extract clusters 
dpMVN_clusters_Aravo <- as.numeric(dp$clusterLabels)
#new <- 1:2
#old <- c(2,1)
#dpMVN_clusters_Aravo[dpMVN_clusters_Aravo %in% old] <- new[match(dpMVN_clusters_Aravo, old, nomatch = 0)]
print(dpMVN_clusters_Aravo)
table(true_lab_Aravo, dpMVN_clusters_Aravo)

calculate_clustering_metricsAravo_2 <- function(dataAravo_2, true_clusters_Aravo, estimated_clusters_Aravo_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Aravo_list)) {
    clusters <- estimated_clusters_Aravo_list[[method_name]]
    
    # Calculate metrics
    table_result_Aravo <- table(true_clusters_Aravo, clusters)
    kappa_result_Aravo <- kappa2(data.frame(rater1 = true_clusters_Aravo, rater2 = clusters))
    ari_result_Aravo <- adjustedRandIndex(true_clusters_Aravo, clusters)
    nmi_result_Aravo <- NMI(true_clusters_Aravo, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Aravo,
      Kappa = kappa_result_Aravo$value,
      ARI = ari_result_Aravo,
      NMI = nmi_result_Aravo,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Aravo <- list(
  KMeans = kmeans_clusters_Aravo,
  CLARA = clara_clusters_Aravo,
  PAM = pam_clusters_Aravo,
  Hierarchical = hclust_clusters_Aravo,
  Mclust = mclust_clusters_Aravo,
  DPMVN = dpMVN_clusters_Aravo,
  True = true_lab_Aravo
)

times_list_Aravo <- list(
  KMeans = kmeans_time_Aravo,
  CLARA = clara_time_Aravo,
  PAM = pam_time_Aravo,
  Hierarchical = hclust_time_Aravo,
  Mclust = mclust_time_Aravo,
  DPMVN = DPMVN_time_Aravo
)

# Call the function
clustering_metricsAravo_2 <- calculate_clustering_metricsAravo_2(dataAravo_2 = Aravo_scaled_data, true_clusters_Aravo = true_lab_Aravo, estimated_clusters_Aravo_list = cluster_result_Aravo, times_list = times_list_Aravo)

# Print the results for each method
for (method_name in names(clustering_metricsAravo_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsAravo_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsAravo_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsAravo_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsAravo_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsAravo_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Aravo_K <- data.frame(clustering_metricsAravo_2$True$Kappa, dp_Gmvlg_metricsAravo_1$dp_Gmvlg$Kappa, clustering_metricsAravo_2$KMeans$Kappa, clustering_metricsAravo_2$CLARA$Kappa, clustering_metricsAravo_2$PAM$Kappa, clustering_metricsAravo_2$Hierarchical$Kappa, clustering_metricsAravo_2$Mclust$Kappa, clustering_metricsAravo_2$DPMVN$Kappa)
row_Aravo_ARI <- data.frame(clustering_metricsAravo_2$True$ARI, dp_Gmvlg_metricsAravo_1$dp_Gmvlg$ARI, clustering_metricsAravo_2$KMeans$ARI, clustering_metricsAravo_2$CLARA$ARI, clustering_metricsAravo_2$PAM$ARI, clustering_metricsAravo_2$Hierarchical$ARI, clustering_metricsAravo_2$Mclust$ARI, clustering_metricsAravo_2$DPMVN$ARI)
row_Aravo_NMI <- data.frame(clustering_metricsAravo_2$True$NMI, dp_Gmvlg_metricsAravo_1$dp_Gmvlg$NMI, clustering_metricsAravo_2$KMeans$NMI, clustering_metricsAravo_2$CLARA$NMI, clustering_metricsAravo_2$PAM$NMI, clustering_metricsAravo_2$Hierarchical$NMI, clustering_metricsAravo_2$Mclust$NMI, clustering_metricsAravo_2$DPMVN$NMI)
row_Aravo_CPU <- data.frame(dp_Gmvlg_metricsAravo_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsAravo_2$KMeans$CPU_RUN_TIME, clustering_metricsAravo_2$CLARA$CPU_RUN_TIME, clustering_metricsAravo_2$PAM$CPU_RUN_TIME, clustering_metricsAravo_2$Hierarchical$CPU_RUN_TIME, clustering_metricsAravo_2$Mclust$CPU_RUN_TIME, clustering_metricsAravo_2$DPMVN$CPU_RUN_TIME)

tableAravo_1 <- clustering_metricsAravo_2$True$Table
tableAravo_2 <- dp_Gmvlg_metricsAravo_1$dp_Gmvlg$Table
tableAravo_3 <- clustering_metricsAravo_2$KMeans$Table
tableAravo_4 <- clustering_metricsAravo_2$CLARA$Table
tableAravo_5 <- clustering_metricsAravo_2$PAM$Table
tableAravo_6 <- clustering_metricsAravo_2$Hierarchical$Table
tableAravo_7 <- clustering_metricsAravo_2$Mclust$Table
tableAravo_8 <- clustering_metricsAravo_2$DPMVN$Table

colnames(row_Aravo_K) <- NULL
colnames(row_Aravo_ARI) <- NULL
colnames(row_Aravo_NMI) <- NULL
colnames(row_Aravo_CPU) <- NULL

row_Aravo_K <- as.matrix(row_Aravo_K)
row_Aravo_ARI <- as.matrix(row_Aravo_ARI)
row_Aravo_NMI <- as.matrix(row_Aravo_NMI)
row_Aravo_CPU <- as.matrix(row_Aravo_CPU)

kappa_table_Aravo <- rbind(row_Aravo_K)
ARI_table_Aravo <- rbind(row_Aravo_ARI)
NMI_table_Aravo <- rbind(row_Aravo_NMI)
CPU_runtime_table_Aravo <- rbind(row_Aravo_CPU)
clustering_table_Aravo <- rbind(tableAravo_1,tableAravo_2,tableAravo_3,tableAravo_4,tableAravo_5,tableAravo_6,tableAravo_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Aravo)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Aravo)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Aravo)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", CPU_runtime_table_Aravo)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Aravo)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/11_ARAVO/Cluster_metricsARAVO.xlsx", overwrite = TRUE)

