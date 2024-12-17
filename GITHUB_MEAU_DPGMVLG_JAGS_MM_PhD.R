#################################################################
#                       MEAU DATASET                           ##
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

# Load the dataset
data("meau", package ="ade4")
help(meau)
meau1 <- data.frame(meau$env,meau$design$season)
names(meau1)[ncol(meau1)] <- "Season"
Meau_data <- as.matrix(meau1[, -c(7,8,11)])
true_lab_Meau=as.numeric(meau1[,11])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Meau_data)

# 2. Check for missing values
sum(is.na(Meau_data))

# 3. Distribution of each feature
par(mfrow=c(3, 3))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(Meau_data)) {
  hist(Meau_data[, i], main=colnames(Meau_data)[i], xlab=colnames(Meau_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Meau_data, main="Meau Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_Meau])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(Meau_data)
Cor_normalized$Season <- meau$design$season

par(mfrow=c(3,3))  # Reset the plotting area for boxplots
for (i in 1:ncol(Meau_data)) {
  boxplot(Meau_data[, i] ~ Cor_normalized$Season, main=colnames(Meau_data)[i], xlab="Habitat", ylab=colnames(Meau_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Meau_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the corvus dataset
Meau_scaled_data <- scale(Meau_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Meau_scaled_data, 2, skewness)
kurtosis_values <- apply(Meau_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Meau_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Meau_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Meau_scaled_data[,i])$out
  print(paste("Feature:", colnames(Meau_scaled_data)[i]))
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
outliers_matrix <- apply(Meau_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Meau_scaled_data, cluster = true_lab_Meau), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Meau_scaled_data)
D <- ncol(Meau_scaled_data)

#Try with different number of clusters
G <- 4  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Meau_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Meau_scaled_data))^(1 / (D - 1)),
  G = G
)

deltaData <- det(cor(Meau_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(4,1)         
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
run_time_Meau <- toc()
codaSamples_Meau = as.mcmc.list( runJagsOut )
summaryChains_Meau <- summary(codaSamples_Meau)

diagMCMC( codaObject=codaSamples_Meau , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Meau , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Meau , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Meau , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Meau , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Meau , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Meau , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Meau , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Meau$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_Meau$statistics),1], P, G)
z_mode_Meau <- apply(matrix(summaryChains_Meau$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_Meau$statistics),1], P, G),1, which.max)
z_mode_Meau

plot(Meau_scaled_data, col= z_mode_Meau, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Meau), col=unique(z_mode_Meau), pch=16, title="Cluster")
table( true_lab_Meau , z_mode_Meau)

# To switch to the same labels in the true clustering
new <- c(1,2)
old <- c(2,1)
z_mode_Meau[z_mode_Meau %in% old] <- new[match(z_mode_Meau,old,nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
z_mode_Meau[z_mode_Meau %in% old] <- new[match(z_mode_Meau,old,nomatch = 0)]
#new <- c(3,4)
#old <- c(4,3)
#z_mode_Meau[z_mode_Meau %in% old] <- new[match(z_mode_Meau,old,nomatch = 0)]
#new <- c(2,4)
#old <- c(4,2)
#z_mode_Meau[z_mode_Meau %in% old] <- new[match(z_mode_Meau,old,nomatch = 0)]

plot(Meau_scaled_data, col= z_mode_Meau, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Meau), col=unique(z_mode_Meau), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table(true_lab_Meau , z_mode_Meau)
kappa2(data.frame(rater1 = true_lab_Meau, rater2 = z_mode_Meau))

calculate_dp_Gmvlg_clustering_metrics <- function(dataMeau_1, true_clusters_Meau) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Meau <- table(true_clusters_Meau, z_mode_Meau)
  kappa_result_Meau <- kappa2(data.frame(rater1 = true_clusters_Meau, rater2 = z_mode_Meau))
  ari_result_Meau <- adjustedRandIndex(true_clusters_Meau, z_mode_Meau)
  nmi_result_Meau <- NMI(true_clusters_Meau, z_mode_Meau)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Meau,
    Kappa = kappa_result_Meau$value,
    ARI = ari_result_Meau,
    NMI = nmi_result_Meau,
    CPU_RUN_TIME = run_time_Meau$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsMeau_1 <- calculate_dp_Gmvlg_clustering_metrics(dataMeau_1 = Meau_scaled_data, true_clusters_Meau = true_lab_Meau)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsMeau_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsMeau_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsMeau_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsMeau_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsMeau_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsMeau_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Meau <- eclust(Meau_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Meau <- toc()
kmeans_clusters_Meau <- kmeans_result_Meau$cluster
table(true_lab_Meau, kmeans_clusters_Meau)
new <- c(1,4)
old <- c(4,1)
kmeans_clusters_Meau[kmeans_clusters_Meau %in% old] <- new[match(kmeans_clusters_Meau,old,nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
kmeans_clusters_Meau[kmeans_clusters_Meau %in% old] <- new[match(kmeans_clusters_Meau,old,nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_Meau[kmeans_clusters_Meau %in% old] <- new[match(kmeans_clusters_Meau,old,nomatch = 0)]
table(true_lab_Meau, kmeans_clusters_Meau)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Meau <- eclust(Meau_scaled_data, "clara", G, graph = FALSE)
clara_time_Meau <- toc()
clara_clusters_Meau <- clara_result_Meau$cluster
table(true_lab_Meau, clara_clusters_Meau)

# PAM clustering
tic("PAM Runtime")
pam_result_Meau <- eclust(Meau_scaled_data, "pam", G, graph = FALSE)
pam_time_Meau <- toc()
pam_clusters_Meau <- pam_result_Meau$cluster
table(true_lab_Meau, pam_clusters_Meau)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Meau <- hclust(dist(Meau_scaled_data), method = "ward.D2")
hclust_time_Meau <- toc()
hclust_clusters_Meau <- cutree(hclust_result_Meau, k = G)
table(true_lab_Meau, hclust_clusters_Meau)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Meau <- Mclust(Meau_scaled_data, G = G)
mclust_time_Meau <- toc()
summary(mclust_result_Meau)
mclust_clusters_Meau <- mclust_result_Meau$classification
table(true_lab_Meau, mclust_clusters_Meau)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Meau_scaled_data, alphaPriors = c(6,19))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Meau <- toc()
# Extract clusters 
dpMVN_clusters_Meau <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Meau)
table(true_lab_Meau, dpMVN_clusters_Meau)
new <- c(1,2)
old <- c(2,1)
dpMVN_clusters_Meau[dpMVN_clusters_Meau %in% old] <- new[match(dpMVN_clusters_Meau, old, nomatch = 0)]
new <- c(1,3)
old <- c(3,1)
dpMVN_clusters_Meau[dpMVN_clusters_Meau %in% old] <- new[match(dpMVN_clusters_Meau, old, nomatch = 0)]
table(true_lab_Meau, dpMVN_clusters_Meau)


calculate_clustering_metricsMeau_2 <- function(dataMeau_2, true_clusters_Meau, estimated_clusters_Meau_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Meau_list)) {
    clusters <- estimated_clusters_Meau_list[[method_name]]
    
    # Calculate metrics
    table_result_Meau <- table(true_clusters_Meau, clusters)
    kappa_result_Meau <- kappa2(data.frame(rater1 = true_clusters_Meau, rater2 = clusters))
    ari_result_Meau <- adjustedRandIndex(true_clusters_Meau, clusters)
    nmi_result_Meau <- NMI(true_clusters_Meau, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Meau,
      Kappa = kappa_result_Meau$value,
      ARI = ari_result_Meau,
      NMI = nmi_result_Meau,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Meau <- list(
  KMeans = kmeans_clusters_Meau,
  CLARA = clara_clusters_Meau,
  PAM = pam_clusters_Meau,
  Hierarchical = hclust_clusters_Meau,
  Mclust = mclust_clusters_Meau,
  DPMVN = dpMVN_clusters_Meau,
  True = true_lab_Meau
)

times_list_Meau <- list(
  KMeans = kmeans_time_Meau,
  CLARA = clara_time_Meau,
  PAM = pam_time_Meau,
  Hierarchical = hclust_time_Meau,
  Mclust = mclust_time_Meau,
  DPMVN = DPMVN_time_Meau
)

# Call the function
clustering_metricsMeau_2 <- calculate_clustering_metricsMeau_2(dataMeau_2 = Meau_scaled_data, true_clusters_Meau = true_lab_Meau, estimated_clusters_Meau_list = cluster_result_Meau, times_list = times_list_Meau)

# Print the results for each method
for (method_name in names(clustering_metricsMeau_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsMeau_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsMeau_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsMeau_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsMeau_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsMeau_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Meau_K <- data.frame(clustering_metricsMeau_2$True$Kappa, dp_Gmvlg_metricsMeau_1$dp_Gmvlg$Kappa, clustering_metricsMeau_2$KMeans$Kappa, clustering_metricsMeau_2$CLARA$Kappa, clustering_metricsMeau_2$PAM$Kappa, clustering_metricsMeau_2$Hierarchical$Kappa, clustering_metricsMeau_2$Mclust$Kappa, clustering_metricsMeau_2$DPMVN$Kappa)
row_Meau_ARI <- data.frame(clustering_metricsMeau_2$True$ARI, dp_Gmvlg_metricsMeau_1$dp_Gmvlg$ARI, clustering_metricsMeau_2$KMeans$ARI, clustering_metricsMeau_2$CLARA$ARI, clustering_metricsMeau_2$PAM$ARI, clustering_metricsMeau_2$Hierarchical$ARI, clustering_metricsMeau_2$Mclust$ARI, clustering_metricsMeau_2$DPMVN$ARI)
row_Meau_NMI <- data.frame(clustering_metricsMeau_2$True$NMI, dp_Gmvlg_metricsMeau_1$dp_Gmvlg$NMI, clustering_metricsMeau_2$KMeans$NMI, clustering_metricsMeau_2$CLARA$NMI, clustering_metricsMeau_2$PAM$NMI, clustering_metricsMeau_2$Hierarchical$NMI, clustering_metricsMeau_2$Mclust$NMI, clustering_metricsMeau_2$DPMVN$NMI)
row_Meau_Cpu <- data.frame(dp_Gmvlg_metricsMeau_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsMeau_2$KMeans$CPU_RUN_TIME, clustering_metricsMeau_2$CLARA$CPU_RUN_TIME, clustering_metricsMeau_2$PAM$CPU_RUN_TIME, clustering_metricsMeau_2$Hierarchical$CPU_RUN_TIME, clustering_metricsMeau_2$Mclust$CPU_RUN_TIME, clustering_metricsMeau_2$DPMVN$CPU_RUN_TIME)

tableMeau_1 <- clustering_metricsMeau_2$True$Table
tableMeau_2 <- dp_Gmvlg_metricsMeau_1$dp_Gmvlg$Table
tableMeau_3 <- clustering_metricsMeau_2$KMeans$Table
tableMeau_4 <- clustering_metricsMeau_2$CLARA$Table
tableMeau_5 <- clustering_metricsMeau_2$PAM$Table
tableMeau_6 <- clustering_metricsMeau_2$Hierarchical$Table
tableMeau_7 <- clustering_metricsMeau_2$Mclust$Table
tableMeau_8 <- clustering_metricsMeau_2$DPMVN$Table

colnames(row_Meau_K) <- NULL
colnames(row_Meau_ARI) <- NULL
colnames(row_Meau_NMI) <- NULL
colnames(row_Meau_Cpu) <- NULL


row_Meau_K <- as.matrix(row_Meau_K)
row_Meau_ARI <- as.matrix(row_Meau_ARI)
row_Meau_NMI <- as.matrix(row_Meau_NMI)
row_Meau_Cpu <- as.matrix(row_Meau_Cpu)

kappa_table_Meau <- rbind(row_Meau_K)
ARI_table_Meau <- rbind(row_Meau_ARI)
NMI_table_Meau <- rbind(row_Meau_NMI)
cpu_runtime_table_Meau <- rbind(row_Meau_Cpu)
clustering_table_Meau <- rbind(tableMeau_1,tableMeau_2,tableMeau_3,tableMeau_4,tableMeau_5,tableMeau_6,tableMeau_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Meau)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Meau)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Meau)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Meau)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Meau)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/05_MEAU DATA/Cluster_metricsMeau.xlsx", overwrite = TRUE)

