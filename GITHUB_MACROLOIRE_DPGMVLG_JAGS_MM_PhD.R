#################################################################
#                       MACROLOIRE DATASET                     ##
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
data("macroloire",package ="ade4")
help(macroloire) 
macro_data <- as.matrix(macroloire$envir[, c(2,3)])
true_lab_macro=as.numeric(macroloire$envir[,5])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(macro_data)

# 2. Check for missing values
sum(is.na(macro_data))

# 3. Distribution of each feature
par(mfrow=c(4, 3))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(macro_data)) {
  hist(macro_data[, i], main=colnames(macro_data)[i], xlab=colnames(macro_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(macro_data, main="Corvus Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_macro])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(macro_data)
Cor_normalized$Order <- macroloire$taxo$Order

par(mfrow=c(4, 3))  # Reset the plotting area for boxplots
for (i in 1:ncol(macro_data)) {
  boxplot(macro_data[, i] ~ Cor_normalized$Order, main=colnames(macro_data)[i], xlab="Habitat", ylab=colnames(macro_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(macro_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Corvus dataset
macro_scaled_data <- scale(macro_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(macro_scaled_data, 2, skewness)
kurtosis_values <- apply(macro_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=1)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=1)

#outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(macro_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(macro_scaled_data[,i])$out
  print(paste("Feature:", colnames(macro_scaled_data)[i]))
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
outliers_matrix <- apply(macro_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = macro_scaled_data, cluster = true_lab_macro), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(macro_scaled_data)
D <- ncol(macro_scaled_data)

#Try with different number of clusters
G <- 3  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = macro_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(macro_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(macro_scaled_data))^(1 / (D - 1))

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
run_time_macro <- toc()
codaSamples_macro = as.mcmc.list( runJagsOut )
summaryChains_macro <- summary(codaSamples_macro)

diagMCMC( codaObject=codaSamples_macro , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_macro , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_macro , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_macro , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_macro , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_macro , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_macro , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_macro , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_macro$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_macro$statistics),1], P, G)
z_mode_macro <- apply(matrix(summaryChains_macro$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_macro$statistics),1], P, G),1, which.max)
z_mode_macro

plot(macro_scaled_data, col= z_mode_macro, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_macro), col=unique(z_mode_macro), pch=16, title="Cluster")
table(true_lab_macro , z_mode_macro)

# To switch to the same labels in the true clustering
new <- c(2,3)
old <- c(3,2)
z_mode_macro[z_mode_macro %in% old] <- new[match(z_mode_macro,old,nomatch = 0)]

plot(macro_scaled_data, col= z_mode_macro, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_macro), col=unique(z_mode_macro), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table(true_lab_macro , z_mode_macro)
kappa2(data.frame(rater1 = true_lab_macro, rater2 = z_mode_macro))

calculate_dp_Gmvlg_clustering_metrics <- function(datamacro_1, true_clusters_macro) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_macro <- table(true_clusters_macro, z_mode_macro)
  kappa_result_macro <- kappa2(data.frame(rater1 = true_clusters_macro, rater2 = z_mode_macro))
  ari_result_macro <- adjustedRandIndex(true_clusters_macro, z_mode_macro)
  nmi_result_macro <- NMI(true_clusters_macro, z_mode_macro)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_macro,
    Kappa = kappa_result_macro$value,
    ARI = ari_result_macro,
    NMI = nmi_result_macro,
    CPU_RUN_TIME = run_time_macro$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsmacro_1 <- calculate_dp_Gmvlg_clustering_metrics(datamacro_1 = macro_scaled_data, true_clusters_macro = true_lab_macro)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsmacro_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsmacro_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsmacro_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsmacro_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsmacro_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsmacro_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_macro <- eclust(macro_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_macro <- toc()
kmeans_clusters_macro <- kmeans_result_macro$cluster
table( true_lab_macro , kmeans_clusters_macro)
new <- c(1,3)
old <- c(3,1)
kmeans_clusters_macro[kmeans_clusters_macro %in% old] <- new[match(kmeans_clusters_macro,old,nomatch = 0)]
table( true_lab_macro , kmeans_clusters_macro)

# CLARA clustering
tic("CLARA Runtime")
clara_result_macro <- eclust(macro_scaled_data, "clara", G, graph = FALSE)
clara_time_macro <- toc()
clara_clusters_macro <- clara_result_macro$cluster
table( true_lab_macro , clara_clusters_macro)

# PAM clustering
tic("PAM Runtime")
pam_result_macro <- eclust(macro_scaled_data, "pam", G, graph = FALSE)
pam_time_macro <- toc()
pam_clusters_macro <- pam_result_macro$cluster
table( true_lab_macro , pam_clusters_macro)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_macro <- hclust(dist(macro_scaled_data), method = "ward.D2")
hclust_time_macro <- toc()
hclust_clusters_macro <- cutree(hclust_result_macro, k = G)
table( true_lab_macro , hclust_clusters_macro)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_macro <- Mclust(macro_scaled_data, G = G)
mclust_time_macro <- toc()
summary(mclust_result_macro)
mclust_clusters_macro <- mclust_result_macro$classification
table( true_lab_macro , mclust_clusters_macro)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(macro_scaled_data, alphaPriors = c(25,20))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_macro <- toc()
dp
# Extract clusters 
dpMVN_clusters_macro <- as.numeric(dp$clusterLabels)
new <- c(1,3)
old <- c(3,1)
dpMVN_clusters_macro[dpMVN_clusters_macro %in% old] <- new[match(dpMVN_clusters_macro, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
dpMVN_clusters_macro[dpMVN_clusters_macro %in% old] <- new[match(dpMVN_clusters_macro, old, nomatch = 0)]
print(dpMVN_clusters_macro)
table(true_lab_macro, dpMVN_clusters_macro)


calculate_clustering_metricsmacro_2 <- function(datamacro_2, true_clusters_macro, estimated_clusters_macro_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_macro_list)) {
    clusters <- estimated_clusters_macro_list[[method_name]]
    
    # Calculate metrics
    table_result_macro <- table(true_clusters_macro, clusters)
    kappa_result_macro <- kappa2(data.frame(rater1 = true_clusters_macro, rater2 = clusters))
    ari_result_macro <- adjustedRandIndex(true_clusters_macro, clusters)
    nmi_result_macro <- NMI(true_clusters_macro, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_macro,
      Kappa = kappa_result_macro$value,
      ARI = ari_result_macro,
      NMI = nmi_result_macro,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_macro <- list(
  KMeans = kmeans_clusters_macro,
  CLARA = clara_clusters_macro,
  PAM = pam_clusters_macro,
  Hierarchical = hclust_clusters_macro,
  Mclust = mclust_clusters_macro,
  DPMVN = dpMVN_clusters_macro,
  True = true_lab_macro
)

times_list_macro <- list(
  KMeans = kmeans_time_macro,
  CLARA = clara_time_macro,
  PAM = pam_time_macro,
  Hierarchical = hclust_time_macro,
  Mclust = mclust_time_macro,
  DPMVN = DPMVN_time_macro
)

# Call the function
clustering_metricsmacro_2 <- calculate_clustering_metricsmacro_2(datamacro_2 = macro_scaled_data, true_clusters_macro = true_lab_macro, estimated_clusters_macro_list = cluster_result_macro, times_list = times_list_macro)

# Print the results for each method
for (method_name in names(clustering_metricsmacro_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsmacro_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsmacro_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsmacro_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsmacro_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsmacro_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_macro_K <- data.frame(clustering_metricsmacro_2$True$Kappa, dp_Gmvlg_metricsmacro_1$dp_Gmvlg$Kappa, clustering_metricsmacro_2$KMeans$Kappa, clustering_metricsmacro_2$CLARA$Kappa, clustering_metricsmacro_2$PAM$Kappa, clustering_metricsmacro_2$Hierarchical$Kappa, clustering_metricsmacro_2$Mclust$Kappa, clustering_metricsmacro_2$DPMVN$Kappa)
row_macro_ARI <- data.frame(clustering_metricsmacro_2$True$ARI, dp_Gmvlg_metricsmacro_1$dp_Gmvlg$ARI, clustering_metricsmacro_2$KMeans$ARI, clustering_metricsmacro_2$CLARA$ARI, clustering_metricsmacro_2$PAM$ARI, clustering_metricsmacro_2$Hierarchical$ARI, clustering_metricsmacro_2$Mclust$ARI, clustering_metricsmacro_2$DPMVN$ARI)
row_macro_NMI <- data.frame(clustering_metricsmacro_2$True$NMI, dp_Gmvlg_metricsmacro_1$dp_Gmvlg$NMI, clustering_metricsmacro_2$KMeans$NMI, clustering_metricsmacro_2$CLARA$NMI, clustering_metricsmacro_2$PAM$NMI, clustering_metricsmacro_2$Hierarchical$NMI, clustering_metricsmacro_2$Mclust$NMI, clustering_metricsmacro_2$DPMVN$NMI)
row_macro_Cpu <- data.frame(dp_Gmvlg_metricsmacro_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsmacro_2$KMeans$CPU_RUN_TIME, clustering_metricsmacro_2$CLARA$CPU_RUN_TIME, clustering_metricsmacro_2$PAM$CPU_RUN_TIME, clustering_metricsmacro_2$Hierarchical$CPU_RUN_TIME, clustering_metricsmacro_2$Mclust$CPU_RUN_TIME, clustering_metricsmacro_2$DPMVN$CPU_RUN_TIME)

tablemacro_1 <- clustering_metricsmacro_2$True$Table
tablemacro_2 <- dp_Gmvlg_metricsmacro_1$dp_Gmvlg$Table
tablemacro_3 <- clustering_metricsmacro_2$KMeans$Table
tablemacro_4 <- clustering_metricsmacro_2$CLARA$Table
tablemacro_5 <- clustering_metricsmacro_2$PAM$Table
tablemacro_6 <- clustering_metricsmacro_2$Hierarchical$Table
tablemacro_7 <- clustering_metricsmacro_2$Mclust$Table
tablemacro_8 <- clustering_metricsmacro_2$DPMVN$Table

colnames(row_macro_K) <- NULL
colnames(row_macro_ARI) <- NULL
colnames(row_macro_NMI) <- NULL
colnames(row_macro_Cpu) <- NULL

row_macro_K <- as.matrix(row_macro_K)
row_macro_ARI <- as.matrix(row_macro_ARI)
row_macro_NMI <- as.matrix(row_macro_NMI)
row_macro_Cpu <- as.matrix(row_macro_Cpu)

kappa_table_macro <- rbind(row_macro_K)
ARI_table_macro <- rbind(row_macro_ARI)
NMI_table_macro <- rbind(row_macro_NMI)
cpu_runtime_table_macro <- rbind(row_macro_Cpu)
clustering_table_macro <- rbind(tablemacro_1,tablemacro_2,tablemacro_3,tablemacro_4,tablemacro_5,tablemacro_6,tablemacro_7, tablemacro_8)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_macro)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_macro)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_macro)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_macro)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_macro)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/07_MACROLOIRE/Cluster_metricsmacro.xlsx", overwrite = TRUE)
