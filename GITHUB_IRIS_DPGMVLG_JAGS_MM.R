#####################################
#             IRIS DATASET         ##
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
library(aricode)

# Load the Iris dataset
data(iris)
iris_data <- as.matrix(iris[, c(1,2,4)])
true_lab_iris=iris[,5]
true_lab_iris=as.numeric(true_lab_iris)

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(iris_data)

# 2. Check for missing values
sum(is.na(iris_data))

# 3. Distribution of each feature
par(mfrow=c(2,2))  # Set up the plotting area
hist(iris_data[,1], main="Sepal Length", xlab="Sepal Length", col="lightblue", border="black")
hist(iris_data[,2], main="Sepal Width", xlab="Sepal Width", col="lightblue", border="black")
#hist(iris_data[,3], main="Petal Length", xlab="Petal Length", col="lightblue", border="black")
hist(iris_data[,3], main="Petal Width", xlab="Petal Width", col="lightblue", border="black")

# 4. Pairwise scatter plots
pairs(iris_data, main="Iris Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[unclass(iris$Species)])

# 5. Boxplots for each feature by species
# Add normalized data to the original data frame
iris_normalized <- as.data.frame(iris_data)
iris_normalized$Species <- iris$Species

par(mfrow=c(2,2))  # Reset the plotting area
boxplot(Sepal.Length ~ Species, data=iris_normalized, main="Sepal Length by Species", xlab="Species", ylab="Sepal Length", col="lightblue")
boxplot(Sepal.Width ~ Species, data=iris_normalized, main="Sepal Width by Species", xlab="Species", ylab="Sepal Width", col="lightblue")
#boxplot(Petal.Length ~ Species, data=iris_normalized, main="Petal Length by Species", xlab="Species", ylab="Petal Length", col="lightblue")
boxplot(Petal.Width ~ Species, data=iris_normalized, main="Petal Width by Species", xlab="Species", ylab="Petal Width", col="lightblue")

# 6. Correlation matrix
cor_matrix <- cor(iris_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1,1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the iris data
iris_scaled_data <- scale(iris_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(iris_scaled_data, 2, skewness)
kurtosis_values <- apply(iris_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

# Plot skewness values
par(mfrow=c(1,2))
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=1)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=1)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(iris_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(iris_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(iris_scaled_data[,i])$out
  print(paste("Feature:", colnames(iris_scaled_data)[i]))
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
outliers_matrix <- apply(iris_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = iris_scaled_data, cluster = true_lab_iris), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(iris_scaled_data)
D <- ncol(iris_scaled_data)

G = 3 # maximum number of clusters

# Set up data for JAGS
data_list <- list(
  x = iris_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(iris_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(iris_scaled_data))^(1 / (D - 1))

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
    alpha[g] ~ dbeta(2,4)
   }
  pi_g[1:G] ~ ddirich(alpha/sum(alpha))
}
"
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
                        burnin=2000 ,
                        inits=inits,
                        sample=1000,
                        thin = 2) 

run_time_iris <- toc() # mins (/60)
codaSamples_iris = as.mcmc.list( runJagsOut )
summaryChains_iris <- summary(codaSamples_iris)

diagMCMC( codaObject=codaSamples_iris , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_iris , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_iris , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_iris , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_iris , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_iris , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_iris , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_iris , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_iris$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_iris$statistics),1], P, G) 
z_mode_iris <- apply(matrix(summaryChains_iris$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_iris$statistics),1], P, G),1, which.max)
z_mode_iris

plot(iris_scaled_data, col= z_mode_iris, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_iris), col=unique(z_mode_iris), pch=16, title="Cluster")
table( true_lab_iris , z_mode_iris)

# To switch to the same labels in the true clustering
new <- c(1,3) #the new vector contains the desired labels to match the true labels.
old <- c(3,1) #The old vector contains the cluster labels assigned by the clustering algorithm
z_mode_iris[z_mode_iris %in% old] <- new[match(z_mode_iris,old,nomatch = 0)]

# To switch to the same labels in the true clustering
new <- c(2,3) #the new vector contains the desired labels to match the true labels.
old <- c(3,2) #The old vector contains the cluster labels assigned by the clustering algorithm
z_mode_iris[z_mode_iris %in% old] <- new[match(z_mode_iris,old,nomatch = 0)]

plot(iris_scaled_data, col= z_mode_iris, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_iris), col=unique(z_mode_iris), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_iris , z_mode_iris)
kappa2(data.frame(rater1 = true_lab_iris, rater2 = z_mode_iris))

calculate_dp_Gmvlg_clustering_metrics <- function(datairis_1, true_clusters_iris) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_iris <- table(true_clusters_iris, z_mode_iris)
  kappa_result_iris <- kappa2(data.frame(rater1 = true_clusters_iris, rater2 = z_mode_iris))
  ari_result_iris <- adjustedRandIndex(true_clusters_iris, z_mode_iris)
  nmi_result_iris <- NMI(true_clusters_iris, z_mode_iris)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_iris,
    Kappa = kappa_result_iris$value,
    ARI = ari_result_iris,
    NMI = nmi_result_iris,
    CPU_RUN_TIME = run_time_iris$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsiris_1 <- calculate_dp_Gmvlg_clustering_metrics(datairis_1 = iris_scaled_data, true_clusters_iris = true_lab_iris)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsiris_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsiris_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsiris_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsiris_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsiris_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsiris_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_iris <- eclust(iris_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_iris <- toc()
kmeans_clusters_iris <- kmeans_result_iris$cluster
table(true_lab_iris, kmeans_clusters_iris)
new <- 1:2
old <- c(2,1)
kmeans_clusters_iris[kmeans_clusters_iris %in% old] <- new[match(kmeans_clusters_iris, old, nomatch = 0)]
table(true_lab_iris, kmeans_clusters_iris)

# CLARA clustering
tic("CLARA Runtime")
clara_result_iris <- eclust(iris_scaled_data, "clara", G, graph = FALSE)
clara_time_iris <- toc()
clara_clusters_iris <- clara_result_iris$cluster
table(true_lab_iris, clara_clusters_iris)

# PAM clustering
tic("PAM Runtime")
pam_result_iris <- eclust(iris_scaled_data, "pam", G, graph = FALSE)
pam_time_iris <- toc()
pam_clusters_iris <- pam_result_iris$cluster
table(true_lab_iris, pam_clusters_iris)
new <- c(2,3)
old <- c(3,2)
pam_clusters_iris[pam_clusters_iris %in% old] <- new[match(pam_clusters_iris, old, nomatch = 0)]
table(true_lab_iris, pam_clusters_iris)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_iris <- hclust(dist(iris_scaled_data), method = "ward.D2")
hclust_time_iris <- toc()
hclust_clusters_iris <- cutree(hclust_result_iris, k = G)
table(true_lab_iris, hclust_clusters_iris)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_iris <- Mclust(iris_scaled_data, G = G)
mclust_time_iris <- toc()
summary(mclust_result_iris)
mclust_clusters_iris <- mclust_result_iris$classification
table(true_lab_iris, mclust_clusters_iris)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(iris_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_iris <- toc()
dp
# Extract clusters 
dpMVN_clusters_iris <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_iris)
table(true_lab_iris, dpMVN_clusters_iris)

calculate_clustering_metricsiris_2 <- function(datairis_2, true_clusters_iris, estimated_clusters_iris_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_iris_list)) {
    clusters <- estimated_clusters_iris_list[[method_name]]
    
    # Calculate metrics
    table_result_iris <- table(true_clusters_iris, clusters)
    kappa_result_iris <- kappa2(data.frame(rater1 = true_clusters_iris, rater2 = clusters))
    ari_result_iris <- adjustedRandIndex(true_clusters_iris, clusters)
    nmi_result_iris <- NMI(true_clusters_iris, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_iris,
      Kappa = kappa_result_iris$value,
      ARI = ari_result_iris,
      NMI = nmi_result_iris,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_iris <- list(
  KMeans = kmeans_clusters_iris,
  CLARA = clara_clusters_iris,
  PAM = pam_clusters_iris,
  Hierarchical = hclust_clusters_iris,
  Mclust = mclust_clusters_iris,
  DPMVN = dpMVN_clusters_iris,
  True = true_lab_iris
)

times_list_iris <- list(
  KMeans = kmeans_time_iris,
  CLARA = clara_time_iris,
  PAM = pam_time_iris,
  Hierarchical = hclust_time_iris,
  Mclust = mclust_time_iris,
  DPMVN = DPMVN_time_iris
)

# Call the function
clustering_metricsiris_2 <- calculate_clustering_metricsiris_2(datairis_2 = iris_scaled_data, true_clusters_iris = true_lab_iris, estimated_clusters_iris_list = cluster_result_iris, times_list = times_list_iris)

# Print the results for each method
for (method_name in names(clustering_metricsiris_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsiris_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsiris_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsiris_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsiris_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsiris_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_iris_K <- data.frame(clustering_metricsiris_2$True$Kappa, dp_Gmvlg_metricsiris_1$dp_Gmvlg$Kappa, clustering_metricsiris_2$KMeans$Kappa, clustering_metricsiris_2$CLARA$Kappa, clustering_metricsiris_2$PAM$Kappa, clustering_metricsiris_2$Hierarchical$Kappa, clustering_metricsiris_2$Mclust$Kappa, clustering_metricsiris_2$DPMVN$Kappa)
row_iris_ARI <- data.frame(clustering_metricsiris_2$True$ARI, dp_Gmvlg_metricsiris_1$dp_Gmvlg$ARI, clustering_metricsiris_2$KMeans$ARI, clustering_metricsiris_2$CLARA$ARI, clustering_metricsiris_2$PAM$ARI, clustering_metricsiris_2$Hierarchical$ARI, clustering_metricsiris_2$Mclust$ARI, clustering_metricsiris_2$DPMVN$ARI)
row_iris_NMI <- data.frame(clustering_metricsiris_2$True$NMI, dp_Gmvlg_metricsiris_1$dp_Gmvlg$NMI, clustering_metricsiris_2$KMeans$NMI, clustering_metricsiris_2$CLARA$NMI, clustering_metricsiris_2$PAM$NMI, clustering_metricsiris_2$Hierarchical$NMI, clustering_metricsiris_2$Mclust$NMI, clustering_metricsiris_2$DPMVN$NMI)
row_iris_Cpu <- data.frame(dp_Gmvlg_metricsiris_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsiris_2$KMeans$CPU_RUN_TIME, clustering_metricsiris_2$CLARA$CPU_RUN_TIME, clustering_metricsiris_2$PAM$CPU_RUN_TIME, clustering_metricsiris_2$Hierarchical$CPU_RUN_TIME, clustering_metricsiris_2$Mclust$CPU_RUN_TIME, clustering_metricsiris_2$DPMVN$CPU_RUN_TIME)

tableiris_1 <- clustering_metricsiris_2$True$Table
tableiris_2 <- dp_Gmvlg_metricsiris_1$dp_Gmvlg$Table
tableiris_3 <- clustering_metricsiris_2$KMeans$Table
tableiris_4 <- clustering_metricsiris_2$CLARA$Table
tableiris_5 <- clustering_metricsiris_2$PAM$Table
tableiris_6 <- clustering_metricsiris_2$Hierarchical$Table
tableiris_7 <- clustering_metricsiris_2$Mclust$Table
tableiris_8 <- clustering_metricsiris_2$DPMVN$Table

colnames(row_iris_K) <- NULL
colnames(row_iris_ARI) <- NULL
colnames(row_iris_NMI) <- NULL
colnames(row_iris_Cpu) <- NULL

row_iris_K <- as.matrix(row_iris_K)
row_iris_ARI <- as.matrix(row_iris_ARI)
row_iris_NMI <- as.matrix(row_iris_NMI)
row_iris_Cpu <- as.matrix(row_iris_Cpu)

kappa_table_iris <- rbind(row_iris_K)
ARI_table_iris <- rbind(row_iris_ARI)
NMI_table_iris <- rbind(row_iris_NMI)
cpu_runtime_table_iris <- rbind(row_iris_Cpu)
clustering_table_iris <- rbind(tableiris_1,tableiris_2,tableiris_3,tableiris_4,tableiris_5,tableiris_6,tableiris_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_iris)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_iris)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_iris)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_iris)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_iris)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/14_IRIS/Cluster_metricsiris.xlsx", overwrite = TRUE)
