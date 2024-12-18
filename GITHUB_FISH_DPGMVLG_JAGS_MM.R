##################################################
#                FISH CATCH DATASET             ##
##################################################

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
library(ade4)
library(moments)
library(aricode)

# Load the Fish Catch dataset
data("fish")
help(fish)
Fish_data <- as.matrix(fish[,c(3,5,6)])
true_lab_Fish=fish[,7]

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Fish_data)

# 2. Check for missing values
sum(is.na(Fish_data))

# 3. Distribution of each feature
par(mfrow=c(1, 3))  # Set up the plotting area for 3 histograms
for (i in 1:ncol(Fish_data)) {
  hist(Fish_data[, i], main=colnames(Fish_data)[i], xlab=colnames(Fish_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Fish_data, main="Fish Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Fish])

# 5. Boxplots for each feature by species
# Add normalized data to the original data frame
fish_normalized <- as.data.frame(Fish_data)
fish_normalized$Species <- fish$Species

par(mfrow=c(1, 3))  # Reset the plotting area for boxplots
for (i in 1:ncol(Fish_data)) {
  boxplot(Fish_data[, i] ~ fish$Species, main=colnames(Fish_data)[i], xlab="Species", ylab=colnames(Fish_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Fish_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Fish Catch data
Fish_scaled_data <- scale(Fish_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Fish_scaled_data, 2, skewness)
kurtosis_values <- apply(Fish_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Fish_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Fish_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Fish_scaled_data[,i])$out
  print(paste("Feature:", colnames(Fish_scaled_data)[i]))
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
outliers_matrix <- apply(Fish_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Fish_scaled_data, cluster = true_lab_Fish), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Fish_scaled_data)
D <- ncol(Fish_scaled_data)

#Try with different number of clusters
G <- 7  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Fish_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Fish_scaled_data))^(1 / (D - 1)),
  G = G
)

deltaData <- det(cor(Fish_scaled_data))^(1 / (D - 1))

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
       lambda_g[j, g] ~ dgamma(1,3) 
    }
     alpha[g] ~ dgamma(2,4)
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
run_time_Fish <- toc()
codaSamples_Fish = as.mcmc.list( runJagsOut )
summaryChains_Fish <- summary(codaSamples_Fish)

diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Fish$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_Fish$statistics),1], P, G)
z_mode_Fish <- apply(matrix(summaryChains_Fish$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_Fish$statistics),1], P, G),1, which.max)
z_mode_Fish

plot(Fish_scaled_data, col= z_mode_Fish, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Fish), col=unique(z_mode_Fish), pch=16, title="Cluster")
table( true_lab_Fish , z_mode_Fish)

# To switch to the same labels in the true clustering
new <- c(1,7)
old <- c(7,1)
z_mode_Fish[z_mode_Fish %in% old] <- new[match(z_mode_Fish,old,nomatch = 0)]
new <- c(2,5)
old <- c(5,2)
z_mode_Fish[z_mode_Fish %in% old] <- new[match(z_mode_Fish,old,nomatch = 0)]
new <- c(3,6)
old <- c(6,3)
z_mode_Fish[z_mode_Fish %in% old] <- new[match(z_mode_Fish,old,nomatch = 0)]
new <- c(5,7)
old <- c(7,5)
z_mode_Fish[z_mode_Fish %in% old] <- new[match(z_mode_Fish,old,nomatch = 0)]
new <- c(5,6)
old <- c(6,5)
z_mode_Fish[z_mode_Fish %in% old] <- new[match(z_mode_Fish,old,nomatch = 0)]

# To switch to the same labels in the true clustering
table( true_lab_Fish , z_mode_Fish)
kappa2(data.frame(rater1 = true_lab_Fish, rater2 = z_mode_Fish))

calculate_dp_Gmvlg_clustering_metrics <- function(dataFish_1, true_clusters_Fish) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Fish <- table(true_clusters_Fish, z_mode_Fish)
  kappa_result_Fish <- kappa2(data.frame(rater1 = true_clusters_Fish, rater2 = z_mode_Fish))
  ari_result_Fish <- adjustedRandIndex(true_clusters_Fish, z_mode_Fish)
  nmi_result_Fish <- NMI(true_clusters_Fish, z_mode_Fish)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Fish,
    Kappa = kappa_result_Fish$value,
    ARI = ari_result_Fish,
    NMI = nmi_result_Fish,
    CPU_RUN_TIME = run_time_Fish$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsFish_1 <- calculate_dp_Gmvlg_clustering_metrics(dataFish_1 = Fish_scaled_data, true_clusters_Fish = true_lab_Fish)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsFish_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Fish <- eclust(Fish_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Fish <- toc()
kmeans_clusters_Fish <- kmeans_result_Fish$cluster
table(true_lab_Fish, kmeans_clusters_Fish)
new <- 1:2
old <- c(2,1)
kmeans_clusters_Fish[kmeans_clusters_Fish %in% old] <- new[match(kmeans_clusters_Fish, old, nomatch = 0)]
new <- c(6,7)
old <- c(7,6)
kmeans_clusters_Fish[kmeans_clusters_Fish %in% old] <- new[match(kmeans_clusters_Fish, old, nomatch = 0)]
table(true_lab_Fish, kmeans_clusters_Fish)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Fish <- eclust(Fish_scaled_data, "clara", G, graph = FALSE)
clara_time_Fish <- toc()
clara_clusters_Fish <- clara_result_Fish$cluster
table(true_lab_Fish, clara_clusters_Fish)
new <- 1:2
old <- c(2,1)
clara_clusters_Fish[clara_clusters_Fish %in% old] <- new[match(clara_clusters_Fish, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
clara_clusters_Fish[clara_clusters_Fish %in% old] <- new[match(clara_clusters_Fish, old, nomatch = 0)]
table(true_lab_Fish, clara_clusters_Fish)

# PAM clustering
tic("PAM Runtime")
pam_result_Fish <- eclust(Fish_scaled_data, "pam", G, graph = FALSE)
pam_time_Fish <- toc()
pam_clusters_Fish <- pam_result_Fish$cluster
table(true_lab_Fish, pam_clusters_Fish)
new <- 1:2
old <- c(2,1)
pam_clusters_Fish[pam_clusters_Fish %in% old] <- new[match(pam_clusters_Fish, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
pam_clusters_Fish[pam_clusters_Fish %in% old] <- new[match(pam_clusters_Fish, old, nomatch = 0)]
table(true_lab_Fish, pam_clusters_Fish)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Fish <- hclust(dist(Fish_scaled_data), method = "ward.D2")
hclust_time_Fish <- toc()
hclust_clusters_Fish <- cutree(hclust_result_Fish, k = G)
table(true_lab_Fish, hclust_clusters_Fish)
new <- c(2,3)
old <- c(3,2)
hclust_clusters_Fish[hclust_clusters_Fish %in% old] <- new[match(hclust_clusters_Fish, old, nomatch = 0)]
table(true_lab_Fish, hclust_clusters_Fish)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Fish <- Mclust(Fish_scaled_data, G = G)
mclust_time_Fish <- toc()
summary(mclust_result_Fish)
mclust_clusters_Fish <- mclust_result_Fish$classification
table(true_lab_Fish, mclust_clusters_Fish)
new <- 1:2
old <- c(2,1)
mclust_clusters_Fish[mclust_clusters_Fish %in% old] <- new[match(mclust_clusters_Fish, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
mclust_clusters_Fish[mclust_clusters_Fish %in% old] <- new[match(mclust_clusters_Fish, old, nomatch = 0)]
new <- c(5,6)
old <- c(6,5)
mclust_clusters_Fish[mclust_clusters_Fish %in% old] <- new[match(mclust_clusters_Fish, old, nomatch = 0)]
new <- c(6,7)
old <- c(7,6)
mclust_clusters_Fish[mclust_clusters_Fish %in% old] <- new[match(mclust_clusters_Fish, old, nomatch = 0)]
table(true_lab_Fish, mclust_clusters_Fish)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Fish_scaled_data, alphaPriors = c(6,19))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Fish <- toc()
# Extract clusters 
dpMVN_clusters_Fish <- as.numeric(dp$clusterLabels)
new <- c(1,5)
old <- c(5,1)
dpMVN_clusters_Fish[dpMVN_clusters_Fish %in% old] <- new[match(dpMVN_clusters_Fish, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
dpMVN_clusters_Fish[dpMVN_clusters_Fish %in% old] <- new[match(dpMVN_clusters_Fish, old, nomatch = 0)]
print(dpMVN_clusters_Fish)
table(true_lab_Fish, dpMVN_clusters_Fish)

calculate_clustering_metricsFish_2 <- function(dataFish_2, true_clusters_Fish, estimated_clusters_Fish_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Fish_list)) {
    clusters <- estimated_clusters_Fish_list[[method_name]]
    
    # Calculate metrics
    table_result_Fish <- table(true_clusters_Fish, clusters)
    kappa_result_Fish <- kappa2(data.frame(rater1 = true_clusters_Fish, rater2 = clusters))
    ari_result_Fish <- adjustedRandIndex(true_clusters_Fish, clusters)
    nmi_result_Fish <- NMI(true_clusters_Fish, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Fish,
      Kappa = kappa_result_Fish$value,
      ARI = ari_result_Fish,
      NMI = nmi_result_Fish,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Fish <- list(
  KMeans = kmeans_clusters_Fish,
  CLARA = clara_clusters_Fish,
  PAM = pam_clusters_Fish,
  Hierarchical = hclust_clusters_Fish,
  Mclust = mclust_clusters_Fish,
  DPMVN = dpMVN_clusters_Fish,
  True = true_lab_Fish
)

times_list_Fish <- list(
  KMeans = kmeans_time_Fish,
  CLARA = clara_time_Fish,
  PAM = pam_time_Fish,
  Hierarchical = hclust_time_Fish,
  Mclust = mclust_time_Fish,
  DPMVN = DPMVN_time_Fish
)

# Call the function
clustering_metricsFish_2 <- calculate_clustering_metricsFish_2(dataFish_2 = Fish_scaled_data, true_clusters_Fish = true_lab_Fish, estimated_clusters_Fish_list = cluster_result_Fish, times_list = times_list_Fish)

# Print the results for each method
for (method_name in names(clustering_metricsFish_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsFish_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsFish_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsFish_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsFish_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsFish_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Fish_K <- data.frame(clustering_metricsFish_2$True$Kappa, dp_Gmvlg_metricsFish_1$dp_Gmvlg$Kappa, clustering_metricsFish_2$KMeans$Kappa, clustering_metricsFish_2$CLARA$Kappa, clustering_metricsFish_2$PAM$Kappa, clustering_metricsFish_2$Hierarchical$Kappa, clustering_metricsFish_2$Mclust$Kappa, clustering_metricsFish_2$DPMVN$Kappa)
row_Fish_ARI <- data.frame(clustering_metricsFish_2$True$ARI, dp_Gmvlg_metricsFish_1$dp_Gmvlg$ARI, clustering_metricsFish_2$KMeans$ARI, clustering_metricsFish_2$CLARA$ARI, clustering_metricsFish_2$PAM$ARI, clustering_metricsFish_2$Hierarchical$ARI, clustering_metricsFish_2$Mclust$ARI, clustering_metricsFish_2$DPMVN$ARI)
row_Fish_NMI <- data.frame(clustering_metricsFish_2$True$NMI, dp_Gmvlg_metricsFish_1$dp_Gmvlg$NMI, clustering_metricsFish_2$KMeans$NMI, clustering_metricsFish_2$CLARA$NMI, clustering_metricsFish_2$PAM$NMI, clustering_metricsFish_2$Hierarchical$NMI, clustering_metricsFish_2$Mclust$NMI, clustering_metricsFish_2$DPMVN$NMI)
row_Fish_Cpu <- data.frame(dp_Gmvlg_metricsFish_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsFish_2$KMeans$CPU_RUN_TIME, clustering_metricsFish_2$CLARA$CPU_RUN_TIME, clustering_metricsFish_2$PAM$CPU_RUN_TIME, clustering_metricsFish_2$Hierarchical$CPU_RUN_TIME, clustering_metricsFish_2$Mclust$CPU_RUN_TIME, clustering_metricsFish_2$DPMVN$CPU_RUN_TIME)

tableFish_1 <- clustering_metricsFish_2$True$Table
tableFish_2 <- dp_Gmvlg_metricsFish_1$dp_Gmvlg$Table
tableFish_3 <- clustering_metricsFish_2$KMeans$Table
tableFish_4 <- clustering_metricsFish_2$CLARA$Table
tableFish_5 <- clustering_metricsFish_2$PAM$Table
tableFish_6 <- clustering_metricsFish_2$Hierarchical$Table
tableFish_7 <- clustering_metricsFish_2$Mclust$Table
tableFish_8 <- clustering_metricsFish_2$DPMVN$Table

colnames(row_Fish_K) <- NULL
colnames(row_Fish_ARI) <- NULL
colnames(row_Fish_NMI) <- NULL
colnames(row_Fish_Cpu) <- NULL

row_Fish_K <- as.matrix(row_Fish_K)
row_Fish_ARI <- as.matrix(row_Fish_ARI)
row_Fish_NMI <- as.matrix(row_Fish_NMI)
row_Fish_Cpu <- as.matrix(row_Fish_Cpu)

kappa_table_Fish <- rbind(row_Fish_K)
ARI_table_Fish <- rbind(row_Fish_ARI)
NMI_table_Fish <- rbind(row_Fish_NMI)
cpu_runtime_table_Fish <- rbind(row_Fish_Cpu)
clustering_table_Fish <- rbind(tableFish_1,tableFish_2,tableFish_3,tableFish_4,tableFish_5,tableFish_6,tableFish_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Fish)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Fish)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Fish)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Fish)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Fish)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/15_FISH CATCH/Cluster_metricsFish.xlsx", overwrite = TRUE)
