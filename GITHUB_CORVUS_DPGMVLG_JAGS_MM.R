#################################################################
#                       CORVUS DATASET                         ##
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
library(ade4)
library(aricode)
library(moments)

# Load the Corvus dataset
data("corvus",package ="ade4")
help(corvus)
Corvus_data <- as.matrix(corvus[, -c(3,4)])
true_lab_Corvus=as.numeric(corvus[,3])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Corvus_data)

# 2. Check for missing values
sum(is.na(Corvus_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(Corvus_data)) {
  hist(Corvus_data[, i], main=colnames(Corvus_data)[i], xlab=colnames(Corvus_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Corvus_data, main="Corvus Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_Corvus])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(Corvus_data)
Cor_normalized$habitat <- corvus$habitat

par(mfrow=c(1, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(Corvus_data)) {
  boxplot(Corvus_data[, i] ~ Cor_normalized$habitat, main=colnames(Corvus_data)[i], xlab="Habitat", ylab=colnames(Corvus_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Corvus_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the corvus dataset
Corvus_scaled_data <- scale(Corvus_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Corvus_scaled_data, 2, skewness)
kurtosis_values <- apply(Corvus_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Corvus_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Corvus_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Corvus_scaled_data[,i])$out
  print(paste("Feature:", colnames(Corvus_scaled_data)[i]))
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
outliers_matrix <- apply(Corvus_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Corvus_scaled_data, cluster = true_lab_Corvus), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Corvus_scaled_data)
D <- ncol(Corvus_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Corvus_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Corvus_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Corvus_scaled_data))^(1 / (D - 1))

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
    alpha[g] ~ dgamma(11,5)
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
run_time_Corvus <- toc()
codaSamples_Corvus = as.mcmc.list( runJagsOut )
summaryChains_Corvus <- summary(codaSamples_Corvus)

diagMCMC( codaObject=codaSamples_Corvus , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Corvus , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Corvus$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Corvus$statistics),1], P, G)
z_mode_Corvus <- apply(matrix(summaryChains_Corvus$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Corvus$statistics),1], P, G),1, which.max)
z_mode_Corvus

plot(Corvus_scaled_data, col= z_mode_Corvus, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Corvus), col=unique(z_mode_Corvus), pch=16, title="Cluster")
table( true_lab_Corvus , z_mode_Corvus)

# To switch to the same labels in the true clustering
new <- 1:2
old <- c(2,1)
z_mode_Corvus[z_mode_Corvus %in% old] <- new[match(z_mode_Corvus,old,nomatch = 0)]

plot(Corvus_scaled_data, col= z_mode_Corvus, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Corvus), col=unique(z_mode_Corvus), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table(true_lab_Corvus , z_mode_Corvus)
kappa2(data.frame(rater1 = true_lab_Corvus, rater2 = z_mode_Corvus))

calculate_dp_Gmvlg_clustering_metrics <- function(dataCorvus_1, true_clusters_Corvus) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Corvus <- table(true_clusters_Corvus, z_mode_Corvus)
  kappa_result_Corvus <- kappa2(data.frame(rater1 = true_clusters_Corvus, rater2 = z_mode_Corvus))
  ari_result_Corvus <- adjustedRandIndex(true_clusters_Corvus, z_mode_Corvus)
  nmi_result_Corvus <- NMI(true_clusters_Corvus, z_mode_Corvus)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Corvus,
    Kappa = kappa_result_Corvus$value,
    ARI = ari_result_Corvus,
    NMI = nmi_result_Corvus,
    CPU_RUN_TIME = run_time_Corvus$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsCorvus_1 <- calculate_dp_Gmvlg_clustering_metrics(dataCorvus_1 = Corvus_scaled_data, true_clusters_Corvus = true_lab_Corvus)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsCorvus_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsCorvus_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsCorvus_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsCorvus_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsCorvus_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsCorvus_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Corvus <- eclust(Corvus_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Corvus <- toc()
kmeans_clusters_Corvus <- kmeans_result_Corvus$cluster
table(true_lab_Corvus, kmeans_clusters_Corvus)
new <- 1:2
old <- c(2,1)
kmeans_clusters_Corvus[kmeans_clusters_Corvus %in% old] <- new[match(kmeans_clusters_Corvus, old, nomatch = 0)]
table(true_lab_Corvus, kmeans_clusters_Corvus)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Corvus <- eclust(Corvus_scaled_data, "clara", G, graph = FALSE)
clara_time_Corvus <- toc()
clara_clusters_Corvus <- clara_result_Corvus$cluster
table(true_lab_Corvus, clara_clusters_Corvus)
new <- 1:2
old <- c(2,1)
clara_clusters_Corvus[clara_clusters_Corvus %in% old] <- new[match(clara_clusters_Corvus, old, nomatch = 0)]
table(true_lab_Corvus, clara_clusters_Corvus)

# PAM clustering
tic("PAM Runtime")
pam_result_Corvus <- eclust(Corvus_scaled_data, "pam", G, graph = FALSE)
pam_time_Corvus <- toc()
pam_clusters_Corvus <- pam_result_Corvus$cluster
table(true_lab_Corvus, pam_clusters_Corvus)
new <- 1:2
old <- c(2,1)
pam_clusters_Corvus[pam_clusters_Corvus %in% old] <- new[match(pam_clusters_Corvus, old, nomatch = 0)]
table(true_lab_Corvus, pam_clusters_Corvus)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Corvus <- hclust(dist(Corvus_scaled_data), method = "ward.D2")
hclust_time_Corvus <- toc()
hclust_clusters_Corvus <- cutree(hclust_result_Corvus, k = G)
table(true_lab_Corvus, hclust_clusters_Corvus)
new <- 1:2
old <- c(2,1)
hclust_clusters_Corvus[hclust_clusters_Corvus %in% old] <- new[match(hclust_clusters_Corvus, old, nomatch = 0)]
table(true_lab_Corvus, hclust_clusters_Corvus)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Corvus <- Mclust(Corvus_scaled_data, G = G)
mclust_time_Corvus <- toc()
summary(mclust_result_Corvus)
mclust_clusters_Corvus <- mclust_result_Corvus$classification
table(true_lab_Corvus, mclust_clusters_Corvus)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Corvus_scaled_data, alphaPriors = c(7,19))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Corvus <- toc()
# Extract clusters 
dpMVN_clusters_Corvus <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Corvus)
table(true_lab_Corvus, dpMVN_clusters_Corvus)
new <- 1:2
old <- c(2,1)
dpMVN_clusters_Corvus[dpMVN_clusters_Corvus %in% old] <- new[match(dpMVN_clusters_Corvus, old, nomatch = 0)]
print(dpMVN_clusters_Corvus)
table(true_lab_Corvus, dpMVN_clusters_Corvus)


calculate_clustering_metricsCorvus_2 <- function(dataCorvus_2, true_clusters_Corvus, estimated_clusters_Corvus_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Corvus_list)) {
    clusters <- estimated_clusters_Corvus_list[[method_name]]
    
    # Calculate metrics
    table_result_Corvus <- table(true_clusters_Corvus, clusters)
    kappa_result_Corvus <- kappa2(data.frame(rater1 = true_clusters_Corvus, rater2 = clusters))
    ari_result_Corvus <- adjustedRandIndex(true_clusters_Corvus, clusters)
    nmi_result_Corvus <- NMI(true_clusters_Corvus, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Corvus,
      Kappa = kappa_result_Corvus$value,
      ARI = ari_result_Corvus,
      NMI = nmi_result_Corvus,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Corvus <- list(
  KMeans = kmeans_clusters_Corvus,
  CLARA = clara_clusters_Corvus,
  PAM = pam_clusters_Corvus,
  Hierarchical = hclust_clusters_Corvus,
  Mclust = mclust_clusters_Corvus,
  DPMVN = dpMVN_clusters_Corvus,
  True = true_lab_Corvus
)

times_list_Pu <- list(
  KMeans = kmeans_time_Corvus,
  CLARA = clara_time_Corvus,
  PAM = pam_time_Corvus,
  Hierarchical = hclust_time_Corvus,
  Mclust = mclust_time_Corvus,
  DPMVN = DPMVN_time_Corvus
)

# Call the function
clustering_metricsCorvus_2 <- calculate_clustering_metricsCorvus_2(dataCorvus_2 = Corvus_scaled_data, true_clusters_Corvus = true_lab_Corvus, estimated_clusters_Corvus_list = cluster_result_Corvus, times_list = times_list_Pu)

# Print the results for each method
for (method_name in names(clustering_metricsCorvus_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsCorvus_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsCorvus_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsCorvus_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsCorvus_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsCorvus_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Corvus_K <- data.frame(clustering_metricsCorvus_2$True$Kappa, dp_Gmvlg_metricsCorvus_1$dp_Gmvlg$Kappa, clustering_metricsCorvus_2$KMeans$Kappa, clustering_metricsCorvus_2$CLARA$Kappa, clustering_metricsCorvus_2$PAM$Kappa, clustering_metricsCorvus_2$Hierarchical$Kappa, clustering_metricsCorvus_2$Mclust$Kappa, clustering_metricsCorvus_2$DPMVN$Kappa)
row_Corvus_ARI <- data.frame(clustering_metricsCorvus_2$True$ARI, dp_Gmvlg_metricsCorvus_1$dp_Gmvlg$ARI, clustering_metricsCorvus_2$KMeans$ARI, clustering_metricsCorvus_2$CLARA$ARI, clustering_metricsCorvus_2$PAM$ARI, clustering_metricsCorvus_2$Hierarchical$ARI, clustering_metricsCorvus_2$Mclust$ARI, clustering_metricsCorvus_2$DPMVN$ARI)
row_Corvus_NMI <- data.frame(clustering_metricsCorvus_2$True$NMI, dp_Gmvlg_metricsCorvus_1$dp_Gmvlg$NMI, clustering_metricsCorvus_2$KMeans$NMI, clustering_metricsCorvus_2$CLARA$NMI, clustering_metricsCorvus_2$PAM$NMI, clustering_metricsCorvus_2$Hierarchical$NMI, clustering_metricsCorvus_2$Mclust$NMI, clustering_metricsCorvus_2$DPMVN$NMI)
row_Corvus_Cpu <- data.frame(dp_Gmvlg_metricsCorvus_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsCorvus_2$KMeans$CPU_RUN_TIME, clustering_metricsCorvus_2$CLARA$CPU_RUN_TIME, clustering_metricsCorvus_2$PAM$CPU_RUN_TIME, clustering_metricsCorvus_2$Hierarchical$CPU_RUN_TIME, clustering_metricsCorvus_2$Mclust$CPU_RUN_TIME, clustering_metricsCorvus_2$DPMVN$CPU_RUN_TIME)

tableCorvus_1 <- clustering_metricsCorvus_2$True$Table
tableCorvus_2 <- dp_Gmvlg_metricsCorvus_1$dp_Gmvlg$Table
tableCorvus_3 <- clustering_metricsCorvus_2$KMeans$Table
tableCorvus_4 <- clustering_metricsCorvus_2$CLARA$Table
tableCorvus_5 <- clustering_metricsCorvus_2$PAM$Table
tableCorvus_6 <- clustering_metricsCorvus_2$Hierarchical$Table
tableCorvus_7 <- clustering_metricsCorvus_2$Mclust$Table
tableCorvus_8 <- clustering_metricsCorvus_2$DPMVN$Table

colnames(row_Corvus_K) <- NULL
colnames(row_Corvus_ARI) <- NULL
colnames(row_Corvus_NMI) <- NULL
colnames(row_Corvus_Cpu) <- NULL

row_Corvus_K <- as.matrix(row_Corvus_K)
row_Corvus_ARI <- as.matrix(row_Corvus_ARI)
row_Corvus_NMI <- as.matrix(row_Corvus_NMI)
row_Corvus_Cpu <- as.matrix(row_Corvus_Cpu)

kappa_table_Corvus <- rbind(row_Corvus_K)
ARI_table_Corvus <- rbind(row_Corvus_ARI)
NMI_table_Corvus <- rbind(row_Corvus_NMI)
cpu_runtime_table_Corvus <- rbind(row_Corvus_Cpu)
clustering_table_Corvus <- rbind(tableCorvus_1,tableCorvus_2,tableCorvus_3,tableCorvus_4,tableCorvus_5,tableCorvus_6,tableCorvus_7, tableCorvus_8)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Corvus)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Corvus)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Corvus)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Corvus)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Corvus)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/06_CORVUS/Cluster_metricsCorvus.xlsx", overwrite = TRUE)

