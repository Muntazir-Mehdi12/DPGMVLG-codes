#################################################################
#                       CarniHerbi49 DATASET                   ##
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
data("carniherbi49",package ="ade4")
carniherbi49                            
help(carniherbi49)
str(carniherbi49)
uniques <- lapply(carniherbi49, unique)
uniques
Carni_data <- as.matrix(carniherbi49$tab2[, -c(1)])
true_lab_Carni=as.numeric(carniherbi49$tab2[,1])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Carni_data)

# 2. Check for missing values
sum(is.na(Carni_data))

# 3. Distribution of each feature
par(mfrow=c(2, 2))  # Adjusted for the number of features
for (i in 1:ncol(Carni_data)) {
  hist(Carni_data[, i], main=colnames(Carni_data)[i], xlab=colnames(Carni_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Carni_data, main="Carniherbi49 Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[true_lab_Carni])

# 5. Boxplots for each feature by 'Weight Category'
# Add normalized data to the original data frame
mt_normalized <- as.data.frame(Carni_data)
mt_normalized$am <- carniherbi49$tab2$clade

par(mfrow=c(2, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(Carni_data)) {
  boxplot(Carni_data[, i] ~ mt_normalized$am, main=colnames(Carni_data)[i], xlab="dietary Category", ylab=colnames(Carni_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Carni_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Corvus dataset
Carni_scaled_data <- scale(Carni_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Carni_scaled_data, 2, skewness)
kurtosis_values <- apply(Carni_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Carni_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Carni_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Carni_scaled_data[,i])$out
  print(paste("Feature:", colnames(Carni_scaled_data)[i]))
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
outliers_matrix <- apply(Carni_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Carni_scaled_data, cluster = true_lab_Carni), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Carni_scaled_data)
D <- ncol(Carni_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Carni_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Carni_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Carni_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(4,4)        
       lambda_g[j, g] ~ dgamma(1,5) 
    }
    alpha[g] <- 5
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
run_time_carni <- toc()
codaSamples_carni = as.mcmc.list( runJagsOut )
summaryChains_carni <- summary(codaSamples_carni)

diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_carni$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_carni$statistics),1], P, G)
z_mode_Carni <- apply(matrix(summaryChains_carni$statistics[(1+4*G+4*G+P+1):nrow(summaryChains_carni$statistics),1], P, G),1, which.max)
z_mode_Carni

plot(Carni_scaled_data, col= z_mode_Carni, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Carni), col=unique(z_mode_Carni), pch=16, title="Cluster")
table( true_lab_Carni , z_mode_Carni)

# To switch to the same labels in the true clustering
new <- c(1,2)
old <- c(2,1)
z_mode_Carni[z_mode_Carni %in% old] <- new[match(z_mode_Carni,old,nomatch = 0)]

plot(Carni_scaled_data, col= z_mode_Carni, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Carni), col=unique(z_mode_Carni), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Carni , z_mode_Carni)
kappa2(data.frame(rater1 = true_lab_Carni, rater2 = z_mode_Carni))

calculate_dp_Gmvlg_clustering_metrics <- function(dataCarni_1, true_clusters_Carni) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Carni <- table(true_clusters_Carni, z_mode_Carni)
  kappa_result_Carni <- kappa2(data.frame(rater1 = true_clusters_Carni, rater2 = z_mode_Carni))
  ari_result_Carni <- adjustedRandIndex(true_clusters_Carni, z_mode_Carni)
  nmi_result_Carni <- NMI(true_clusters_Carni, z_mode_Carni)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Carni,
    Kappa = kappa_result_Carni$value,
    ARI = ari_result_Carni,
    NMI = nmi_result_Carni,
    CPU_RUN_TIME = run_time_carni$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsCarni_1 <- calculate_dp_Gmvlg_clustering_metrics(dataCarni_1 = Carni_scaled_data, true_clusters_Carni = true_lab_Carni)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsCarni_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsCarni_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsCarni_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsCarni_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsCarni_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsCarni_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Carni <- eclust(Carni_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Carni <- toc()
kmeans_clusters_Carni <- kmeans_result_Carni$cluster
table(true_lab_Carni, kmeans_clusters_Carni)
new <- 1:2
old <- c(2,1)
kmeans_clusters_Carni[kmeans_clusters_Carni %in% old] <- new[match(kmeans_clusters_Carni, old, nomatch = 0)]
table(true_lab_Carni, kmeans_clusters_Carni)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Carni <- eclust(Carni_scaled_data, "clara", G, graph = FALSE)
clara_time_Carni <- toc()
clara_clusters_Carni <- clara_result_Carni$cluster
table(true_lab_Carni, clara_clusters_Carni)

# PAM clustering
tic("PAM Runtime")
pam_result_Carni <- eclust(Carni_scaled_data, "pam", G, graph = FALSE)
pam_time_Carni <- toc()
pam_clusters_Carni <- pam_result_Carni$cluster
table(true_lab_Carni, pam_clusters_Carni)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Carni <- hclust(dist(Carni_scaled_data), method = "ward.D2")
hclust_time_Carni <- toc()
hclust_clusters_Carni <- cutree(hclust_result_Carni, k = G)
table(true_lab_Carni, hclust_clusters_Carni)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Carni <- Mclust(Carni_scaled_data, G = G)
mclust_time_Carni <- toc()
summary(mclust_result_Carni)
mclust_clusters_Carni <- mclust_result_Carni$classification
table(true_lab_Carni, mclust_clusters_Carni)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Carni_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Carni <- toc()
# Extract clusters 
dpMVN_clusters_Carni <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Carni)
table(true_lab_Carni, dpMVN_clusters_Carni)

calculate_clustering_metricsCarni_2 <- function(dataCarni_2, true_clusters_Carni, estimated_clusters_Carni_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Carni_list)) {
    clusters <- estimated_clusters_Carni_list[[method_name]]
    
    # Calculate metrics
    table_result_Carni <- table(true_clusters_Carni, clusters)
    kappa_result_Carni <- kappa2(data.frame(rater1 = true_clusters_Carni, rater2 = clusters))
    ari_result_Carni <- adjustedRandIndex(true_clusters_Carni, clusters)
    nmi_result_Carni <- NMI(true_clusters_Carni, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Carni,
      Kappa = kappa_result_Carni$value,
      ARI = ari_result_Carni,
      NMI = nmi_result_Carni,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Carni <- list(
  KMeans = kmeans_clusters_Carni,
  CLARA = clara_clusters_Carni,
  PAM = pam_clusters_Carni,
  Hierarchical = hclust_clusters_Carni,
  Mclust = mclust_clusters_Carni,
  DPMVN = dpMVN_clusters_Carni,
  True = true_lab_Carni
)

times_list_Pu <- list(
  KMeans = kmeans_time_Carni,
  CLARA = clara_time_Carni,
  PAM = pam_time_Carni,
  Hierarchical = hclust_time_Carni,
  Mclust = mclust_time_Carni,
  DPMVN = DPMVN_time_Carni
)

# Call the function
clustering_metricsCarni_2 <- calculate_clustering_metricsCarni_2(dataCarni_2 = Carni_scaled_data, true_clusters_Carni = true_lab_Carni, estimated_clusters_Carni_list = cluster_result_Carni, times_list = times_list_Pu)

# Print the results for each method
for (method_name in names(clustering_metricsCarni_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsCarni_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsCarni_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsCarni_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsCarni_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsCarni_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Carni_K <- data.frame(clustering_metricsCarni_2$True$Kappa, dp_Gmvlg_metricsCarni_1$dp_Gmvlg$Kappa, clustering_metricsCarni_2$KMeans$Kappa, clustering_metricsCarni_2$CLARA$Kappa, clustering_metricsCarni_2$PAM$Kappa, clustering_metricsCarni_2$Hierarchical$Kappa, clustering_metricsCarni_2$Mclust$Kappa, clustering_metricsCarni_2$DPMVN$Kappa)
row_Carni_ARI <- data.frame(clustering_metricsCarni_2$True$ARI, dp_Gmvlg_metricsCarni_1$dp_Gmvlg$ARI, clustering_metricsCarni_2$KMeans$ARI, clustering_metricsCarni_2$CLARA$ARI, clustering_metricsCarni_2$PAM$ARI, clustering_metricsCarni_2$Hierarchical$ARI, clustering_metricsCarni_2$Mclust$ARI, clustering_metricsCarni_2$DPMVN$ARI)
row_Carni_NMI <- data.frame(clustering_metricsCarni_2$True$NMI, dp_Gmvlg_metricsCarni_1$dp_Gmvlg$NMI, clustering_metricsCarni_2$KMeans$NMI, clustering_metricsCarni_2$CLARA$NMI, clustering_metricsCarni_2$PAM$NMI, clustering_metricsCarni_2$Hierarchical$NMI, clustering_metricsCarni_2$Mclust$NMI, clustering_metricsCarni_2$DPMVN$NMI)
row_Carni_Cpu <- data.frame(dp_Gmvlg_metricsCarni_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsCarni_2$KMeans$CPU_RUN_TIME, clustering_metricsCarni_2$CLARA$CPU_RUN_TIME, clustering_metricsCarni_2$PAM$CPU_RUN_TIME, clustering_metricsCarni_2$Hierarchical$CPU_RUN_TIME, clustering_metricsCarni_2$Mclust$CPU_RUN_TIME, clustering_metricsCarni_2$DPMVN$CPU_RUN_TIME)

tableCarni_1 <- clustering_metricsCarni_2$True$Table
tableCarni_2 <- dp_Gmvlg_metricsCarni_1$dp_Gmvlg$Table
tableCarni_3 <- clustering_metricsCarni_2$KMeans$Table
tableCarni_4 <- clustering_metricsCarni_2$CLARA$Table
tableCarni_5 <- clustering_metricsCarni_2$PAM$Table
tableCarni_6 <- clustering_metricsCarni_2$Hierarchical$Table
tableCarni_7 <- clustering_metricsCarni_2$Mclust$Table
tableCarni_8 <- clustering_metricsCarni_2$DPMVN$Table

colnames(row_Carni_K) <- NULL
colnames(row_Carni_ARI) <- NULL
colnames(row_Carni_NMI) <- NULL
colnames(row_Carni_Cpu) <- NULL

row_Carni_K <- as.matrix(row_Carni_K)
row_Carni_ARI <- as.matrix(row_Carni_ARI)
row_Carni_NMI <- as.matrix(row_Carni_NMI)
row_Carni_Cpu <- as.matrix(row_Carni_Cpu)

kappa_table_Carni <- rbind(row_Carni_K)
ARI_table_Carni <- rbind(row_Carni_ARI)
NMI_table_Carni <- rbind(row_Carni_NMI)
cpu_runtime_table_Carni <- rbind(row_Carni_Cpu)
clustering_table_Carni <- rbind(tableCarni_1,tableCarni_2,tableCarni_3,tableCarni_4,tableCarni_5,tableCarni_6,tableCarni_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Carni)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Carni)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Carni)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Carni)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Carni)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/09_CARNIHERBI/Cluster_metricsCarni.xlsx", overwrite = TRUE)

