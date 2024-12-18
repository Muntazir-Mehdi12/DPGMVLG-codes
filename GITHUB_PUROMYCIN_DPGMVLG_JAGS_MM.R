#################################################################
#                     Puromycin DATASET                        ##
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
library(titanic)
library(rattle)
library(kernlab)
library(pgmm)
library(factoextra)
library(gridExtra)
library(clustertend)
library(hopkins)
library(clValid)
library(clusterSim)
library(dirichletprocess)
library(fossil)
library(corrplot)
library(aricode)
library(moments)

# Load the Puromycin dataset
data("Puromycin", package = "datasets")
help(Puromycin)
Pu_data <- as.matrix(Puromycin[,c(-3)]) # Convert all columns to numeric
true_lab_Pu <- as.numeric(Puromycin$state)    # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Pu_data)

# 2. Check for missing values
sum(is.na(Pu_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Adjusted for the number of features
for (i in 1:ncol(Pu_data)) {
  hist(Pu_data[, i], main=colnames(Pu_data)[i], xlab=colnames(Pu_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Pu_data, main="Puromycin Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[true_lab_Pu])

# 5. Boxplots for each feature by 'Type'
# Add normalized data to the original data frame
Pu_normalized <- as.data.frame(Pu_data)
Pu_normalized$type <- as.factor(Puromycin$state)

par(mfrow=c(1, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(Pu_data)) {
  boxplot(Pu_data[, i] ~ Pu_normalized$type, main=colnames(Pu_data)[i], xlab="Type", ylab=colnames(Pu_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Pu_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Puromycin dataset
Pu_scaled_data <- scale(Pu_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Pu_scaled_data, 2, skewness)
kurtosis_values <- apply(Pu_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Pu_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Pu_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Pu_scaled_data[,i])$out
  print(paste("Feature:", colnames(Pu_scaled_data)[i]))
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
outliers_matrix <- apply(Pu_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Pu_scaled_data, cluster = true_lab_Pu), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Pu_scaled_data)
D <- ncol(Pu_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Pu_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Pu_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Pu_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(25,20)          
       lambda_g[j, g] ~ dgamma(11,3) 
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
run_time_Pu <- toc()  # Capture the runtime
codaSamples_Pu = as.mcmc.list( runJagsOut )
summaryChains_Pu <- summary(codaSamples_Pu)

diagMCMC( codaObject=codaSamples_Pu , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Pu , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Pu , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Pu , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Pu , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Pu , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Pu , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Pu , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Pu$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Pu$statistics),1], P, G)
z_mode_Pu <- apply(matrix(summaryChains_Pu$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Pu$statistics),1], P, G),1, which.max)
z_mode_Pu

plot(Pu_scaled_data, col= z_mode_Pu, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Pu), col=unique(z_mode_Pu), pch=16, title="Cluster")
table( true_lab_Pu , z_mode_Pu)

# To switch to the same labels in the true clustering
#new <- c(1,2)
#old <- c(2,1)
#z_mode_Pu[z_mode_Pu %in% old] <- new[match(z_mode_Pu,old,nomatch = 0)]

plot(Pu_scaled_data, col= z_mode_Pu, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Pu), col=unique(z_mode_Pu), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Pu , z_mode_Pu)
kappa2(data.frame(rater1 = true_lab_Pu, rater2 = z_mode_Pu))

calculate_dp_Gmvlg_clustering_metrics <- function(dataPu_1, true_clusters_Pu) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Pu <- table(true_clusters_Pu, z_mode_Pu)
  kappa_result_Pu <- kappa2(data.frame(rater1 = true_clusters_Pu, rater2 = z_mode_Pu))
  ari_result_Pu <- adjustedRandIndex(true_clusters_Pu, z_mode_Pu)
  nmi_result_Pu <- NMI(true_clusters_Pu, z_mode_Pu)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Pu,
    Kappa = kappa_result_Pu$value,
    ARI = ari_result_Pu,
    NMI = nmi_result_Pu,
    CPU_RUN_TIME = run_time_Pu$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsPu_1 <- calculate_dp_Gmvlg_clustering_metrics(dataPu_1 = Pu_scaled_data, true_clusters_Pu = true_lab_Pu)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsPu_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Pu <- eclust(Pu_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Pu <- toc()
kmeans_clusters_Pu <- kmeans_result_Pu$cluster
table(true_lab_Pu, kmeans_clusters_Pu)
new <- 1:2
old <- c(2,1)
kmeans_clusters_Pu[kmeans_clusters_Pu %in% old] <- new[match(kmeans_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, kmeans_clusters_Pu)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Pu <- eclust(Pu_scaled_data, "clara", G, graph = FALSE)
clara_time_Pu <- toc()
clara_clusters_Pu <- clara_result_Pu$cluster
table(true_lab_Pu, clara_clusters_Pu)
new <- 1:2
old <- c(2,1)
clara_clusters_Pu[clara_clusters_Pu %in% old] <- new[match(clara_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, clara_clusters_Pu)

# PAM clustering
tic("PAM Runtime")
pam_result_Pu <- eclust(Pu_scaled_data, "pam", G, graph = FALSE)
pam_time_Pu <- toc()
pam_clusters_Pu <- pam_result_Pu$cluster
table(true_lab_Pu, pam_clusters_Pu)
new <- 1:2
old <- c(2,1)
pam_clusters_Pu[pam_clusters_Pu %in% old] <- new[match(pam_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, pam_clusters_Pu)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Pu <- hclust(dist(Pu_scaled_data), method = "ward.D2")
hclust_time_Pu <- toc()
hclust_clusters_Pu <- cutree(hclust_result_Pu, k = G)
table(true_lab_Pu, hclust_clusters_Pu)
new <- 1:2
old <- c(2,1)
hclust_clusters_Pu[hclust_clusters_Pu %in% old] <- new[match(hclust_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, hclust_clusters_Pu)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Pu <- Mclust(Pu_scaled_data, G = G)
mclust_time_Pu <- toc()
summary(mclust_result_Pu)
mclust_clusters_Pu <- mclust_result_Pu$classification
table(true_lab_Pu, mclust_clusters_Pu)
new <- 1:2
old <- c(2,1)
mclust_clusters_Pu[mclust_clusters_Pu %in% old] <- new[match(mclust_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, mclust_clusters_Pu)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Pu_scaled_data, alphaPriors = c(1,1))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Pu <- toc()
# Extract clusters 
dpMVN_clusters_Pu <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Pu)
table(true_lab_Pu, dpMVN_clusters_Pu)
new <- 1:2
old <- c(2,1)
dpMVN_clusters_Pu[dpMVN_clusters_Pu %in% old] <- new[match(dpMVN_clusters_Pu, old, nomatch = 0)]
table(true_lab_Pu, dpMVN_clusters_Pu)
kappa2(data.frame(rater1 = true_lab_Pu, rater2 = dpMVN_clusters_Pu))


calculate_clustering_metricsPu_2 <- function(dataPu_2, true_clusters_Pu, estimated_clusters_Pu_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Pu_list)) {
    clusters <- estimated_clusters_Pu_list[[method_name]]
    
    # Calculate metrics
    table_result_Pu <- table(true_clusters_Pu, clusters)
    kappa_result_Pu <- kappa2(data.frame(rater1 = true_clusters_Pu, rater2 = clusters))
    ari_result_Pu <- adjustedRandIndex(true_clusters_Pu, clusters)
    nmi_result_Pu <- NMI(true_clusters_Pu, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Pu,
      Kappa = kappa_result_Pu$value,
      ARI = ari_result_Pu,
      NMI = nmi_result_Pu,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Pu <- list(
  KMeans = kmeans_clusters_Pu,
  CLARA = clara_clusters_Pu,
  PAM = pam_clusters_Pu,
  Hierarchical = hclust_clusters_Pu,
  Mclust = mclust_clusters_Pu,
  DPMVN = dpMVN_clusters_Pu,
  True = true_lab_Pu
)

times_list_Pu <- list(
  KMeans = kmeans_time_Pu,
  CLARA = clara_time_Pu,
  PAM = pam_time_Pu,
  Hierarchical = hclust_time_Pu,
  Mclust = mclust_time_Pu,
  DPMVN = DPMVN_time_Pu
)

# Call the function
clustering_metricsPu_2 <- calculate_clustering_metricsPu_2(dataPu_2 = Pu_scaled_data, true_clusters_Pu = true_lab_Pu, estimated_clusters_Pu_list = cluster_result_Pu, times_list = times_list_Pu)

# Print the results for each method
for (method_name in names(clustering_metricsPu_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsPu_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsPu_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsPu_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsPu_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsPu_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Pu_K <- data.frame(clustering_metricsPu_2$True$Kappa, dp_Gmvlg_metricsPu_1$dp_Gmvlg$Kappa, clustering_metricsPu_2$KMeans$Kappa, clustering_metricsPu_2$CLARA$Kappa, clustering_metricsPu_2$PAM$Kappa, clustering_metricsPu_2$Hierarchical$Kappa, clustering_metricsPu_2$Mclust$Kappa, clustering_metricsPu_2$DPMVN$Kappa)
row_Pu_ARI <- data.frame(clustering_metricsPu_2$True$ARI, dp_Gmvlg_metricsPu_1$dp_Gmvlg$ARI, clustering_metricsPu_2$KMeans$ARI, clustering_metricsPu_2$CLARA$ARI, clustering_metricsPu_2$PAM$ARI, clustering_metricsPu_2$Hierarchical$ARI, clustering_metricsPu_2$Mclust$ARI, clustering_metricsPu_2$DPMVN$ARI)
row_Pu_NMI <- data.frame(clustering_metricsPu_2$True$NMI, dp_Gmvlg_metricsPu_1$dp_Gmvlg$NMI, clustering_metricsPu_2$KMeans$NMI, clustering_metricsPu_2$CLARA$NMI, clustering_metricsPu_2$PAM$NMI, clustering_metricsPu_2$Hierarchical$NMI, clustering_metricsPu_2$Mclust$NMI, clustering_metricsPu_2$DPMVN$NMI)
row_Pu_Cpu <- data.frame(dp_Gmvlg_metricsPu_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsPu_2$KMeans$CPU_RUN_TIME, clustering_metricsPu_2$CLARA$CPU_RUN_TIME, clustering_metricsPu_2$PAM$CPU_RUN_TIME, clustering_metricsPu_2$Hierarchical$CPU_RUN_TIME, clustering_metricsPu_2$Mclust$CPU_RUN_TIME, clustering_metricsPu_2$DPMVN$CPU_RUN_TIME)

tablePu_1 <- clustering_metricsPu_2$True$Table
tablePu_2 <- dp_Gmvlg_metricsPu_1$dp_Gmvlg$Table
tablePu_3 <- clustering_metricsPu_2$KMeans$Table
tablePu_4 <- clustering_metricsPu_2$CLARA$Table
tablePu_5 <- clustering_metricsPu_2$PAM$Table
tablePu_6 <- clustering_metricsPu_2$Hierarchical$Table
tablePu_7 <- clustering_metricsPu_2$Mclust$Table
tablePu_8 <- clustering_metricsPu_2$DPMVN$Table

colnames(row_Pu_K) <- NULL
colnames(row_Pu_ARI) <- NULL
colnames(row_Pu_NMI) <- NULL
colnames(row_Pu_Cpu) <- NULL

row_Pu_K <- as.matrix(row_Pu_K)
row_Pu_ARI <- as.matrix(row_Pu_ARI)
row_Pu_NMI <- as.matrix(row_Pu_NMI)
row_Pu_Cpu <- as.matrix(row_Pu_Cpu)

kappa_table_Pu <- rbind(row_Pu_K)
ARI_table_Pu <- rbind(row_Pu_ARI)
NMI_table_Pu <- rbind(row_Pu_NMI)
cpu_runtime_table_Pu <- rbind(row_Pu_Cpu)
clustering_table_Pu <- rbind(tablePu_1,tablePu_2,tablePu_3,tablePu_4,tablePu_5,tablePu_6,tablePu_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Pu)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Pu)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Pu)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Pu)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Pu)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/03_Puromycin/Cluster_metricsPu.xlsx", overwrite = TRUE)
