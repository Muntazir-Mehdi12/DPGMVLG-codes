#################################################################
#                         DUNE DATASET                         ##
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

# Load the Dune dataset
data("dunedata", package ="ade4")
dunedata$envir 
Dune_data <- as.matrix(dunedata$envir[, -c(4,5)])
true_lab_Dune=as.numeric(dunedata$envir[,5])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Dune_data)

# 2. Check for missing values
sum(is.na(Dune_data))

# 3. Distribution of each feature
par(mfrow=c(2, 2))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(Dune_data)) {
  hist(Dune_data[, i], main=colnames(Dune_data)[i], xlab=colnames(Dune_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Dune_data, main="Dune Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_Dune])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(Dune_data)
Cor_normalized$management <- dunedata$envir$management

par(mfrow=c(2, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(Dune_data)) {
  boxplot(Dune_data[, i] ~ Cor_normalized$management, main=colnames(Dune_data)[i], xlab="management", ylab=colnames(Dune_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Dune_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Dune dataset
Dune_scaled_data <- scale(Dune_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Dune_scaled_data, 2, skewness)
kurtosis_values <- apply(Dune_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Dune_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Dune_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Dune_scaled_data[,i])$out
  print(paste("Feature:", colnames(Dune_scaled_data)[i]))
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
outliers_matrix <- apply(Dune_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Dune_scaled_data, cluster = true_lab_Dune), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Dune_scaled_data)
D <- ncol(Dune_scaled_data)

#Try with different number of clusters
G <- 4  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Dune_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Dune_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Dune_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(4,2)   # dgamma(25,20)       
       lambda_g[j, g] ~ dgamma(1,6) # dgamma(8,3)
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
run_time_Dune<- toc()
codaSamples_Dune = as.mcmc.list(runJagsOut)
summaryChains_Dune <- summary(codaSamples_Dune)

diagMCMC( codaObject=codaSamples_Dune , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Dune , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Dune , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Dune , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Dune , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Dune , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Dune , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Dune , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Dune$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_Dune$statistics),1], P, G)
z_mode_Dune <- apply(matrix(summaryChains_Dune$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_Dune$statistics),1], P, G),1, which.max)
z_mode_Dune

plot(Dune_scaled_data, col= z_mode_Dune, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Dune), col=unique(z_mode_Dune), pch=16, title="Cluster")
table(true_lab_Dune , z_mode_Dune)

# To switch to the same labels in the true clustering
new <- c(1,3)
old <- c(3,1)
z_mode_Dune[z_mode_Dune %in% old] <- new[match(z_mode_Dune,old,nomatch = 0)]

plot(Dune_scaled_data, col= z_mode_Dune, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Dune), col=unique(z_mode_Dune), pch=16, title="Cluster")
table(true_lab_Dune , z_mode_Dune)

calculate_dp_Gmvlg_clustering_metrics <- function(dataDune_1, true_clusters_Dune) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Dune <- table(true_clusters_Dune, z_mode_Dune)
  kappa_result_Dune <- kappa2(data.frame(rater1 = true_clusters_Dune, rater2 = z_mode_Dune))
  ari_result_Dune <- adjustedRandIndex(true_clusters_Dune, z_mode_Dune)
  nmi_result_Dune <- NMI(true_clusters_Dune, z_mode_Dune)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Dune,
    Kappa = kappa_result_Dune$value,
    ARI = ari_result_Dune,
    NMI = nmi_result_Dune,
    CPU_RUN_TIME = run_time_Dune$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsDune_1 <- calculate_dp_Gmvlg_clustering_metrics(dataDune_1 = Dune_scaled_data, true_clusters_Dune = true_lab_Dune)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsDune_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsDune_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsDune_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsDune_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsDune_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsDune_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Dune <- eclust(Dune_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Dune <- toc()
kmeans_clusters_Dune <- kmeans_result_Dune$cluster
table(true_lab_Dune, kmeans_clusters_Dune)
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_Dune[kmeans_clusters_Dune %in% old] <- new[match(kmeans_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, kmeans_clusters_Dune)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Dune <- eclust(Dune_scaled_data, "clara", G, graph = FALSE)
clara_time_Dune <- toc()
clara_clusters_Dune <- clara_result_Dune$cluster
table(true_lab_Dune, clara_clusters_Dune)
new <- c(3,4)
old <- c(4,3)
clara_clusters_Dune[clara_clusters_Dune %in% old] <- new[match(clara_clusters_Dune, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
clara_clusters_Dune[clara_clusters_Dune %in% old] <- new[match(clara_clusters_Dune, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
clara_clusters_Dune[clara_clusters_Dune %in% old] <- new[match(clara_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, clara_clusters_Dune)

# PAM clustering
tic("PAM Runtime")
pam_result_Dune <- eclust(Dune_scaled_data, "pam", G, graph = FALSE)
pam_time_Dune <- toc()
pam_clusters_Dune <- pam_result_Dune$cluster
table(true_lab_Dune, pam_clusters_Dune)
new <- c(3,4)
old <- c(4,3)
pam_clusters_Dune[pam_clusters_Dune %in% old] <- new[match(pam_clusters_Dune, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
pam_clusters_Dune[pam_clusters_Dune %in% old] <- new[match(pam_clusters_Dune, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
pam_clusters_Dune[pam_clusters_Dune %in% old] <- new[match(pam_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, pam_clusters_Dune)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Dune <- hclust(dist(Dune_scaled_data), method = "ward.D2")
hclust_time_Dune <- toc()
hclust_clusters_Dune <- cutree(hclust_result_Dune, k = G)
table(true_lab_Dune, hclust_clusters_Dune)
new <- c(3,4)
old <- c(4,3)
hclust_clusters_Dune[hclust_clusters_Dune %in% old] <- new[match(hclust_clusters_Dune, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
hclust_clusters_Dune[hclust_clusters_Dune %in% old] <- new[match(hclust_clusters_Dune, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
hclust_clusters_Dune[hclust_clusters_Dune %in% old] <- new[match(hclust_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, hclust_clusters_Dune)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Dune <- Mclust(Dune_scaled_data, G = G)
mclust_time_Dune <- toc()
summary(mclust_result_Dune)
mclust_clusters_Dune <- mclust_result_Dune$classification
table(true_lab_Dune, mclust_clusters_Dune)
new <- c(3,4)
old <- c(4,3)
mclust_clusters_Dune[mclust_clusters_Dune %in% old] <- new[match(mclust_clusters_Dune, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
mclust_clusters_Dune[mclust_clusters_Dune %in% old] <- new[match(mclust_clusters_Dune, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
mclust_clusters_Dune[mclust_clusters_Dune %in% old] <- new[match(mclust_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, mclust_clusters_Dune)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Dune_scaled_data,alphaPriors = c(6,19))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Dune <- toc()
dp$clusterLabels
# Extract clusters 
dpMVN_clusters_Dune <- as.numeric(dp$clusterLabels)
table(true_lab_Dune, dpMVN_clusters_Dune)
new <- c(1,4)
old <- c(4,1)
dpMVN_clusters_Dune[dpMVN_clusters_Dune %in% old] <- new[match(dpMVN_clusters_Dune, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
dpMVN_clusters_Dune[dpMVN_clusters_Dune %in% old] <- new[match(dpMVN_clusters_Dune, old, nomatch = 0)]
new <- c(1,2)
old <- c(2,1)
dpMVN_clusters_Dune[dpMVN_clusters_Dune %in% old] <- new[match(dpMVN_clusters_Dune, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
dpMVN_clusters_Dune[dpMVN_clusters_Dune %in% old] <- new[match(dpMVN_clusters_Dune, old, nomatch = 0)]
table(true_lab_Dune, dpMVN_clusters_Dune)

calculate_clustering_metricsDune_2 <- function(dataPu_2, true_clusters_Dune, estimated_clusters_Dune_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Dune_list)) {
    clusters <- estimated_clusters_Dune_list[[method_name]]
    
    # Calculate metrics
    table_result_Dune <- table(true_clusters_Dune, clusters)
    kappa_result_Dune <- kappa2(data.frame(rater1 = true_clusters_Dune, rater2 = clusters))
    ari_result_Dune <- adjustedRandIndex(true_clusters_Dune, clusters)
    nmi_result_Dune <- NMI(true_clusters_Dune, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Dune,
      Kappa = kappa_result_Dune$value,
      ARI = ari_result_Dune,
      NMI = nmi_result_Dune,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Dune <- list(
  KMeans = kmeans_clusters_Dune,
  CLARA = clara_clusters_Dune,
  PAM = pam_clusters_Dune,
  Hierarchical = hclust_clusters_Dune,
  Mclust = mclust_clusters_Dune,
  DPMVN = dpMVN_clusters_Dune,
  True = true_lab_Dune
)

times_list_Dune <- list(
  KMeans = kmeans_time_Dune,
  CLARA = clara_time_Dune,
  PAM = pam_time_Dune,
  Hierarchical = hclust_time_Dune,
  Mclust = mclust_time_Dune,
  DPMVN = DPMVN_time_Dune
)

# Call the function
clustering_metricsDune_2 <- calculate_clustering_metricsDune_2(dataPu_2 = Dune_scaled_data, true_clusters_Dune = true_lab_Dune, estimated_clusters_Dune_list = cluster_result_Dune, times_list = times_list_Dune)

# Print the results for each method
for (method_name in names(clustering_metricsDune_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsDune_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsDune_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsDune_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsDune_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsDune_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Dune_K <- data.frame(clustering_metricsDune_2$True$Kappa, dp_Gmvlg_metricsDune_1$dp_Gmvlg$Kappa, clustering_metricsDune_2$KMeans$Kappa, clustering_metricsDune_2$CLARA$Kappa, clustering_metricsDune_2$PAM$Kappa, clustering_metricsDune_2$Hierarchical$Kappa, clustering_metricsDune_2$Mclust$Kappa, clustering_metricsDune_2$DPMVN$Kappa)
row_Dune_ARI <- data.frame(clustering_metricsDune_2$True$ARI, dp_Gmvlg_metricsDune_1$dp_Gmvlg$ARI, clustering_metricsDune_2$KMeans$ARI, clustering_metricsDune_2$CLARA$ARI, clustering_metricsDune_2$PAM$ARI, clustering_metricsDune_2$Hierarchical$ARI, clustering_metricsDune_2$Mclust$ARI, clustering_metricsDune_2$DPMVN$ARI)
row_Dune_NMI <- data.frame(clustering_metricsDune_2$True$NMI, dp_Gmvlg_metricsDune_1$dp_Gmvlg$NMI, clustering_metricsDune_2$KMeans$NMI, clustering_metricsDune_2$CLARA$NMI, clustering_metricsDune_2$PAM$NMI, clustering_metricsDune_2$Hierarchical$NMI, clustering_metricsDune_2$Mclust$NMI, clustering_metricsDune_2$DPMVN$NMI)
row_Dune_Cpu <- data.frame(dp_Gmvlg_metricsDune_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsDune_2$KMeans$CPU_RUN_TIME, clustering_metricsDune_2$CLARA$CPU_RUN_TIME, clustering_metricsDune_2$PAM$CPU_RUN_TIME, clustering_metricsDune_2$Hierarchical$CPU_RUN_TIME, clustering_metricsDune_2$Mclust$CPU_RUN_TIME, clustering_metricsDune_2$DPMVN$CPU_RUN_TIME)

tableDune_1 <- clustering_metricsDune_2$True$Table
tableDune_2 <- dp_Gmvlg_metricsDune_1$dp_Gmvlg$Table
tableDune_3 <- clustering_metricsDune_2$KMeans$Table
tableDune_4 <- clustering_metricsDune_2$CLARA$Table
tableDune_5 <- clustering_metricsDune_2$PAM$Table
tableDune_6 <- clustering_metricsDune_2$Hierarchical$Table
tableDune_7 <- clustering_metricsDune_2$Mclust$Table
tableDune_8 <- clustering_metricsDune_2$DPMVN$Table

colnames(row_Dune_K) <- NULL
colnames(row_Dune_ARI) <- NULL
colnames(row_Dune_NMI) <- NULL
colnames(row_Dune_Cpu) <- NULL

row_Dune_K <- as.matrix(row_Dune_K)
row_Dune_ARI <- as.matrix(row_Dune_ARI)
row_Dune_NMI <- as.matrix(row_Dune_NMI)
row_Dune_Cpu <- as.matrix(row_Dune_Cpu)

kappa_table_Dune <- rbind(row_Dune_K)
ARI_table_Dune <- rbind(row_Dune_ARI)
NMI_table_Dune <- rbind(row_Dune_NMI)
cpu_runtime_table_Dune <- rbind(row_Dune_Cpu)
clustering_table_Dune <- rbind(tableDune_1,tableDune_2,tableDune_3,tableDune_4,tableDune_5,tableDune_6,tableDune_7, tableDune_8)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Dune)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Dune)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Dune)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Dune)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Dune)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/02_DUNE DATA/Cluster_metricsDune.xlsx", overwrite = TRUE)
