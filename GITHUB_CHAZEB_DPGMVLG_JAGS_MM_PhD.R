#################################################################
#                       Chaz DATASET                         ##
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
data("chazeb",package ="ade4")
help(chazeb)
meau1 <- data.frame(chazeb$tab,chazeb$cla)
names(meau1)[ncol(meau1)] <- "cla"
Chaz_data <- as.matrix(meau1[, -c(4,7)])
true_lab_Chaz=as.numeric(meau1[,7])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Chaz_data)

# 2. Check for missing values
sum(is.na(Chaz_data))

# 3. Distribution of each feature
par(mfrow=c(2, 3))  # Set up the plotting area for 2 histograms
for (i in 1:ncol(Chaz_data)) {
  hist(Chaz_data[, i], main=colnames(Chaz_data)[i], xlab=colnames(Chaz_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Chaz_data, main="Fruits Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow")[true_lab_Chaz])

# 5. Boxplots for each feature by Diet type
# Add normalized data to the original data frame
Cor_normalized <- as.data.frame(Chaz_data)
Cor_normalized$type <- chazeb$cla

par(mfrow=c(2, 3))  # Reset the plotting area for boxplots
for (i in 1:ncol(Chaz_data)) {
  boxplot(Chaz_data[, i] ~ Cor_normalized$type, main=colnames(Chaz_data)[i], xlab="type", ylab=colnames(Chaz_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Chaz_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the corvus dataset
Chaz_scaled_data <- scale(Chaz_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Chaz_scaled_data, 2, skewness)
kurtosis_values <- apply(Chaz_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Chaz_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Chaz_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Chaz_scaled_data[,i])$out
  print(paste("Feature:", colnames(Chaz_scaled_data)[i]))
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
outliers_matrix <- apply(Chaz_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Chaz_scaled_data, cluster = true_lab_Chaz), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Chaz_scaled_data)
D <- ncol(Chaz_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Chaz_scaled_data,
  #N = N, 
  P = P,
  D = D,
  deltaData = det(cor(Chaz_scaled_data))^(1 / (D - 1)),
  G = G # maximum number of clusters
)

deltaData <- det(cor(Chaz_scaled_data))^(1 / (D - 1))

model_string <- "
Data {
  C <- 1000000000
  for (i in 1:P) {
     zeros[i] <- 0
  }
  v_g <- 1.4262          # uniform prior for conditional
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
    
      z[i, g] <- ifelse(cluster[i] == g, 1, 0) #is an indicator variable that equals 1 if observation i is assigned to cluster g.
   }
   cluster[i] ~ dcat(pi_g) #Each observation i is assigned to one of the G clusters based on the categorical distribution with probabilities pi_g.
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
run_time_Chaz <- toc()
codaSamples_Chaz = as.mcmc.list( runJagsOut )
summaryChains_Chaz <- summary(codaSamples_Chaz)

diagMCMC( codaObject=codaSamples_Chaz , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Chaz , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Chaz$statistics[(1+5*G+5*G+P+1):nrow(summaryChains_Chaz$statistics),1], P, G)
z_mode_Chaz <- apply(matrix(summaryChains_Chaz$statistics[(1+5*G+5*G+P+1):nrow(summaryChains_Chaz$statistics),1], P, G),1, which.max)
z_mode_Chaz

plot(Chaz_scaled_data, col= z_mode_Chaz, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Chaz), col=unique(z_mode_Chaz), pch=16, title="Cluster")
table( true_lab_Chaz , z_mode_Chaz)

# To switch to the same labels in the true clustering
new <- 1:2
old <- c(2,1)
z_mode_Chaz[z_mode_Chaz %in% old] <- new[match(z_mode_Chaz,old,nomatch = 0)]

plot(Chaz_scaled_data, col= z_mode_Chaz, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Chaz), col=unique(z_mode_Chaz), pch=16, title="Cluster")
table( true_lab_Chaz , z_mode_Chaz)
kappa2(data.frame(rater1 = true_lab_Chaz, rater2 = z_mode_Chaz))

calculate_dp_Gmvlg_clustering_metrics <- function(dataChaz_1, true_clusters_Chaz) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Chaz <- table(true_clusters_Chaz, z_mode_Chaz)
  kappa_result_Chaz <- kappa2(data.frame(rater1 = true_clusters_Chaz, rater2 = z_mode_Chaz))
  ari_result_Chaz <- adjustedRandIndex(true_clusters_Chaz, z_mode_Chaz)
  nmi_result_Chaz <- NMI(true_clusters_Chaz, z_mode_Chaz)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Chaz,
    Kappa = kappa_result_Chaz$value,
    ARI = ari_result_Chaz,
    NMI = nmi_result_Chaz,
    CPU_RUN_TIME = run_time_Chaz$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsChaz_1 <- calculate_dp_Gmvlg_clustering_metrics(dataChaz_1 = Chaz_scaled_data, true_clusters_Chaz = true_lab_Chaz)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsChaz_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsChaz_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsChaz_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsChaz_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsChaz_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsChaz_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Chaz <- eclust(Chaz_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Chaz <- toc()
kmeans_clusters_Chaz <- kmeans_result_Chaz$cluster
table(true_lab_Chaz, kmeans_clusters_Chaz)
#new <- 1:2
#old <- c(2,1)
#kmeans_clusters_Chaz[kmeans_clusters_Chaz %in% old] <- new[match(kmeans_clusters_Chaz, old, nomatch = 0)]
table(true_lab_Chaz, kmeans_clusters_Chaz)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Chaz <- eclust(Chaz_scaled_data, "clara", G, graph = FALSE)
clara_time_Chaz <- toc()
clara_clusters_Chaz <- clara_result_Chaz$cluster
table(true_lab_Chaz, clara_clusters_Chaz)
#new <- 1:2
#old <- c(2,1)
#clara_clusters_Chaz[clara_clusters_Chaz %in% old] <- new[match(clara_clusters_Chaz, old, nomatch = 0)]
table(true_lab_Chaz, clara_clusters_Chaz)

# PAM clustering
tic("PAM Runtime")
pam_result_Chaz <- eclust(Chaz_scaled_data, "pam", G, graph = FALSE)
pam_time_Chaz <- toc()
pam_clusters_Chaz <- pam_result_Chaz$cluster
table(true_lab_Chaz, pam_clusters_Chaz)
#new <- 1:2
#old <- c(2,1)
#pam_clusters_Chaz[pam_clusters_Chaz %in% old] <- new[match(pam_clusters_Chaz, old, nomatch = 0)]
table(true_lab_Chaz, pam_clusters_Chaz)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Chaz <- hclust(dist(Chaz_scaled_data), method = "ward.D")
hclust_time_Chaz <- toc()
hclust_clusters_Chaz <- cutree(hclust_result_Chaz, k = G)
table(true_lab_Chaz, hclust_clusters_Chaz)
new <- 1:2
old <- c(2,1)
hclust_clusters_Chaz[hclust_clusters_Chaz %in% old] <- new[match(hclust_clusters_Chaz, old, nomatch = 0)]
table(true_lab_Chaz, hclust_clusters_Chaz)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Chaz <- Mclust(Chaz_scaled_data, G = G)
mclust_time_Chaz <- toc()
summary(mclust_result_Chaz)
mclust_clusters_Chaz <- mclust_result_Chaz$classification
table(true_lab_Chaz, mclust_clusters_Chaz)
new <- 1:2
old <- c(2,1)
mclust_clusters_Chaz[mclust_clusters_Chaz %in% old] <- new[match(mclust_clusters_Chaz, old, nomatch = 0)]
table(true_lab_Chaz, mclust_clusters_Chaz)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Chaz_scaled_data, alphaPriors = c(2,4))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Chaz <- toc()
# Extract clusters 
dpMVN_clusters_Chaz <- as.numeric(dp$clusterLabels)
table(true_lab_Chaz, dpMVN_clusters_Chaz)
new <- 1:2
old <- c(2,1)
dpMVN_clusters_Chaz[dpMVN_clusters_Chaz %in% old] <- new[match(dpMVN_clusters_Chaz, old, nomatch = 0)]
print(dpMVN_clusters_Chaz)
table(true_lab_Chaz, dpMVN_clusters_Chaz)

calculate_clustering_metricsChaz_2 <- function(dataChaz_2, true_clusters_Chaz, estimated_clusters_Chaz_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Chaz_list)) {
    clusters <- estimated_clusters_Chaz_list[[method_name]]
    
    # Calculate metrics
    table_result_Chaz <- table(true_clusters_Chaz, clusters)
    kappa_result_Chaz <- kappa2(data.frame(rater1 = true_clusters_Chaz, rater2 = clusters))
    ari_result_Chaz <- adjustedRandIndex(true_clusters_Chaz, clusters)
    nmi_result_Chaz <- NMI(true_clusters_Chaz, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Chaz,
      Kappa = kappa_result_Chaz$value,
      ARI = ari_result_Chaz,
      NMI = nmi_result_Chaz,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Chaz <- list(
  KMeans = kmeans_clusters_Chaz,
  CLARA = clara_clusters_Chaz,
  PAM = pam_clusters_Chaz,
  Hierarchical = hclust_clusters_Chaz,
  Mclust = mclust_clusters_Chaz,
  DPMVN = dpMVN_clusters_Chaz,
  True = true_lab_Chaz
)

times_list_Chaz <- list(
  KMeans = kmeans_time_Chaz,
  CLARA = clara_time_Chaz,
  PAM = pam_time_Chaz,
  Hierarchical = hclust_time_Chaz,
  Mclust = mclust_time_Chaz,
  DPMVN = DPMVN_time_Chaz
)

# Call the function
clustering_metricsChaz_2 <- calculate_clustering_metricsChaz_2(dataChaz_2 = Chaz_scaled_data, true_clusters_Chaz = true_lab_Chaz, estimated_clusters_Chaz_list = cluster_result_Chaz, times_list = times_list_Chaz)

# Print the results for each method
for (method_name in names(clustering_metricsChaz_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsChaz_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsChaz_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsChaz_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsChaz_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsChaz_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Chaz_K <- data.frame(clustering_metricsChaz_2$True$Kappa, dp_Gmvlg_metricsChaz_1$dp_Gmvlg$Kappa, clustering_metricsChaz_2$KMeans$Kappa, clustering_metricsChaz_2$CLARA$Kappa, clustering_metricsChaz_2$PAM$Kappa, clustering_metricsChaz_2$Hierarchical$Kappa, clustering_metricsChaz_2$Mclust$Kappa, clustering_metricsChaz_2$DPMVN$Kappa)
row_Chaz_ARI <- data.frame(clustering_metricsChaz_2$True$ARI, dp_Gmvlg_metricsChaz_1$dp_Gmvlg$ARI, clustering_metricsChaz_2$KMeans$ARI, clustering_metricsChaz_2$CLARA$ARI, clustering_metricsChaz_2$PAM$ARI, clustering_metricsChaz_2$Hierarchical$ARI, clustering_metricsChaz_2$Mclust$ARI, clustering_metricsChaz_2$DPMVN$ARI)
row_Chaz_NMI <- data.frame(clustering_metricsChaz_2$True$NMI, dp_Gmvlg_metricsChaz_1$dp_Gmvlg$NMI, clustering_metricsChaz_2$KMeans$NMI, clustering_metricsChaz_2$CLARA$NMI, clustering_metricsChaz_2$PAM$NMI, clustering_metricsChaz_2$Hierarchical$NMI, clustering_metricsChaz_2$Mclust$NMI, clustering_metricsChaz_2$DPMVN$NMI)
row_Chaz_Cpu <- data.frame(dp_Gmvlg_metricsChaz_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsChaz_2$KMeans$CPU_RUN_TIME, clustering_metricsChaz_2$CLARA$CPU_RUN_TIME, clustering_metricsChaz_2$PAM$CPU_RUN_TIME, clustering_metricsChaz_2$Hierarchical$CPU_RUN_TIME, clustering_metricsChaz_2$Mclust$CPU_RUN_TIME, clustering_metricsChaz_2$DPMVN$CPU_RUN_TIME)

tableChaz_1 <- clustering_metricsChaz_2$True$Table
tableChaz_2 <- dp_Gmvlg_metricsChaz_1$dp_Gmvlg$Table
tableChaz_3 <- clustering_metricsChaz_2$KMeans$Table
tableChaz_4 <- clustering_metricsChaz_2$CLARA$Table
tableChaz_5 <- clustering_metricsChaz_2$PAM$Table
tableChaz_6 <- clustering_metricsChaz_2$Hierarchical$Table
tableChaz_7 <- clustering_metricsChaz_2$Mclust$Table

colnames(row_Chaz_K) <- NULL
colnames(row_Chaz_ARI) <- NULL
colnames(row_Chaz_NMI) <- NULL
colnames(row_Chaz_Cpu) <- NULL

row_Chaz_K <- as.matrix(row_Chaz_K)
row_Chaz_ARI <- as.matrix(row_Chaz_ARI)
row_Chaz_NMI <- as.matrix(row_Chaz_NMI)
row_Chaz_Cpu <- as.matrix(row_Chaz_Cpu)

kappa_table_Chaz <- rbind(row_Chaz_K)
ARI_table_Chaz <- rbind(row_Chaz_ARI)
NMI_table_Chaz <- rbind(row_Chaz_NMI)
cpu_runtime_table_Chaz <- rbind(row_Chaz_Cpu)
clustering_table_Chaz <- rbind(tableChaz_1,tableChaz_2,tableChaz_3,tableChaz_4,tableChaz_5,tableChaz_6,tableChaz_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Chaz)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Chaz)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Chaz)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Chaz)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Chaz)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/04_CHAZEB/Cluster_metricsChaz.xlsx", overwrite = TRUE)
