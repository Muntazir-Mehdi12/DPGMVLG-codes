#################################################################
#                       Elisa DATASET                          ##
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

# Load the dataset
data("elisa",package ="isdals")
help("elisa",package ="isdals")
Elisa_data <- as.matrix(elisa[, -c(1)])
true_lab_Elisa=as.numeric(elisa[,1])

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Elisa_data)

# 2. Check for missing values
sum(is.na(Elisa_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Set up the plotting area for 9 histograms
for (i in 1:2) {
  hist(Elisa_data[, i], main=colnames(Elisa_data)[i], xlab=colnames(Elisa_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Elisa_data, main="Elisa Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue", "yellow", "purple", "orange", "pink")[true_lab_Elisa])

# 5. Boxplots for each feature by 'Type'
# Add normalized data to the original data frame
Elisa_normalized <- as.data.frame(Elisa_data)
Elisa_normalized$Type <- elisa$type

par(mfrow=c(1, 2))  # Reset the plotting area for boxplots
for (i in 1:2) {
  boxplot(Elisa_data[, i] ~ elisa$type, main=colnames(Elisa_data)[i], xlab="Type", ylab=colnames(Elisa_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Elisa_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the dataset
Elisa_scaled_data <- scale(Elisa_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Elisa_scaled_data, 2, skewness)
kurtosis_values <- apply(Elisa_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Elisa_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Elisa_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Elisa_scaled_data[,i])$out
  print(paste("Feature:", colnames(Elisa_scaled_data)[i]))
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
outliers_matrix <- apply(Elisa_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Elisa_scaled_data, cluster = true_lab_Elisa), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Elisa_scaled_data)
D <- ncol(Elisa_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Elisa_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Elisa_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Elisa_scaled_data))^(1 / (D - 1))

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
       lambda_g[j, g] ~ dgamma(1,6)     
    }
    alpha[g] ~ dgamma(2,6)
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
run_time_Elisa <- toc()
codaSamples_Elisa = as.mcmc.list( runJagsOut )
summaryChains_Elisa <- summary(codaSamples_Elisa)

diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[2,2]" )
graphics.off()

z_mode_Elisa <- apply(matrix(summaryChains_Elisa$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Elisa$statistics),1], P, G),1, which.max)
z_mode_Elisa

plot(Elisa_scaled_data, col= z_mode_Elisa, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Elisa), col=unique(z_mode_Elisa), pch=16, title="Cluster")
table( true_lab_Elisa , z_mode_Elisa)

calculate_dp_Gmvlg_clustering_metrics <- function(dataElisa_1, true_clusters_Elisa) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Elisa <- table(true_clusters_Elisa, z_mode_Elisa)
  kappa_result_Elisa <- kappa2(data.frame(rater1 = true_clusters_Elisa, rater2 = z_mode_Elisa))
  ari_result_Elisa <- adjustedRandIndex(true_clusters_Elisa, z_mode_Elisa)
  nmi_result_Elisa <- NMI(true_clusters_Elisa, z_mode_Elisa)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Elisa,
    Kappa = kappa_result_Elisa$value,
    ARI = ari_result_Elisa,
    NMI = nmi_result_Elisa,
    CPU_RUN_TIME = run_time_Elisa$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsElisa_1 <- calculate_dp_Gmvlg_clustering_metrics(dataElisa_1 = Elisa_scaled_data, true_clusters_Elisa = true_lab_Elisa)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsElisa_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsElisa_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsElisa_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsElisa_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsElisa_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsElisa_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Elisa <- eclust(Elisa_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Elisa <- toc()
kmeans_clusters_Elisa <- kmeans_result_Elisa$cluster
table(true_lab_Elisa, kmeans_clusters_Elisa)
#Label switching
new <- 1:2
old <- c(2,1)
kmeans_clusters_Elisa[kmeans_clusters_Elisa %in% old] <- new[match(kmeans_clusters_Elisa, old, nomatch = 0)]
table(true_lab_Elisa, kmeans_clusters_Elisa)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Elisa <- eclust(Elisa_scaled_data, "clara", G, graph = FALSE)
clara_time_Elisa <- toc()
clara_clusters_Elisa <- clara_result_Elisa$cluster
table(true_lab_Elisa, clara_clusters_Elisa)

# PAM clustering
tic("PAM Runtime")
pam_result_Elisa <- eclust(Elisa_scaled_data, "pam", G, graph = FALSE)
pam_time_Elisa <- toc()
pam_clusters_Elisa <- pam_result_Elisa$cluster
table(true_lab_Elisa, pam_clusters_Elisa)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Elisa <- hclust(dist(Elisa_scaled_data), method = "ward.D2")
hclust_time_Elisa <- toc()
hclust_clusters_Elisa <- cutree(hclust_result_Elisa, k = G)
table(true_lab_Elisa, hclust_clusters_Elisa)
new <- 1:2
old <- c(2,1)
hclust_clusters_Elisa[hclust_clusters_Elisa %in% old] <- new[match(hclust_clusters_Elisa, old, nomatch = 0)]
table(true_lab_Elisa, hclust_clusters_Elisa)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Elisa <- Mclust(Elisa_scaled_data, G = G)
mclust_time_Elisa <- toc()
summary(mclust_result_Elisa)
mclust_clusters_Elisa <- mclust_result_Elisa$classification
table(true_lab_Elisa, mclust_clusters_Elisa)
new <- 1:2
old <- c(2,1)
mclust_clusters_Elisa[mclust_clusters_Elisa %in% old] <- new[match(mclust_clusters_Elisa, old, nomatch = 0)]
table(true_lab_Elisa, mclust_clusters_Elisa)

# Bayesian Clustering using Dirichlet Process Multivariate normal distributions (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Elisa_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Elisa <- toc()
# Extract clusters 
dpMVN_clusters_Elisa <- as.numeric(dp$clusterLabels)
table(true_lab_Elisa, dpMVN_clusters_Elisa)

calculate_clustering_metricsElisa_2 <- function(dataElisa_2, true_clusters_Elisa, estimated_clusters_Elisa_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Elisa_list)) {
    clusters <- estimated_clusters_Elisa_list[[method_name]]
    
    # Calculate metrics
    table_result_Elisa <- table(true_clusters_Elisa, clusters)
    kappa_result_Elisa <- kappa2(data.frame(rater1 = true_clusters_Elisa, rater2 = clusters))
    ari_result_Elisa <- adjustedRandIndex(true_clusters_Elisa, clusters)
    nmi_result_Elisa <- NMI(true_clusters_Elisa, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Elisa,
      Kappa = kappa_result_Elisa$value,
      ARI = ari_result_Elisa,
      NMI = nmi_result_Elisa,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Elisa <- list(
  KMeans = kmeans_clusters_Elisa,
  CLARA = clara_clusters_Elisa,
  PAM = pam_clusters_Elisa,
  Hierarchical = hclust_clusters_Elisa,
  Mclust = mclust_clusters_Elisa,
  DPMVN = dpMVN_clusters_Elisa,
  True = true_lab_Elisa
)

times_list_Elisa <- list(
  KMeans = kmeans_time_Elisa,
  CLARA = clara_time_Elisa,
  PAM = pam_time_Elisa,
  Hierarchical = hclust_time_Elisa,
  Mclust = mclust_time_Elisa,
  DPMVN = DPMVN_time_Elisa
)

# Call the function
clustering_metricsElisa_2 <- calculate_clustering_metricsElisa_2(dataElisa_2 = Elisa_scaled_data, true_clusters_Elisa = true_lab_Elisa, estimated_clusters_Elisa_list = cluster_result_Elisa, times_list = times_list_Elisa)

# Print the results for each method
for (method_name in names(clustering_metricsElisa_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsElisa_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsElisa_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsElisa_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsElisa_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsElisa_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Elisa_K <- data.frame(clustering_metricsElisa_2$True$Kappa, dp_Gmvlg_metricsElisa_1$dp_Gmvlg$Kappa, clustering_metricsElisa_2$KMeans$Kappa, clustering_metricsElisa_2$CLARA$Kappa, clustering_metricsElisa_2$PAM$Kappa, clustering_metricsElisa_2$Hierarchical$Kappa, clustering_metricsElisa_2$Mclust$Kappa, clustering_metricsElisa_2$DPMVN$Kappa)
row_Elisa_ARI <- data.frame(clustering_metricsElisa_2$True$ARI, dp_Gmvlg_metricsElisa_1$dp_Gmvlg$ARI, clustering_metricsElisa_2$KMeans$ARI, clustering_metricsElisa_2$CLARA$ARI, clustering_metricsElisa_2$PAM$ARI, clustering_metricsElisa_2$Hierarchical$ARI, clustering_metricsElisa_2$Mclust$ARI, clustering_metricsElisa_2$DPMVN$ARI)
row_Elisa_NMI <- data.frame(clustering_metricsElisa_2$True$NMI, dp_Gmvlg_metricsElisa_1$dp_Gmvlg$NMI, clustering_metricsElisa_2$KMeans$NMI, clustering_metricsElisa_2$CLARA$NMI, clustering_metricsElisa_2$PAM$NMI, clustering_metricsElisa_2$Hierarchical$NMI, clustering_metricsElisa_2$Mclust$NMI, clustering_metricsElisa_2$DPMVN$NMI)
row_Elisa_Cpu <- data.frame(dp_Gmvlg_metricsElisa_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsElisa_2$KMeans$CPU_RUN_TIME, clustering_metricsElisa_2$CLARA$CPU_RUN_TIME, clustering_metricsElisa_2$PAM$CPU_RUN_TIME, clustering_metricsElisa_2$Hierarchical$CPU_RUN_TIME, clustering_metricsElisa_2$Mclust$CPU_RUN_TIME, clustering_metricsElisa_2$DPMVN$CPU_RUN_TIME)

tableElisa_1 <- clustering_metricsElisa_2$True$Table
tableElisa_2 <- dp_Gmvlg_metricsElisa_1$dp_Gmvlg$Table
tableElisa_3 <- clustering_metricsElisa_2$KMeans$Table
tableElisa_4 <- clustering_metricsElisa_2$CLARA$Table
tableElisa_5 <- clustering_metricsElisa_2$PAM$Table
tableElisa_6 <- clustering_metricsElisa_2$Hierarchical$Table
tableElisa_7 <- clustering_metricsElisa_2$Mclust$Table
tableElisa_8 <- clustering_metricsElisa_2$DPMVN$Table

colnames(row_Elisa_K) <- NULL
colnames(row_Elisa_ARI) <- NULL
colnames(row_Elisa_NMI) <- NULL
colnames(row_Elisa_Cpu) <- NULL

row_Elisa_K <- as.matrix(row_Elisa_K)
row_Elisa_ARI <- as.matrix(row_Elisa_ARI)
row_Elisa_NMI <- as.matrix(row_Elisa_NMI)
row_Elisa_Cpu <- as.matrix(row_Elisa_Cpu)

kappa_table_Elisa <- rbind(row_Elisa_K)
ARI_table_Elisa <- rbind(row_Elisa_ARI)
NMI_table_Elisa <- rbind(row_Elisa_NMI)
cpu_runtime_table_Elisa <- rbind(row_Elisa_Cpu)
clustering_table_Elisa <- rbind(tableElisa_1,tableElisa_2,tableElisa_3,tableElisa_4,tableElisa_5,tableElisa_6,tableElisa_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Elisa)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Elisa)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Elisa)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Elisa)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Elisa)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/01_Elisa/Cluster_metricsElisa.xlsx", overwrite = TRUE)
