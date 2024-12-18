#################################################################
#                        CATS DATASET                          ##
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
library(moments)
library(aricode)

# Load the Cats dataset
data("cats")
help(cats)
Cat_data <- as.matrix(cats[,-1]) # Convert all columns to numeric
true_lab_Cat <- as.numeric(cats$Sex)    # Use the 'class' column as the true labels

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(cat_data)

# 2. Check for missing values
sum(is.na(cat_data))

# 3. Distribution of each feature
par(mfrow=c(1, 2))  # Adjusted for the number of features (cats dataset has 4 features)
for (i in 1:ncol(cat_data)) {
  hist(cat_data[, i], main=colnames(cat_data)[i], xlab=colnames(cat_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(cat_data, main="Cats Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Cat])

# 5. Boxplots for each feature by 'Sex'
# Add normalized data to the original data frame
cats_normalized <- as.data.frame(cat_data)
cats_normalized$Sex <- as.factor(cats$Sex)

par(mfrow=c(1, 2))  # Reset the plotting area for boxplots
for (i in 1:ncol(cat_data)) {
  boxplot(cat_data[, i] ~ cats_normalized$Sex, main=colnames(cat_data)[i], xlab="Sex", ylab=colnames(cat_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(cat_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Cats dataset
Cat_scaled_data <- scale(Cat_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Cat_scaled_data, 2, skewness)
kurtosis_values <- apply(Cat_scaled_data, 2, kurtosis)

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
combined_data <- as.vector(as.matrix(Cat_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Cat_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Cat_scaled_data[,i])$out
  print(paste("Feature:", colnames(Cat_scaled_data)[i]))
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
outliers_matrix <- apply(Cat_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Cat_scaled_data, cluster = true_lab_Cat), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Cat_scaled_data)
D <- ncol(Cat_scaled_data)

#Try with different number of clusters
G <- 2  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Cat_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Cat_scaled_data))^(1 / (D - 1)),
  G = G 
)

deltaData <- det(cor(Cat_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(10,13)        
       lambda_g[j, g] ~ dgamma(13,13)
    }
   }
   alpha[1] ~ dgamma(20,22)
   alpha[2] ~ dgamma(23,24)
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
run_time_Cat <- toc()
codaSamples_Cat = as.mcmc.list( runJagsOut )
summaryChains_Cat <- summary(codaSamples_Cat)

diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[2,2]" )
graphics.off()

matrix(summaryChains_Cat$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Cat$statistics),1], P, G)
z_mode_Cat <- apply(matrix(summaryChains_Cat$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Cat$statistics),1], P, G),1, which.max)
z_mode_Cat

plot(Cat_scaled_data, col= z_mode_Cat, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Cat), col=unique(z_mode_Cat), pch=16, title="Cluster")
table( true_lab_Cat , z_mode_Cat)

# To switch to the same labels in the true clustering
new <- 1:2
old <- c(2,1)
z_mode_Cat[z_mode_Cat %in% old] <- new[match(z_mode_Cat,old,nomatch = 0)]

plot(Cat_scaled_data, col= z_mode_Cat, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Cat), col=unique(z_mode_Cat), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Cat , z_mode_Cat)
kappa2(data.frame(rater1 = true_lab_Cat, rater2 = z_mode_Cat))

calculate_dp_Gmvlg_clustering_metrics <- function(dataCat_1, true_clusters_Cat) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Cat <- table(true_clusters_Cat, z_mode_Cat)
  kappa_result_Cat <- kappa2(data.frame(rater1 = true_clusters_Cat, rater2 = z_mode_Cat))
  ari_result_Cat <- adjustedRandIndex(true_clusters_Cat, z_mode_Cat)
  nmi_result_Cat <- NMI(true_clusters_Cat, z_mode_Cat)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Cat,
    Kappa = kappa_result_Cat$value,
    ARI = ari_result_Cat,
    NMI = nmi_result_Cat,
    CPU_RUN_TIME = run_time_Cat$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsCat_1 <- calculate_dp_Gmvlg_clustering_metrics(dataCat_1 = Cat_scaled_data, true_clusters_Cat = true_lab_Cat)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsCat_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Cat <- eclust(Cat_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Cat <- toc()
kmeans_clusters_Cat <- kmeans_result_Cat$cluster
table(true_lab_Cat, kmeans_clusters_Cat)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Cat <- eclust(Cat_scaled_data, "clara", G, graph = FALSE)
clara_time_Cat <- toc()
clara_clusters_Cat <- clara_result_Cat$cluster
table(true_lab_Cat, clara_clusters_Cat)

# PAM clustering
tic("PAM Runtime")
pam_result_Cat <- eclust(Cat_scaled_data, "pam", G, graph = FALSE)
pam_time_Cat <- toc()
pam_clusters_Cat <- pam_result_Cat$cluster
table(true_lab_Cat, pam_clusters_Cat)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Cat <- hclust(dist(Cat_scaled_data), method = "ward.D2")
hclust_time_Cat <- toc()
hclust_clusters_Cat <- cutree(hclust_result_Cat, k = G)
table(true_lab_Cat, hclust_clusters_Cat)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Cat <- Mclust(Cat_scaled_data, G = G)
mclust_time_Cat <- toc()
summary(mclust_result_Cat)
mclust_clusters_Cat <- mclust_result_Cat$classification
table(true_lab_Cat, mclust_clusters_Cat)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Cat_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Cat <- toc()
# Extract clusters 
dpMVN_clusters_Cat <- as.numeric(dp$clusterLabels)
new <- 1:2
old <- c(2,1)
dpMVN_clusters_Cat[dpMVN_clusters_Cat %in% old] <- new[match(dpMVN_clusters_Cat, old, nomatch = 0)]
print(dpMVN_clusters_Cat)
table(true_lab_Cat, dpMVN_clusters_Cat)

calculate_clustering_metricsCat_2 <- function(dataCat_2, true_clusters_Cat, estimated_clusters_Cat_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Cat_list)) {
    clusters <- estimated_clusters_Cat_list[[method_name]]
    
    # Calculate metrics
    table_result_Cat <- table(true_clusters_Cat, clusters)
    kappa_result_Cat <- kappa2(data.frame(rater1 = true_clusters_Cat, rater2 = clusters))
    ari_result_Cat <- adjustedRandIndex(true_clusters_Cat, clusters)
    nmi_result_Cat <- NMI(true_clusters_Cat, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Cat,
      Kappa = kappa_result_Cat$value,
      ARI = ari_result_Cat,
      NMI = nmi_result_Cat,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Cat <- list(
  KMeans = kmeans_clusters_Cat,
  CLARA = clara_clusters_Cat,
  PAM = pam_clusters_Cat,
  Hierarchical = hclust_clusters_Cat,
  Mclust = mclust_clusters_Cat,
  DPMVN = dpMVN_clusters_Cat,
  True = true_lab_Cat
)

times_list_Cat <- list(
  KMeans = kmeans_time_Cat,
  CLARA = clara_time_Cat,
  PAM = pam_time_Cat,
  Hierarchical = hclust_time_Cat,
  Mclust = mclust_time_Cat,
  DPMVN = DPMVN_time_Cat
)

# Call the function
clustering_metricsCat_2 <- calculate_clustering_metricsCat_2(dataCat_2 = Cat_scaled_data, true_clusters_Cat = true_lab_Cat, estimated_clusters_Cat_list = cluster_result_Cat, times_list = times_list_Cat)

# Print the results for each method
for (method_name in names(clustering_metricsCat_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsCat_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsCat_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsCat_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsCat_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsCat_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Cat_K <- data.frame(clustering_metricsCat_2$True$Kappa, dp_Gmvlg_metricsCat_1$dp_Gmvlg$Kappa, clustering_metricsCat_2$KMeans$Kappa, clustering_metricsCat_2$CLARA$Kappa, clustering_metricsCat_2$PAM$Kappa, clustering_metricsCat_2$Hierarchical$Kappa, clustering_metricsCat_2$Mclust$Kappa, clustering_metricsCat_2$DPMVN$Kappa)
row_Cat_ARI <- data.frame(clustering_metricsCat_2$True$ARI, dp_Gmvlg_metricsCat_1$dp_Gmvlg$ARI, clustering_metricsCat_2$KMeans$ARI, clustering_metricsCat_2$CLARA$ARI, clustering_metricsCat_2$PAM$ARI, clustering_metricsCat_2$Hierarchical$ARI, clustering_metricsCat_2$Mclust$ARI, clustering_metricsCat_2$DPMVN$ARI)
row_Cat_NMI <- data.frame(clustering_metricsCat_2$True$NMI, dp_Gmvlg_metricsCat_1$dp_Gmvlg$NMI, clustering_metricsCat_2$KMeans$NMI, clustering_metricsCat_2$CLARA$NMI, clustering_metricsCat_2$PAM$NMI, clustering_metricsCat_2$Hierarchical$NMI, clustering_metricsCat_2$Mclust$NMI, clustering_metricsCat_2$DPMVN$NMI)
row_Cat_Cpu <- data.frame(dp_Gmvlg_metricsCat_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsCat_2$KMeans$CPU_RUN_TIME, clustering_metricsCat_2$CLARA$CPU_RUN_TIME, clustering_metricsCat_2$PAM$CPU_RUN_TIME, clustering_metricsCat_2$Hierarchical$CPU_RUN_TIME, clustering_metricsCat_2$Mclust$CPU_RUN_TIME, clustering_metricsCat_2$DPMVN$CPU_RUN_TIME)

tableCat_1 <- clustering_metricsCat_2$True$Table
tableCat_2 <- dp_Gmvlg_metricsCat_1$dp_Gmvlg$Table
tableCat_3 <- clustering_metricsCat_2$KMeans$Table
tableCat_4 <- clustering_metricsCat_2$CLARA$Table
tableCat_5 <- clustering_metricsCat_2$PAM$Table
tableCat_6 <- clustering_metricsCat_2$Hierarchical$Table
tableCat_7 <- clustering_metricsCat_2$Mclust$Table
tableCat_8 <- clustering_metricsCat_2$DPMVN$Table

colnames(row_Cat_K) <- NULL
colnames(row_Cat_ARI) <- NULL
colnames(row_Cat_NMI) <- NULL
colnames(row_Cat_Cpu) <- NULL

row_Cat_K <- as.matrix(row_Cat_K)
row_Cat_ARI <- as.matrix(row_Cat_ARI)
row_Cat_NMI <- as.matrix(row_Cat_NMI)
row_Cat_Cpu <- as.matrix(row_Cat_Cpu)

kappa_table_Cat <- rbind(row_Cat_K)
ARI_table_Cat <- rbind(row_Cat_ARI)
NMI_table_Cat <- rbind(row_Cat_NMI)
cpu_runtime_table_Cat <- rbind(row_Cat_Cpu)
clustering_table_Cat <- rbind(tableCat_1,tableCat_2,tableCat_3,tableCat_4,tableCat_5,tableCat_6,tableCat_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Cat)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Cat)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Cat)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Cat)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Cat)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/13_CAT/Cluster_metricsCat.xlsx", overwrite = TRUE)
