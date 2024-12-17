#################################################################
#                       Ecomor DATASET                         ##
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
data(ecomor, package = "ade4")
ecomor                            
help(ecomor)
str(ecomor)
uniques <- lapply(ecomor, unique)
uniques
Ecomor_data <- as.matrix(ecomor$morpho)
true_lab_Ecomor=as.numeric(ecomor$taxo$Ordre)

###########################################################################################################
#                                          Descriptive analysis
###########################################################################################################

# 1. Summary statistics
summary(Ecomor_data)

# 2. Check for missing values
sum(is.na(Ecomor_data))

# 3. Distribution of each feature
par(mfrow=c(2, 4))  # Set up the plotting area for 8 histograms (adjusted for the number of features)
for (i in 1:ncol(Ecomor_data)) {
  hist(Ecomor_data[, i], main=colnames(Ecomor_data)[i], xlab=colnames(Ecomor_data)[i], col="lightblue", border="black")
}

# 4. Pairwise scatter plots
pairs(Ecomor_data, main="ECOMOR Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Ecomor])

# 5. Boxplots for each feature by 'diabetes' status
# Add normalized data to the original data frame
pid_normalized <- as.data.frame(Ecomor_data)
pid_normalized$diabetes <- ecomor$taxo$Ordre

par(mfrow=c(2, 4))  # Reset the plotting area for boxplots
for (i in 1:ncol(Ecomor_data)) {
  boxplot(Ecomor_data[, i] ~ ecomor$taxo$Ordre, main=colnames(Ecomor_data)[i], xlab="Order", ylab=colnames(Ecomor_data)[i], col="lightblue")
}

# 6. Correlation matrix
cor_matrix <- cor(Ecomor_data)
print(cor_matrix)

# 7. Visualizing correlations with heatmap
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

# Normalize the Corvus dataset
Ecomor_scaled_data <- scale(Ecomor_data)

# 8. Calculate skewness and kurtosis for each feature
skewness_values <- apply(Ecomor_scaled_data, 2, skewness)
kurtosis_values <- apply(Ecomor_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

par(mfrow=c(1, 2))
# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=2)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=2)

#Combined skewness and kurtosis
combined_data <- as.vector(as.matrix(Ecomor_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#Outliers
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(Ecomor_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Ecomor_scaled_data[,i])$out
  print(paste("Feature:", colnames(Ecomor_scaled_data)[i]))
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
outliers_matrix <- apply(Ecomor_scaled_data, 2, detect_outlier)
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
fviz_cluster(list(data = Ecomor_scaled_data, cluster = true_lab_Ecomor), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

# Define the number of rows and columns
P <- nrow(Ecomor_scaled_data)
D <- ncol(Ecomor_scaled_data)

#Try with different number of clusters
G <- 7  #maximum no. of clusters

# Set up data for JAGS
data_list <- list(
  x = Ecomor_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Ecomor_scaled_data))^(1 / (D - 1)),
  G = G
)

deltaData <- det(cor(Ecomor_scaled_data))^(1 / (D - 1))

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
       mu_g[j, g] ~ dgamma(11,13)         
       lambda_g[j, g] ~ dgamma(1,4)
    }
    alpha[g] <- 10
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

matrix(summaryChains_carni$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_carni$statistics),1], P, G)
z_mode_Ecomor <- apply(matrix(summaryChains_carni$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_carni$statistics),1], P, G),1, which.max)
z_mode_Ecomor

plot(Ecomor_scaled_data, col= z_mode_Ecomor, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ecomor), col=unique(z_mode_Ecomor), pch=16, title="Cluster")
table( true_lab_Ecomor , z_mode_Ecomor)

# To switch to the same labels in the true clustering
new <- c(2,5)
old <- c(5,2)
z_mode_Ecomor[z_mode_Ecomor %in% old] <- new[match(z_mode_Ecomor,old,nomatch = 0)]
new <- c(1,7)
old <- c(7,1)
z_mode_Ecomor[z_mode_Ecomor %in% old] <- new[match(z_mode_Ecomor,old,nomatch = 0)]
new <- c(3,4)
old <- c(4,3)
z_mode_Ecomor[z_mode_Ecomor %in% old] <- new[match(z_mode_Ecomor,old,nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
z_mode_Ecomor[z_mode_Ecomor %in% old] <- new[match(z_mode_Ecomor,old,nomatch = 0)]

plot(Ecomor_scaled_data, col= z_mode_Ecomor, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ecomor), col=unique(z_mode_Ecomor), pch=16, title="Cluster")

# To switch to the same labels in the true clustering
table( true_lab_Ecomor , z_mode_Ecomor)
kappa2(data.frame(rater1 = true_lab_Ecomor, rater2 = z_mode_Ecomor))

calculate_dp_Gmvlg_clustering_metrics <- function(dataEcomor_1, true_clusters_Ecomor) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Now use z_mode to calculate clustering metrics
  
  # Calculate clustering metrics and the contingency table
  table_result_Ecomor <- table(true_clusters_Ecomor, z_mode_Ecomor)
  kappa_result_Ecomor <- kappa2(data.frame(rater1 = true_clusters_Ecomor, rater2 = z_mode_Ecomor))
  ari_result_Ecomor <- adjustedRandIndex(true_clusters_Ecomor, z_mode_Ecomor)
  nmi_result_Ecomor <- NMI(true_clusters_Ecomor, z_mode_Ecomor)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result_Ecomor,
    Kappa = kappa_result_Ecomor$value,
    ARI = ari_result_Ecomor,
    NMI = nmi_result_Ecomor,
    CPU_RUN_TIME = run_time_carni$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

# Call the function
dp_Gmvlg_metricsEcomor_1 <- calculate_dp_Gmvlg_clustering_metrics(dataEcomor_1 = Ecomor_scaled_data, true_clusters_Ecomor = true_lab_Ecomor)

# Print the results for each method
for (dp_Gmvlg in names(dp_Gmvlg_metricsEcomor_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

######################################################################################################################
#                                 All Classical methods and One Bayesian method (DPMVN)
######################################################################################################################

# K-means clustering
set.seed(123)
tic("K-means Runtime")
kmeans_result_Ecomor <- eclust(Ecomor_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Ecomor <- toc()
kmeans_clusters_Ecomor <- kmeans_result_Ecomor$cluster
table(true_lab_Ecomor, kmeans_clusters_Ecomor)
new <- c(5,7)
old <- c(7,5)
kmeans_clusters_Ecomor[kmeans_clusters_Ecomor %in% old] <- new[match(kmeans_clusters_Ecomor, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
kmeans_clusters_Ecomor[kmeans_clusters_Ecomor %in% old] <- new[match(kmeans_clusters_Ecomor, old, nomatch = 0)]
new <- c(1,7)
old <- c(7,1)
kmeans_clusters_Ecomor[kmeans_clusters_Ecomor %in% old] <- new[match(kmeans_clusters_Ecomor, old, nomatch = 0)]
new <- c(3,7)
old <- c(7,3)
kmeans_clusters_Ecomor[kmeans_clusters_Ecomor %in% old] <- new[match(kmeans_clusters_Ecomor, old, nomatch = 0)]
table(true_lab_Ecomor, kmeans_clusters_Ecomor)

# CLARA clustering
tic("CLARA Runtime")
clara_result_Ecomor <- eclust(Ecomor_scaled_data, "clara", G, graph = FALSE)
clara_time_Ecomor <- toc()
clara_clusters_Ecomor <- clara_result_Ecomor$cluster
table(true_lab_Ecomor, clara_clusters_Ecomor)
new <- c(5,6)
old <- c(6,5)
clara_clusters_Ecomor[clara_clusters_Ecomor %in% old] <- new[match(clara_clusters_Ecomor, old, nomatch = 0)]
new <- c(3,6)
old <- c(6,3)
clara_clusters_Ecomor[clara_clusters_Ecomor %in% old] <- new[match(clara_clusters_Ecomor, old, nomatch = 0)]
table(true_lab_Ecomor, clara_clusters_Ecomor)


# PAM clustering
tic("PAM Runtime")
pam_result_Ecomor <- eclust(Ecomor_scaled_data, "pam", G, graph = FALSE)
pam_time_Ecomor <- toc()
pam_clusters_Ecomor <- pam_result_Ecomor$cluster
table(true_lab_Ecomor, pam_clusters_Ecomor)

# Hierarchical clustering
tic("Hierarchical Runtime")
hclust_result_Ecomor <- hclust(dist(Ecomor_scaled_data), method = "ward.D2")
hclust_time_Ecomor <- toc()
hclust_clusters_Ecomor <- cutree(hclust_result_Ecomor, k = G)
table(true_lab_Ecomor, hclust_clusters_Ecomor)
new <- c(3,4)
old <- c(4,3)
hclust_clusters_Ecomor[hclust_clusters_Ecomor %in% old] <- new[match(hclust_clusters_Ecomor, old, nomatch = 0)]
new <- c(5,6)
old <- c(6,5)
hclust_clusters_Ecomor[hclust_clusters_Ecomor %in% old] <- new[match(hclust_clusters_Ecomor, old, nomatch = 0)]
table(true_lab_Ecomor, hclust_clusters_Ecomor)

# Model-based clustering
tic("Model-based Runtime")
mclust_result_Ecomor <- Mclust(Ecomor_scaled_data, G = G)
mclust_time_Ecomor <- toc()
summary(mclust_result_Ecomor)
mclust_clusters_Ecomor <- mclust_result_Ecomor$classification
table(true_lab_Ecomor, mclust_clusters_Ecomor)
new <- c(4,7)
old <- c(7,4)
mclust_clusters_Ecomor[mclust_clusters_Ecomor %in% old] <- new[match(mclust_clusters_Ecomor, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
mclust_clusters_Ecomor[mclust_clusters_Ecomor %in% old] <- new[match(mclust_clusters_Ecomor, old, nomatch = 0)]
new <- c(2,4)
old <- c(4,2)
mclust_clusters_Ecomor[mclust_clusters_Ecomor %in% old] <- new[match(mclust_clusters_Ecomor, old, nomatch = 0)]
table(true_lab_Ecomor, mclust_clusters_Ecomor)

# Bayesian Clustering using Dirichlet Process Gaussian Mixture Model (DPMVN)
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Ecomor_scaled_data, alphaPriors =  c(2,50))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Ecomor <- toc()
dp
# Extract clusters 
dpMVN_clusters_Ecomor <- as.numeric(dp$clusterLabels)
new <- c(1,5)
old <- c(5,1)
dpMVN_clusters_Ecomor[dpMVN_clusters_Ecomor %in% old] <- new[match(dpMVN_clusters_Ecomor, old, nomatch = 0)]
new <- c(1,3)
old <- c(3,1)
dpMVN_clusters_Ecomor[dpMVN_clusters_Ecomor %in% old] <- new[match(dpMVN_clusters_Ecomor, old, nomatch = 0)]
new <- c(2,9)
old <- c(9,2)
dpMVN_clusters_Ecomor[dpMVN_clusters_Ecomor %in% old] <- new[match(dpMVN_clusters_Ecomor, old, nomatch = 0)]
new <- c(3,13)
old <- c(13,3)
dpMVN_clusters_Ecomor[dpMVN_clusters_Ecomor %in% old] <- new[match(dpMVN_clusters_Ecomor, old, nomatch = 0)]
new <- c(7,9)
old <- c(9,7)
dpMVN_clusters_Ecomor[dpMVN_clusters_Ecomor %in% old] <- new[match(dpMVN_clusters_Ecomor, old, nomatch = 0)]
table(true_lab_Ecomor, dpMVN_clusters_Ecomor)

calculate_clustering_metricsEcomor_2 <- function(dataEcomor_2, true_clusters_Ecomor, estimated_clusters_Ecomor_list, times_list) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each set of estimated clusters
  for (method_name in names(estimated_clusters_Ecomor_list)) {
    clusters <- estimated_clusters_Ecomor_list[[method_name]]
    
    # Calculate metrics
    table_result_Ecomor <- table(true_clusters_Ecomor, clusters)
    kappa_result_Ecomor <- kappa2(data.frame(rater1 = true_clusters_Ecomor, rater2 = clusters))
    ari_result_Ecomor <- adjustedRandIndex(true_clusters_Ecomor, clusters)
    nmi_result_Ecomor <- NMI(true_clusters_Ecomor, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_Ecomor,
      Kappa = kappa_result_Ecomor$value,
      ARI = ari_result_Ecomor,
      NMI = nmi_result_Ecomor,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  
  # Return the list of results
  return(results_list)
}

# Define your cluster results and times as a list
cluster_result_Ecomor <- list(
  KMeans = kmeans_clusters_Ecomor,
  CLARA = clara_clusters_Ecomor,
  PAM = pam_clusters_Ecomor,
  Hierarchical = hclust_clusters_Ecomor,
  Mclust = mclust_clusters_Ecomor,
  DPMVN = dpMVN_clusters_Ecomor,
  True = true_lab_Ecomor
)

times_list_Pu <- list(
  KMeans = kmeans_time_Ecomor,
  CLARA = clara_time_Ecomor,
  PAM = pam_time_Ecomor,
  Hierarchical = hclust_time_Ecomor,
  Mclust = mclust_time_Ecomor,
  DPMVN = DPMVN_time_Ecomor
)

# Call the function
clustering_metricsEcomor_2 <- calculate_clustering_metricsEcomor_2(dataEcomor_2 = Ecomor_scaled_data, true_clusters_Ecomor = true_lab_Ecomor, estimated_clusters_Ecomor_list = cluster_result_Ecomor, times_list = times_list_Pu)

# Print the results for each method
for (method_name in names(clustering_metricsEcomor_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsEcomor_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsEcomor_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsEcomor_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsEcomor_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsEcomor_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#################################################################################################################
#                                                 Save Results in Excel
#################################################################################################################

row_Ecomor_K <- data.frame(clustering_metricsEcomor_2$True$Kappa, dp_Gmvlg_metricsEcomor_1$dp_Gmvlg$Kappa, clustering_metricsEcomor_2$KMeans$Kappa, clustering_metricsEcomor_2$CLARA$Kappa, clustering_metricsEcomor_2$PAM$Kappa, clustering_metricsEcomor_2$Hierarchical$Kappa, clustering_metricsEcomor_2$Mclust$Kappa, clustering_metricsEcomor_2$DPMVN$Kappa)
row_Ecomor_ARI <- data.frame(clustering_metricsEcomor_2$True$ARI, dp_Gmvlg_metricsEcomor_1$dp_Gmvlg$ARI, clustering_metricsEcomor_2$KMeans$ARI, clustering_metricsEcomor_2$CLARA$ARI, clustering_metricsEcomor_2$PAM$ARI, clustering_metricsEcomor_2$Hierarchical$ARI, clustering_metricsEcomor_2$Mclust$ARI, clustering_metricsEcomor_2$DPMVN$ARI)
row_Ecomor_NMI <- data.frame(clustering_metricsEcomor_2$True$NMI, dp_Gmvlg_metricsEcomor_1$dp_Gmvlg$NMI, clustering_metricsEcomor_2$KMeans$NMI, clustering_metricsEcomor_2$CLARA$NMI, clustering_metricsEcomor_2$PAM$NMI, clustering_metricsEcomor_2$Hierarchical$NMI, clustering_metricsEcomor_2$Mclust$NMI, clustering_metricsEcomor_2$DPMVN$NMI)
row_Ecomor_Cpu <- data.frame(dp_Gmvlg_metricsEcomor_1$dp_Gmvlg$CPU_RUN_TIME, clustering_metricsEcomor_2$KMeans$CPU_RUN_TIME, clustering_metricsEcomor_2$CLARA$CPU_RUN_TIME, clustering_metricsEcomor_2$PAM$CPU_RUN_TIME, clustering_metricsEcomor_2$Hierarchical$CPU_RUN_TIME, clustering_metricsEcomor_2$Mclust$CPU_RUN_TIME, clustering_metricsEcomor_2$DPMVN$CPU_RUN_TIME)

tableEcomor_1 <- clustering_metricsEcomor_2$True$Table
tableEcomor_2 <- dp_Gmvlg_metricsEcomor_1$dp_Gmvlg$Table
tableEcomor_3 <- clustering_metricsEcomor_2$KMeans$Table
tableEcomor_4 <- clustering_metricsEcomor_2$CLARA$Table
tableEcomor_5 <- clustering_metricsEcomor_2$PAM$Table
tableEcomor_6 <- clustering_metricsEcomor_2$Hierarchical$Table
tableEcomor_7 <- clustering_metricsEcomor_2$Mclust$Table
tableEcomor_8 <- clustering_metricsEcomor_2$DPMVN$Table

colnames(row_Ecomor_K) <- NULL
colnames(row_Ecomor_ARI) <- NULL
colnames(row_Ecomor_NMI) <- NULL
colnames(row_Ecomor_Cpu) <- NULL

row_Ecomor_K <- as.matrix(row_Ecomor_K)
row_Ecomor_ARI <- as.matrix(row_Ecomor_ARI)
row_Ecomor_NMI <- as.matrix(row_Ecomor_NMI)
row_Ecomor_Cpu <- as.matrix(row_Ecomor_Cpu)

kappa_table_Ecomor <- rbind(row_Ecomor_K)
ARI_table_Ecomor <- rbind(row_Ecomor_ARI)
NMI_table_Ecomor <- rbind(row_Ecomor_NMI)
cpu_runtime_table_Ecomor <- rbind(row_Ecomor_Cpu)
clustering_table_Ecomor <- rbind(tableEcomor_1,tableEcomor_2,tableEcomor_3,tableEcomor_4,tableEcomor_5,tableEcomor_6,tableEcomor_7)

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each metric as a separate sheet in the workbook
addWorksheet(wb, "Kappa")
writeData(wb, sheet = "Kappa", kappa_table_Ecomor)

addWorksheet(wb, "ARI")
writeData(wb, sheet = "ARI", ARI_table_Ecomor)

addWorksheet(wb, "NMI")
writeData(wb, sheet = "NMI", NMI_table_Ecomor)

addWorksheet(wb, "CPU Run Time")
writeData(wb, sheet = "CPU Run Time", cpu_runtime_table_Ecomor)

addWorksheet(wb, "clustering table")
writeData(wb, sheet = "clustering table", clustering_table_Ecomor)

# Save the workbook to an Excel file
saveWorkbook(wb, "C:/Users/S4031636/OneDrive - RMIT University/PhD Topic Clustering/PAPER 1 DPG-MVLG via JAGS/Numerical Results/REAL DATASETS/12_ECOMOR/Cluster_metricsECOMOR.xlsx", overwrite = TRUE)
