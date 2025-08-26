#################################################################
#                       Wine DATASET                           ##
#################################################################

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.

#--------------------------#
# Load necessary libraries
#--------------------------#
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
library(genie)
library(gtools) 
library(combinat)
library(clue)
library(irr)

#---------------------------------------------------#
#             Load the Wine Dataset
#---------------------------------------------------#
data("wine",package = "rattle")
wine_data <- as.matrix(wine[,-1]) 
true_lab_wine <- as.numeric(wine$Type)    

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                          DPGMVLG
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#--------------------------------------#
#      Scaling the data
#--------------------------------------#
wine_scaled_data <- scale(wine_data)

#--------------------------------------#
#      Plotting true labels
#--------------------------------------#
fviz_cluster(list(data = wine_scaled_data, cluster = true_lab_wine), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#--------------------------------------#
#  Number of rows and columns
#--------------------------------------#
P <- nrow(wine_scaled_data)
D <- ncol(wine_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 3  #maximum no. of clusters

#--------------------------------------#
#      Setting up data for JAGS
#--------------------------------------#
data_list <- list(
  x = wine_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(wine_scaled_data))^(1 / (D - 1)),
  G = G 
)

#-------------------------------------#
#           JAGS MODEL
#-------------------------------------#
model_string <- "
Data {
  C <- 1000000000
  for (i in 1:P) {
     zeros[i] <- 0
  }
  v_g <- 3.991474683 
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
       mu_g[j, g] ~ dgamma(0.756493454,14.25714222)  
       lambda_g[j, g] ~ dgamma(5.448501062,12.85948424)
    }
    alpha[g] ~ dbeta(7.4189,4.1379) 
   }
   pi_g[1:G] ~ ddirch(alpha/sum(alpha))
}
"
#-----------------------------------------------#
#            Writing the model to a file
#-----------------------------------------------#
writeLines(model_string, con = "TEMPmodel.txt")

#-----------------------------------------------#
#             Parameters to monitor
#-----------------------------------------------#
params <- c("delta_g","mu_g", "lambda_g", "cluster", "z")

#-----------------------------------------------#
#               Reproducibility
#-----------------------------------------------#
inits <- list(list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 22021),
              list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 32019))

#-----------------------------------------------#
# Tracking the Run time and Running the JAGS model
#-----------------------------------------------#
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
run_time_wine <- toc()
codaSamples_wine = as.mcmc.list( runJagsOut )
summaryChains_wine <- summary(codaSamples_wine)

#----------------------------------------------#
#               Diagnostics
#----------------------------------------------#
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_wine , parName="mu_g[2,2]" )
graphics.off()

#---------------------------------------------#
#             Summary Chains
#---------------------------------------------#
z_mode_wine <- apply(matrix(summaryChains_wine$statistics[(1+13*G+13*G+P+1):nrow(summaryChains_wine$statistics),1], P, G),1, which.max)
z_mode_wine

#--------------------------------------------#
#         Plotting Predicted labels
#--------------------------------------------#
fviz_cluster(list(data = wine_scaled_data, cluster = z_mode_wine), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

plot(wine_scaled_data, col= z_mode_wine, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_wine), col=unique(z_mode_wine), pch=16, title="Cluster")
table( true_lab_wine , z_mode_wine)

#--------------------------------------------#
#               Relabeling 
#--------------------------------------------#
new <- c(1,3)
old <- c(3,1)
z_mode_wine[z_mode_wine %in% old] <- new[match(z_mode_wine,old,nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
z_mode_wine[z_mode_wine %in% old] <- new[match(z_mode_wine,old,nomatch = 0)]

plot(wine_scaled_data, col= z_mode_wine, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_wine), col=unique(z_mode_wine), pch=16, title="Cluster")

table( true_lab_wine , true_lab_wine)
table( true_lab_wine , z_mode_wine)
kappa2(data.frame(rater1 = true_lab_wine, rater2 = z_mode_wine))

#--------------------------------------------#
#           DPGMVLG Results
#--------------------------------------------#
calculate_dp_Gmvlg_clustering_metrics <- function(data_1, true_clusters, z_mode) {

  results_list <- list()

  table_result <- table(true_clusters, z_mode)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = z_mode))
  ari_result <- adjustedRandIndex(true_clusters, z_mode)
  nmi_result <- NMI(true_clusters, z_mode)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result,
    CPU_RUN_TIME = run_time_wine$callback_msg  # Store the runtime
  )
  return(results_list)
}

#--------------------#
# Call the function
#--------------------#
dp_Gmvlg_metricswine_1 <- calculate_dp_Gmvlg_clustering_metrics(data_1 = wine_scaled_data, true_clusters = true_lab_wine, z_mode = z_mode_wine)

#--------------------#
#       Results 
#--------------------#
for (dp_Gmvlg in names(dp_Gmvlg_metricswine_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricswine_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                   All Classical methods and One Bayesian method (DPMVN) --- Except Genie
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#                            K-means clustering
#-------------------------------------------------------------------------#
set.seed(123)
tic("K-means Runtime")
kmeans_result_wine <- eclust(wine_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_wine <- toc()
kmeans_clusters_wine <- kmeans_result_wine$cluster
table(true_lab_wine, kmeans_clusters_wine)

new <- 1:2
old <- c(2,1)
kmeans_clusters_wine[kmeans_clusters_wine %in% old] <- new[match(kmeans_clusters_wine, old, nomatch = 0)]
new <- c(2,3)
old <- c(3,2)
kmeans_clusters_wine[kmeans_clusters_wine %in% old] <- new[match(kmeans_clusters_wine, old, nomatch = 0)]
table(true_lab_wine, kmeans_clusters_wine)

#-------------------------------------------------------------------------#
#                            CLARA clustering
#-------------------------------------------------------------------------#
tic("CLARA Runtime")
clara_result_wine <- eclust(wine_scaled_data, "clara", G, graph = FALSE)
clara_time_wine <- toc()
clara_clusters_wine <- clara_result_wine$cluster
table(true_lab_wine, clara_clusters_wine)

#-------------------------------------------------------------------------#
#                            PAM clustering
#-------------------------------------------------------------------------#
tic("PAM Runtime")
pam_result_wine <- eclust(wine_scaled_data, "pam", G, graph = FALSE)
pam_time_wine <- toc()
pam_clusters_wine <- pam_result_wine$cluster
table(true_lab_wine, pam_clusters_wine)

#-------------------------------------------------------------------------#
#                       Hierarchical clustering
#-------------------------------------------------------------------------#
tic("Hierarchical Runtime")
hclust_result_wine <- hclust(dist(wine_scaled_data), method = "ward.D2")
hclust_time_wine <- toc()
hclust_clusters_wine <- cutree(hclust_result_wine, k = G)
table(true_lab_wine, hclust_clusters_wine)

#------------------------------------------------------------------------#
#                       Model-based clustering
#------------------------------------------------------------------------#
tic("Model-based Runtime")
mclust_result_wine <- Mclust(wine_scaled_data, G = G)
mclust_time_wine <- toc()
summary(mclust_result_wine)
mclust_clusters_wine <- mclust_result_wine$classification
table(true_lab_wine, mclust_clusters_wine)

#-------------------------------------------------------------------------------------#
# Bayesian Clustering using Dirichlet Process Multivariate Normal Distribution (DPMVN)
#-------------------------------------------------------------------------------------#
# Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(wine_scaled_data, alphaPriors = c(50,30))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_wine <- toc()

# Extract clusters 
dpMVN_clusters_wine <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_wine)
table(true_lab_wine, dpMVN_clusters_wine)

#--------------------------------------------------------------------------------------#
#                 All Classical Except Genie and one Bayesian (DPMVN) Results Together
#--------------------------------------------------------------------------------------#
calculate_clustering_metricswine_2 <- function(datawine_2, true_clusters_wine, estimated_clusters_wine_list, times_list) {

  results_list <- list()

  for (method_name in names(estimated_clusters_wine_list)) {
    clusters <- estimated_clusters_wine_list[[method_name]]

    table_result_wine <- table(true_clusters_wine, clusters)
    kappa_result_wine <- kappa2(data.frame(rater1 = true_clusters_wine, rater2 = clusters))
    ari_result_wine <- adjustedRandIndex(true_clusters_wine, clusters)
    nmi_result_wine <- NMI(true_clusters_wine, clusters)
    
    # Store results
    results_list[[method_name]] <- list(
      Table = table_result_wine,
      Kappa = kappa_result_wine$value,
      ARI = ari_result_wine,
      NMI = nmi_result_wine,
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  # Store the runtime
    )
  }
  return(results_list)
}

cluster_result_wine <- list(
  KMeans = kmeans_clusters_wine,
  CLARA = clara_clusters_wine,
  PAM = pam_clusters_wine,
  Hierarchical = hclust_clusters_wine,
  Mclust = mclust_clusters_wine,
  DPMVN = dpMVN_clusters_wine,
  True = true_lab_wine
)

times_list_wine <- list(
  KMeans = kmeans_time_wine,
  CLARA = clara_time_wine,
  PAM = pam_time_wine,
  Hierarchical = hclust_time_wine,
  Mclust = mclust_time_wine,
  DPMVN = DPMVN_time_wine
)

#-------------------#
# Call the function
#-------------------#
clustering_metricswine_2 <- calculate_clustering_metricswine_2(datawine_2 = wine_scaled_data, true_clusters_wine = true_lab_wine, estimated_clusters_wine_list = cluster_result_wine, times_list = times_list_wine)

#--------------------#
#       Results 
#--------------------#
for (method_name in names(clustering_metricswine_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricswine_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricswine_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricswine_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricswine_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricswine_2[[method_name]]$CPU_RUN_TIME, "\n")
}

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                          GENIE METHOD
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

library(genieclust)

#--------------------------------------------#
#          label mapping function
#--------------------------------------------#
relabel_clusters_optimal <- function(pred_clusters, true_labels) {
  
  conf_mat <- table(pred_clusters, true_labels)
  
  cost_mat <- max(conf_mat) - conf_mat
  
  assignment <- solve_LSAP(cost_mat)
  
  cluster_levels <- as.integer(rownames(conf_mat))
  mapping <- setNames(as.integer(colnames(conf_mat)[assignment]), cluster_levels)
  
  remapped <- mapping[as.character(pred_clusters)]
  
  return(list(
    remapped_labels = remapped,
    mapping = mapping,
    confusion = table(remapped, true_labels)
  ))
}

#-----------------------------------------------#
#             Gini values to test
#-----------------------------------------------#
gini_values <- seq(0.001, 1, 0.1)
kappa_scores <- numeric(length(gini_values))
names(kappa_scores) <- gini_values

k <- length(unique(true_lab_wine))  

for (i in seq_along(gini_values)) {
  gini <- gini_values[i]
  cat("\nTrying gini threshold:", gini, "\n")
  
  tic("Genie Runtime")
  res <- try(gclust(wine_scaled_data, gini_threshold = gini), silent = TRUE)
  run_time_wine <- toc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    cat(" Error or invalid clustering result\n")
    kappa_scores[i] <- NA
  } else {
    cluster_labels <- cutree(res, k = k)

    relabeled <- relabel_clusters_optimal(cluster_labels, true_lab_wine)
    aligned_labels <- relabeled$remapped_labels
    
    cat("  Confusion matrix (after relabeling):\n")
    print(relabeled$confusion)

    kappa_result <- kappa2(data.frame(rater1 = true_lab_wine, rater2 = aligned_labels))
    kappa_scores[i] <- kappa_result$value
    cat(" Kappa:", kappa_result$value, "\n")
  }
}

# final Kappa scores
cat("\nFinal Kappa scores by Gini threshold:\n")
print(kappa_scores)

#-------------------------------------------------------#
#           Running Genie with best Gini Coefficient
#-------------------------------------------------------#
tic("Genie Runtime")
res <- gclust(wine_scaled_data, gini_threshold = 0.001)
run_time_wine_Genie <- toc()

k <- length(unique(true_lab_wine))

aligned_labels <- cutree(res, k = k)

table( true_lab_wine , aligned_labels)
kappa2(data.frame(rater1 = true_lab_wine, rater2 = aligned_labels))

#---------------------------------------------------------------------------------#
#                            GENIE RESULTS
#---------------------------------------------------------------------------------#
calculate_Genie_clustering_metrics <- function(data_1, true_clusters, cluster_labels) {

  results_list <- list()

  table_result <- table(true_clusters, cluster_labels)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = cluster_labels))
  ari_result <- adjustedRandIndex(true_clusters, cluster_labels)
  nmi_result <- NMI(true_clusters, cluster_labels)
  
  # Store results
  results_list[["Genie"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result,
    CPU_RUN_TIME = run_time_wine_Genie$callback_msg  # Store the runtime
  )
  return(results_list)
}

#------------------#
# Call the function
#------------------#
Genie_metricswine_1 <- calculate_Genie_clustering_metrics(data_1 = wine_scaled_data, true_clusters = true_lab_wine, cluster_labels = aligned_labels)

#------------------#
#     Results 
#------------------#
for (Genie in names(Genie_metricswine_1)) {
  cat("\n", Genie, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(Genie_metricswine_1[[Genie]]$Table)
  cat("Kappa:", Genie_metricswine_1[[Genie]]$Kappa, "\n")
  cat("Adjusted Rand Index:", Genie_metricswine_1[[Genie]]$ARI, "\n")
  cat("NMI:", Genie_metricswine_1[[Genie]]$NMI, "\n")
  cat(" CPU RUN TIME:", Genie_metricswine_1[[Genie]]$CPU_RUN_TIME, "\n")
}

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                          Descriptive analysis
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------#
#           1. Summary statistics
#-------------------------------------------------------#
summary(wine_data)

#-------------------------------------------------------#
#           2. Check for missing values
#-------------------------------------------------------#
sum(is.na(wine_data))

#-------------------------------------------------------#
#           3. Distribution of each feature
#-------------------------------------------------------#
par(mfrow=c(4, 4)) 

for (i in 1:ncol(wine_data)) {
  hist(wine_data[, i], main=colnames(wine_data)[i], xlab=colnames(wine_data)[i], col="lightblue", border="black")
}

#-------------------------------------------------------#
#           4. Pairwise scatter plots
#-------------------------------------------------------#
pairs(wine_data, main="Wine Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3", "blue")[true_lab_wine])

#-------------------------------------------------------#
#           5. Boxplots for each feature 
#-------------------------------------------------------#
# Normalized data to the original data frame
wine_normalized <- as.data.frame(wine_data)
wine_normalized$Type <- as.factor(wine$Type)

par(mfrow=c(4, 4)) 

for (i in 1:ncol(wine_data)) {
  boxplot(wine_data[, i] ~ wine_normalized$Type, main=colnames(wine_data)[i], xlab="Type", ylab=colnames(wine_data)[i], col="lightblue")
}

#-------------------------------------------------------#
#           6. Correlation matrix
#-------------------------------------------------------#
cor_matrix <- cor(wine_data)
print(cor_matrix)

#-------------------------------------------------------#
#           7. Visualizing correlations with heatmap
#-------------------------------------------------------#
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

#-------------------------------------------------------#
#          Normalize the Dataset
#-------------------------------------------------------#
wine_scaled_data <- scale(wine_data)

#--------------------------------------------------------------#
#           8. Skewness and kurtosis for each feature
#--------------------------------------------------------------#
skewness_values <- apply(wine_scaled_data, 2, skewness)
kurtosis_values <- apply(wine_scaled_data, 2, kurtosis)

# Print skewness and kurtosis values
print("Skewness values:")
print(skewness_values)

print("Kurtosis values:")
print(kurtosis_values)

# Plot skewness values
barplot(skewness_values, main="Skewness of Each Feature", xlab="Features", ylab="Skewness", col="lightblue", las=2)

# Plot kurtosis values
barplot(kurtosis_values, main="Kurtosis of Each Feature", xlab="Features", ylab="Kurtosis", col="lightblue", las=2)

#----------------------------------------------------------#
#             Combined skewness and kurtosis
#----------------------------------------------------------#
combined_data <- as.vector(as.matrix(wine_scaled_data))
total_skewness <- skewness(combined_data)
total_skewness

total_kurtosis <- kurtosis(combined_data)
total_kurtosis

#---------------------------------------------------------#
#                        Outliers
#---------------------------------------------------------#
total_outliers1 <- 0
outliers1 <- list()
outlier_rows1 <- c()
for (i in 1:ncol(wine_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(wine_scaled_data[,i])$out
  print(paste("Feature:", colnames(wine_scaled_data)[i]))
  print(paste("Outliers:", length(outliers1[[i]])))
  total_outliers1 <- total_outliers1 + length(outliers1[[i]])
  outlier_rows1 <- unique(c(outlier_rows1, outliers1[[i]])) 
}

# Total number of outliers after the loop
print(paste("Total outliers:", total_outliers1))

# Number of rows with outliers
print(paste("Number of rows with outliers:", length(outlier_rows1)))

# Detect_outlier function
detect_outlier <- function(x) {
  
  # the first quantile
  Quantile1 <- quantile(x, probs = 0.25)
  
  # the third quantile
  Quantile3 <- quantile(x, probs = 0.75)
  
  # interquartile range (IQR)
  IQR <- Quantile3 - Quantile1

  outliers <- x > Quantile3 + (IQR * 1.5) | x < Quantile1 - (IQR * 1.5)
  return(outliers)
}

# Detect_outlier function to each column
outliers_matrix <- apply(wine_scaled_data, 2, detect_outlier)
outliers_matrix

# Number of rows that have at least one TRUE (outlier)
rows_with_outliers <- rowSums(outliers_matrix) > 0
num_rows_with_outliers <- sum(rows_with_outliers)

# Result
print(paste("Number of rows with at least one outlier:", num_rows_with_outliers))
