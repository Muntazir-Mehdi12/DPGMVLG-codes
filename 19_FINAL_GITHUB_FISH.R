##################################################
#                FISH CATCH DATASET             ##
##################################################

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
#             Load the Fish Dataset
#---------------------------------------------------#
data("fish", package = "rrcov")
Fish_data <- as.matrix(fish[,c(3,5,6)])
true_lab_Fish=fish[,7]

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
Fish_scaled_data <- scale(Fish_data)

#--------------------------------------#
#      Plotting true labels
#--------------------------------------#
fviz_cluster(list(data = Fish_scaled_data, cluster = true_lab_Fish), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#--------------------------------------#
#  Number of rows and columns
#--------------------------------------#
P <- nrow(Fish_scaled_data)
D <- ncol(Fish_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 7  #maximum no. of clusters

#--------------------------------------#
#      Setting up data for JAGS
#--------------------------------------#
data_list <- list(
  x = Fish_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Fish_scaled_data))^(1 / (D - 1)),
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
  v_g <- 2.918224499999 #2.918224499999
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
       mu_g[j, g] ~ dgamma(4.764564856,9.96240223)     #dgamma(4.764564856,9.96240223)
       lambda_g[j, g] ~ dgamma(8.21516072,8.2906433)  #dgamma(8.21516072,8.2906433)
    }
     alpha[g] ~ dgamma(1,1)  #dgamma(1,1)
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
run_time_Fish <- toc()
codaSamples_Fish = as.mcmc.list( runJagsOut )
summaryChains_Fish <- summary(codaSamples_Fish)

#----------------------------------------------#
#               Diagnostics
#----------------------------------------------#
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Fish , parName="mu_g[2,2]" )
graphics.off()

#---------------------------------------------#
#             Summary Chains
#---------------------------------------------#
z_mode_Fish <- apply(matrix(summaryChains_Fish$statistics[(1+3*G+3*G+P+1):nrow(summaryChains_Fish$statistics),1], P, G),1, which.max)
z_mode_Fish

#--------------------------------------------#
#         Plotting Predicted labels
#--------------------------------------------#
fviz_cluster(list(data = Fish_scaled_data, cluster = z_mode_Fish), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")
fviz_cluster(list(data = Fish_scaled_data, cluster = true_lab_Fish), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

plot(Fish_scaled_data, col= z_mode_Fish, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Fish), col=unique(z_mode_Fish), pch=16, title="Cluster")
table( true_lab_Fish , z_mode_Fish)
table( true_lab_Fish , true_lab_Fish)

#--------------------------------------------#
#           Relabeling using a function
#--------------------------------------------#
library(clue)

# label mapping function
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

# relabeling function
relabeled <- relabel_clusters_optimal(z_mode_Fish, true_lab_Fish)
z_mode_Fish <- relabeled$remapped_labels

cat("  Confusion matrix (after relabeling):\n")
print(relabeled$confusion)

# Cohen's Kappa
kappa_result <- kappa2(data.frame(rater1 = true_lab_Fish, rater2 = z_mode_Fish))
kappa_scores <- kappa_result$value
cat("  ✅ Kappa:", kappa_result$value, "\n")

table( true_lab_Fish , z_mode_Fish)
kappa2(data.frame(rater1 = true_lab_Fish, rater2 = z_mode_Fish))

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
    CPU_RUN_TIME = run_time_Fish$callback_msg  
  )
  
  return(results_list)
}

#--------------------#
# Call the function
#--------------------#
dp_Gmvlg_metricsFish_1 <- calculate_dp_Gmvlg_clustering_metrics(data_1 = Fish_scaled_data, true_clusters = true_lab_Fish, z_mode = z_mode_Fish)

#-------------------#
#    Results
#-------------------#
for (dp_Gmvlg in names(dp_Gmvlg_metricsFish_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsFish_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
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
kmeans_result_Fish <- eclust(Fish_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Fish <- toc()
kmeans_clusters_Fish <- kmeans_result_Fish$cluster
table(true_lab_Fish, kmeans_clusters_Fish)

# relabeling function 
relabeled_kmeans <- relabel_clusters_optimal(kmeans_clusters_Fish, true_lab_Fish)
kmeans_clusters_Fish <- relabeled_kmeans$remapped_labels

#-------------------------------------------------------------------------#
#                            CLARA clustering
#-------------------------------------------------------------------------#
tic("CLARA Runtime")
clara_result_Fish <- eclust(Fish_scaled_data, "clara", G, graph = FALSE)
clara_time_Fish <- toc()
clara_clusters_Fish <- clara_result_Fish$cluster
table(true_lab_Fish, clara_clusters_Fish)

# relabeling function 
relabeled_clara <- relabel_clusters_optimal(clara_clusters_Fish, true_lab_Fish)
clara_clusters_Fish <- relabeled_clara$remapped_labels

#-------------------------------------------------------------------------#
#                            PAM clustering
#-------------------------------------------------------------------------#
tic("PAM Runtime")
pam_result_Fish <- eclust(Fish_scaled_data, "pam", G, graph = FALSE)
pam_time_Fish <- toc()
pam_clusters_Fish <- pam_result_Fish$cluster
table(true_lab_Fish, pam_clusters_Fish)

# relabeling function 
relabeled_pam <- relabel_clusters_optimal(pam_clusters_Fish, true_lab_Fish)
pam_clusters_Fish <- relabeled_pam$remapped_labels

#-------------------------------------------------------------------------#
#                       Hierarchical clustering
#-------------------------------------------------------------------------#
tic("Hierarchical Runtime")
hclust_result_Fish <- hclust(dist(Fish_scaled_data), method = "ward.D2")
hclust_time_Fish <- toc()
hclust_clusters_Fish <- cutree(hclust_result_Fish, k = G)
table(true_lab_Fish, hclust_clusters_Fish)

# relabeling function
relabeled_hclust <- relabel_clusters_optimal(hclust_clusters_Fish, true_lab_Fish)
hclust_clusters_Fish <- relabeled_hclust$remapped_labels

#------------------------------------------------------------------------#
#                       Model-based clustering
#------------------------------------------------------------------------#
tic("Model-based Runtime")
mclust_result_Fish <- Mclust(Fish_scaled_data, G = G)
mclust_time_Fish <- toc()
summary(mclust_result_Fish)
mclust_clusters_Fish <- mclust_result_Fish$classification
table(true_lab_Fish, mclust_clusters_Fish)

# relabeling function 
relabeled_mclust <- relabel_clusters_optimal(mclust_clusters_Fish, true_lab_Fish)
mclust_clusters_Fish <- relabeled_mclust$remapped_labels

#-------------------------------------------------------------------------------------#
# Bayesian Clustering using Dirichlet Process Multivariate Normal Distribution (DPMVN)
#-------------------------------------------------------------------------------------#
# Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Fish_scaled_data, alphaPriors = c(6,19))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Fish <- toc()

# Extract clusters 
dpMVN_clusters_Fish <- as.numeric(dp$clusterLabels)
table(true_lab_Fish, dpMVN_clusters_Fish)

# relabeling function 
relabeled_dpMVN <- relabel_clusters_optimal(dpMVN_clusters_Fish, true_lab_Fish)
dpMVN_clusters_Fish <- relabeled_dpMVN$remapped_labels

#--------------------------------------------------------------------------------------#
#                 All Classical Except Genie and one Bayesian (DPMVN) Results Together
#--------------------------------------------------------------------------------------#
calculate_clustering_metricsFish_2 <- function(dataFish_2, true_clusters_Fish, estimated_clusters_Fish_list, times_list) {

  results_list <- list()

  for (method_name in names(estimated_clusters_Fish_list)) {
    clusters <- estimated_clusters_Fish_list[[method_name]]

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
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg  
    )
  }
  return(results_list)
}

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

#-------------------#
# Call the function
#-------------------#
clustering_metricsFish_2 <- calculate_clustering_metricsFish_2(dataFish_2 = Fish_scaled_data, true_clusters_Fish = true_lab_Fish, estimated_clusters_Fish_list = cluster_result_Fish, times_list = times_list_Fish)

#-------------------#
#    Results
#-------------------#
for (method_name in names(clustering_metricsFish_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsFish_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsFish_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsFish_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsFish_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsFish_2[[method_name]]$CPU_RUN_TIME, "\n")
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

k <- length(unique(true_lab_Fish)) 

for (i in seq_along(gini_values)) {
  gini <- gini_values[i]
  cat("\nTrying gini threshold:", gini, "\n")
  
  tic("Genie Runtime")
  res <- try(gclust(Fish_scaled_data, gini_threshold = gini), silent = TRUE)
  run_time_Fish <- toc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    cat(" Error or invalid clustering result\n")
    kappa_scores[i] <- NA
  } else {
    cluster_labels <- cutree(res, k = k)

    relabeled <- relabel_clusters_optimal(cluster_labels, true_lab_Fish)
    aligned_labels <- relabeled$remapped_labels
    
    cat("  Confusion matrix (after relabeling):\n")
    print(relabeled$confusion)

    kappa_result <- kappa2(data.frame(rater1 = true_lab_Fish, rater2 = aligned_labels))
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
res <- gclust(Fish_scaled_data, gini_threshold = 0.901)
run_time_Fish_Genie <- toc()

k <- length(unique(true_lab_Fish))

aligned_labels <- cutree(res, k = k)

table( true_lab_Fish , aligned_labels)
kappa2(data.frame(rater1 = true_lab_Fish, rater2 = aligned_labels))

#-------------------------------------------------------#
#            Relabeling
#-------------------------------------------------------#
relabeled_Genie <- relabel_clusters_optimal(aligned_labels, true_lab_Fish)
aligned_labels <- relabeled_Genie$remapped_labels

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
    CPU_RUN_TIME = run_time_Fish_Genie$callback_msg  
  )
  
  return(results_list)
}

#------------------#
# Call the function
#------------------#
Genie_metricsFish_1 <- calculate_Genie_clustering_metrics(data_1 = Fish_scaled_data, true_clusters = true_lab_Fish, cluster_labels = aligned_labels)

#------------------#
#     Results 
#------------------#
for (Genie in names(Genie_metricsFish_1)) {
  cat("\n", Genie, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(Genie_metricsFish_1[[Genie]]$Table)
  cat("Kappa:", Genie_metricsFish_1[[Genie]]$Kappa, "\n")
  cat("Adjusted Rand Index:", Genie_metricsFish_1[[Genie]]$ARI, "\n")
  cat("NMI:", Genie_metricsFish_1[[Genie]]$NMI, "\n")
  cat(" CPU RUN TIME:", Genie_metricsFish_1[[Genie]]$CPU_RUN_TIME, "\n")
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
summary(Fish_data)

#-------------------------------------------------------#
#           2. Check for missing values
#-------------------------------------------------------#
sum(is.na(Fish_data))

#-------------------------------------------------------#
#           3. Distribution of each feature
#-------------------------------------------------------#
par(mfrow=c(1, 3))  

for (i in 1:ncol(Fish_data)) {
  hist(Fish_data[, i], main=colnames(Fish_data)[i], xlab=colnames(Fish_data)[i], col="lightblue", border="black")
}

#-------------------------------------------------------#
#           4. Pairwise scatter plots
#-------------------------------------------------------#
pairs(Fish_data, main="Fish Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Fish])

#-------------------------------------------------------#
#           5. Boxplots for each feature 
#-------------------------------------------------------#
# Normalized data to the original data frame
fish_normalized <- as.data.frame(Fish_data)
fish_normalized$Species <- fish$Species

par(mfrow=c(1, 3))  

for (i in 1:ncol(Fish_data)) {
  boxplot(Fish_data[, i] ~ fish$Species, main=colnames(Fish_data)[i], xlab="Species", ylab=colnames(Fish_data)[i], col="lightblue")
}

#-------------------------------------------------------#
#           6. Correlation matrix
#-------------------------------------------------------#
cor_matrix <- cor(Fish_data)
print(cor_matrix)

#-------------------------------------------------------#
#           7. Visualizing correlations with heatmap
#-------------------------------------------------------#
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

#-------------------------------------------------------#
#              Normalize the Dataset
#-------------------------------------------------------#
Fish_scaled_data <- scale(Fish_data)

#--------------------------------------------------------------#
#           8. Skewness and kurtosis for each feature
#--------------------------------------------------------------#
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

#----------------------------------------------------------#
#             Combined skewness and kurtosis
#----------------------------------------------------------#
combined_data <- as.vector(as.matrix(Fish_scaled_data))
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
for (i in 1:ncol(Fish_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Fish_scaled_data[,i])$out
  print(paste("Feature:", colnames(Fish_scaled_data)[i]))
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
outliers_matrix <- apply(Fish_scaled_data, 2, detect_outlier)
outliers_matrix

# Number of rows that have at least one TRUE (outlier)
rows_with_outliers <- rowSums(outliers_matrix) > 0
num_rows_with_outliers <- sum(rows_with_outliers)

# Result
print(paste("Number of rows with at least one outlier:", num_rows_with_outliers))

