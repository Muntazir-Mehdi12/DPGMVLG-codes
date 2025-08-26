#################################################################
#                        CATS DATASET                          ##
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
#             Load the Cats Dataset
#---------------------------------------------------#
data("cats")
help(cats)
Cat_data <- as.matrix(cats[,-1]) # Convert all columns to numeric
true_lab_Cat <- as.numeric(cats$Sex)    # Use the 'class' column as the true labels

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
Cat_scaled_data <- scale(Cat_data)

#--------------------------------------#
#      Plotting true labels
#--------------------------------------#
fviz_cluster(list(data = Cat_scaled_data, cluster = true_lab_Cat), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#--------------------------------------#
#    Number of rows and columns
#--------------------------------------#
P <- nrow(Cat_scaled_data)
D <- ncol(Cat_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 2  #maximum no. of clusters

#--------------------------------------#
#      Setting up data for JAGS
#--------------------------------------#
data_list <- list(
  x = Cat_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Cat_scaled_data))^(1 / (D - 1)),
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
  v_g <- 4.39942111478560
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
       mu_g[j, g] ~ dgamma(10.9929535878589,9.76244791357312)         
       lambda_g[j, g] ~ dgamma(11.0907771108439,9.82007480405737)     
    }
    alpha[g] ~ dgamma(0.99,9.46) 
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
run_time_Cat <- toc()
codaSamples_Cat = as.mcmc.list( runJagsOut )
summaryChains_Cat <- summary(codaSamples_Cat)

#----------------------------------------------#
#               Diagnostics
#----------------------------------------------#
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Cat , parName="mu_g[2,2]" )
graphics.off()

#---------------------------------------------#
#             Summary Chains
#---------------------------------------------#
matrix(summaryChains_Cat$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Cat$statistics),1], P, G)
z_mode_Cat <- apply(matrix(summaryChains_Cat$statistics[(1+2*G+2*G+P+1):nrow(summaryChains_Cat$statistics),1], P, G),1, which.max)
z_mode_Cat

#--------------------------------------------#
#         Plotting Predicted labels
#--------------------------------------------#
fviz_cluster(list(data = Cat_scaled_data, cluster = z_mode_Cat), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

plot(Cat_scaled_data, col= z_mode_Cat, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Cat), col=unique(z_mode_Cat), pch=16, title="Cluster")
table( true_lab_Cat , z_mode_Cat)

#--------------------------------------------#
#               Relabeling 
#--------------------------------------------#
new <- 1:2
old <- c(2,1)
z_mode_Cat[z_mode_Cat %in% old] <- new[match(z_mode_Cat,old,nomatch = 0)]

plot(Cat_scaled_data, col= z_mode_Cat, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Cat), col=unique(z_mode_Cat), pch=16, title="Cluster")

table( true_lab_Cat , true_lab_Cat)
table( true_lab_Cat , z_mode_Cat)
kappa2(data.frame(rater1 = true_lab_Cat, rater2 = z_mode_Cat))

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
    CPU_RUN_TIME = run_time_Cat$callback_msg  
  )
  
  return(results_list)
}

#--------------------#
# Call the function
#--------------------#
dp_Gmvlg_metricsCat_1 <- calculate_dp_Gmvlg_clustering_metrics(data_1 = Cat_scaled_data, true_clusters = true_lab_Cat, z_mode = z_mode_Cat)

#-------------------#
#    Results
#-------------------#
for (dp_Gmvlg in names(dp_Gmvlg_metricsCat_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsCat_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
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
kmeans_result_Cat <- eclust(Cat_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Cat <- toc()
kmeans_clusters_Cat <- kmeans_result_Cat$cluster
table(true_lab_Cat, kmeans_clusters_Cat)

#-------------------------------------------------------------------------#
#                            CLARA clustering
#-------------------------------------------------------------------------#
tic("CLARA Runtime")
clara_result_Cat <- eclust(Cat_scaled_data, "clara", G, graph = FALSE)
clara_time_Cat <- toc()
clara_clusters_Cat <- clara_result_Cat$cluster
table(true_lab_Cat, clara_clusters_Cat)

#-------------------------------------------------------------------------#
#                            PAM clustering
#-------------------------------------------------------------------------#
tic("PAM Runtime")
pam_result_Cat <- eclust(Cat_scaled_data, "pam", G, graph = FALSE)
pam_time_Cat <- toc()
pam_clusters_Cat <- pam_result_Cat$cluster
table(true_lab_Cat, pam_clusters_Cat)

#-------------------------------------------------------------------------#
#                       Hierarchical clustering
#-------------------------------------------------------------------------#
tic("Hierarchical Runtime")
hclust_result_Cat <- hclust(dist(Cat_scaled_data), method = "ward.D2")
hclust_time_Cat <- toc()
hclust_clusters_Cat <- cutree(hclust_result_Cat, k = G)
table(true_lab_Cat, hclust_clusters_Cat)

#------------------------------------------------------------------------#
#                       Model-based clustering
#------------------------------------------------------------------------#
tic("Model-based Runtime")
mclust_result_Cat <- Mclust(Cat_scaled_data, G = G)
mclust_time_Cat <- toc()
summary(mclust_result_Cat)
mclust_clusters_Cat <- mclust_result_Cat$classification
table(true_lab_Cat, mclust_clusters_Cat)

#-------------------------------------------------------------------------------------#
# Bayesian Clustering using Dirichlet Process Multivariate Normal Distribution (DPMVN)
#-------------------------------------------------------------------------------------#
# Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Cat_scaled_data)
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Cat <- toc()

# Extract clusters 
dpMVN_clusters_Cat <- as.numeric(dp$clusterLabels)
print(dpMVN_clusters_Cat)
table(true_lab_Cat, dpMVN_clusters_Cat)

#--------------------------------------------------------------------------------------#
#                 All Classical Except Genie and one Bayesian (DPMVN) Results Together
#--------------------------------------------------------------------------------------#
calculate_clustering_metricsCat_2 <- function(dataCat_2, true_clusters_Cat, estimated_clusters_Cat_list, times_list) {

  results_list <- list()

  for (method_name in names(estimated_clusters_Cat_list)) {
    clusters <- estimated_clusters_Cat_list[[method_name]]

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
      CPU_RUN_TIME = times_list[[method_name]]$callback_msg 
    )
  }

  return(results_list)
}

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

#-------------------#
# Call the function
#-------------------#
clustering_metricsCat_2 <- calculate_clustering_metricsCat_2(dataCat_2 = Cat_scaled_data, true_clusters_Cat = true_lab_Cat, estimated_clusters_Cat_list = cluster_result_Cat, times_list = times_list_Cat)

#-------------------#
#    Results
#-------------------#
for (method_name in names(clustering_metricsCat_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsCat_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsCat_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsCat_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsCat_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsCat_2[[method_name]]$CPU_RUN_TIME, "\n")
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

k <- length(unique(true_lab_Cat)) 

for (i in seq_along(gini_values)) {
  gini <- gini_values[i]
  cat("\nTrying gini threshold:", gini, "\n")
  
  tic("Genie Runtime")
  res <- try(gclust(Cat_scaled_data, gini_threshold = gini), silent = TRUE)
  run_time_Cat <- toc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    cat(" Error or invalid clustering result\n")
    kappa_scores[i] <- NA
  } else {
    cluster_labels <- cutree(res, k = k)

    relabeled <- relabel_clusters_optimal(cluster_labels, true_lab_Cat)
    aligned_labels <- relabeled$remapped_labels
    
    cat("  Confusion matrix (after relabeling):\n")
    print(relabeled$confusion)

    kappa_result <- kappa2(data.frame(rater1 = true_lab_Cat, rater2 = aligned_labels))
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
res <- gclust(Cat_scaled_data, gini_threshold = 0.001)
run_time_Cat_Genie <- toc()

k <- length(unique(true_lab_Cat))

aligned_labels <- cutree(res, k = k)

table( true_lab_Cat , aligned_labels)
kappa2(data.frame(rater1 = true_lab_Cat, rater2 = aligned_labels))

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
    CPU_RUN_TIME = run_time_Cat_Genie$callback_msg  # Store the runtime
  )
  
  return(results_list)
}

#------------------#
# Call the function
#------------------#
Genie_metricsCat_1 <- calculate_Genie_clustering_metrics(data_1 = Cat_scaled_data, true_clusters = true_lab_Cat, cluster_labels = aligned_labels)

#------------------#
#     Results 
#------------------#
for (Genie in names(Genie_metricsCat_1)) {
  cat("\n", Genie, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(Genie_metricsCat_1[[Genie]]$Table)
  cat("Kappa:", Genie_metricsCat_1[[Genie]]$Kappa, "\n")
  cat("Adjusted Rand Index:", Genie_metricsCat_1[[Genie]]$ARI, "\n")
  cat("NMI:", Genie_metricsCat_1[[Genie]]$NMI, "\n")
  cat(" CPU RUN TIME:", Genie_metricsCat_1[[Genie]]$CPU_RUN_TIME, "\n")
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
summary(Cat_data)

#-------------------------------------------------------#
#           2. Check for missing values
#-------------------------------------------------------#
sum(is.na(Cat_data))

#-------------------------------------------------------#
#           3. Distribution of each feature
#-------------------------------------------------------#
par(mfrow=c(1, 2))  

for (i in 1:ncol(Cat_data)) {
  hist(Cat_data[, i], main=colnames(Cat_data)[i], xlab=colnames(Cat_data)[i], col="lightblue", border="black")
}

#-------------------------------------------------------#
#           4. Pairwise scatter plots
#-------------------------------------------------------#
pairs(Cat_data, main="Cats Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Cat])

#-------------------------------------------------------#
#           5. Boxplots for each feature 
#-------------------------------------------------------#
# Normalized data to the original data frame
cats_normalized <- as.data.frame(Cat_data)
cats_normalized$Sex <- as.factor(cats$Sex)

par(mfrow=c(1, 2)) 

for (i in 1:ncol(Cat_data)) {
  boxplot(Cat_data[, i] ~ cats_normalized$Sex, main=colnames(Cat_data)[i], xlab="Sex", ylab=colnames(Cat_data)[i], col="lightblue")
}

#-------------------------------------------------------#
#           6. Correlation matrix
#-------------------------------------------------------#
cor_matrix <- cor(Cat_data)
print(cor_matrix)

#-------------------------------------------------------#
#           7. Visualizing correlations with heatmap
#-------------------------------------------------------#
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

#-------------------------------------------------------#
#          Normalize the Dataset
#-------------------------------------------------------#
Cat_scaled_data <- scale(Cat_data)

#--------------------------------------------------------------#
#           8. Skewness and kurtosis for each feature
#--------------------------------------------------------------#
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

#----------------------------------------------------------#
#             Combined skewness and kurtosis
#----------------------------------------------------------#
combined_data <- as.vector(as.matrix(Cat_scaled_data))
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
for (i in 1:ncol(Cat_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Cat_scaled_data[,i])$out
  print(paste("Feature:", colnames(Cat_scaled_data)[i]))
  print(paste("Outliers:", length(outliers1[[i]])))
  total_outliers1 <- total_outliers1 + length(outliers1[[i]])
  outlier_rows1 <- unique(c(outlier_rows1, outliers1[[i]]))  # Collect unique rows with outliers
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
outliers_matrix <- apply(Cat_scaled_data, 2, detect_outlier)
outliers_matrix

# Number of rows that have at least one TRUE (outlier)
rows_with_outliers <- rowSums(outliers_matrix) > 0
num_rows_with_outliers <- sum(rows_with_outliers)

# Result
print(paste("Number of rows with at least one outlier:", num_rows_with_outliers))
