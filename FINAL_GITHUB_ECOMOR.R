#################################################################
#                       Ecomor DATASET                         ##
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

#---------------------------------------------------#
#             Load the Ecomor Dataset
#---------------------------------------------------#
data(ecomor, package = "ade4")
ecomor                            
help(ecomor)
Ecomor_data <- as.matrix(ecomor$morpho)
true_lab_Ecomor=as.numeric(ecomor$taxo$Ordre)

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
#         Scaling the data
#--------------------------------------#
Ecomor_scaled_data <- scale(Ecomor_data)

#--------------------------------------#
#      Plotting true labels
#--------------------------------------#
fviz_cluster(list(data = Ecomor_scaled_data, cluster = true_lab_Ecomor), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#--------------------------------------#
#    Number of rows and columns
#--------------------------------------#
P <- nrow(Ecomor_scaled_data)
D <- ncol(Ecomor_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 7  #maximum no. of clusters

#--------------------------------------#
#      Set up data for JAGS
#--------------------------------------#
data_list <- list(
  x = Ecomor_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Ecomor_scaled_data))^(1 / (D - 1)),
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
  v_g <- 8.4240999999
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
       mu_g[j, g] ~ dgamma(0.595908182,16.19557011)        
       lambda_g[j, g] ~ dgamma(3.367255795,6.34330397) 
    }
    alpha[g] ~ dgamma(6,14) 
   }
   pi_g[1:G] ~ ddirch(alpha/sum(alpha))
}
"
#-----------------------------------------------#
#            Write the model to a file
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
run_time_Ecomor <- toc()
codaSamples_carni = as.mcmc.list( runJagsOut )
summaryChains_carni <- summary(codaSamples_carni)

#----------------------------------------------#
#               Diagnostics
#----------------------------------------------#
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_carni , parName="mu_g[2,2]" )
graphics.off()

#---------------------------------------------#
#             Summary Chains
#---------------------------------------------#
z_mode_Ecomor <- apply(matrix(summaryChains_carni$statistics[(1+8*G+8*G+P+1):nrow(summaryChains_carni$statistics),1], P, G),1, which.max)
z_mode_Ecomor

#--------------------------------------------#
#         Plotting Predicted labels
#--------------------------------------------#
fviz_cluster(list(data = Ecomor_scaled_data, cluster = z_mode_Ecomor), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")
fviz_cluster(list(data = Ecomor_scaled_data, cluster = true_lab_Ecomor), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

plot(Ecomor_scaled_data, col= z_mode_Ecomor, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ecomor), col=unique(z_mode_Ecomor), pch=16, title="Cluster")
table( true_lab_Ecomor , z_mode_Ecomor)
table( true_lab_Ecomor , true_lab_Ecomor)

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
relabeled <- relabel_clusters_optimal(z_mode_Ecomor, true_lab_Ecomor)
z_mode_Ecomor <- relabeled$remapped_labels

cat("  Confusion matrix (after relabeling):\n")
print(relabeled$confusion)

# Cohen's Kappa
kappa_result <- kappa2(data.frame(rater1 = true_lab_Ecomor, rater2 = z_mode_Ecomor))
kappa_scores <- kappa_result$value
cat("Kappa:", kappa_result$value, "\n")

plot(Ecomor_scaled_data, col= z_mode_Ecomor, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Ecomor), col=unique(z_mode_Ecomor), pch=16, title="Cluster")
table( true_lab_Ecomor , true_lab_Ecomor)
table( true_lab_Ecomor , z_mode_Ecomor)
kappa2(data.frame(rater1 = true_lab_Ecomor, rater2 = z_mode_Ecomor))

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
    CPU_RUN_TIME = run_time_Ecomor$callback_msg  # Store the runtime
  )
  return(results_list)
}

#--------------------#
# Call the function
#--------------------#
dp_Gmvlg_metricsEcomor_1 <- calculate_dp_Gmvlg_clustering_metrics(data_1 = Ecomor_scaled_data, true_clusters = true_lab_Ecomor, z_mode = z_mode_Ecomor)

#--------------------#
#     Results 
#--------------------#
for (dp_Gmvlg in names(dp_Gmvlg_metricsEcomor_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsEcomor_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
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
kmeans_result_Ecomor <- eclust(Ecomor_scaled_data, "kmeans", G, graph = FALSE)
kmeans_time_Ecomor <- toc()
kmeans_clusters_Ecomor <- kmeans_result_Ecomor$cluster
table(true_lab_Ecomor, kmeans_clusters_Ecomor)

# relabeling function
relabeled_kmeans <- relabel_clusters_optimal(kmeans_clusters_Ecomor, true_lab_Ecomor)
kmeans_clusters_Ecomor <- relabeled_kmeans$remapped_labels

#-------------------------------------------------------------------------#
#                            CLARA clustering
#-------------------------------------------------------------------------#
tic("CLARA Runtime")
clara_result_Ecomor <- eclust(Ecomor_scaled_data, "clara", G, graph = FALSE)
clara_time_Ecomor <- toc()
clara_clusters_Ecomor <- clara_result_Ecomor$cluster
table(true_lab_Ecomor, clara_clusters_Ecomor)

# relabeling function
relabeled_clara <- relabel_clusters_optimal(clara_clusters_Ecomor, true_lab_Ecomor)
clara_clusters_Ecomor <- relabeled_clara$remapped_labels

#-------------------------------------------------------------------------#
#                            PAM clustering
#-------------------------------------------------------------------------#
tic("PAM Runtime")
pam_result_Ecomor <- eclust(Ecomor_scaled_data, "pam", G, graph = FALSE)
pam_time_Ecomor <- toc()
pam_clusters_Ecomor <- pam_result_Ecomor$cluster
table(true_lab_Ecomor, pam_clusters_Ecomor)

# relabeling function
relabeled_pam <- relabel_clusters_optimal(pam_clusters_Ecomor, true_lab_Ecomor)
pam_clusters_Ecomor <- relabeled_pam$remapped_labels

#-------------------------------------------------------------------------#
#                       Hierarchical clustering
#-------------------------------------------------------------------------#
tic("Hierarchical Runtime")
hclust_result_Ecomor <- hclust(dist(Ecomor_scaled_data), method = "ward.D2")
hclust_time_Ecomor <- toc()
hclust_clusters_Ecomor <- cutree(hclust_result_Ecomor, k = G)
table(true_lab_Ecomor, hclust_clusters_Ecomor)

# relabeling function
relabeled_hclust <- relabel_clusters_optimal(hclust_clusters_Ecomor, true_lab_Ecomor)
hclust_clusters_Ecomor <- relabeled_hclust$remapped_labels

#------------------------------------------------------------------------#
#                       Model-based clustering
#------------------------------------------------------------------------#
tic("Model-based Runtime")
mclust_result_Ecomor <- Mclust(Ecomor_scaled_data, G = G)
mclust_time_Ecomor <- toc()
summary(mclust_result_Ecomor)
mclust_clusters_Ecomor <- mclust_result_Ecomor$classification
table(true_lab_Ecomor, mclust_clusters_Ecomor)

# relabeling function
relabeled_mclust <- relabel_clusters_optimal(mclust_clusters_Ecomor, true_lab_Ecomor)
mclust_clusters_Ecomor <- relabeled_mclust$remapped_labels

#-------------------------------------------------------------------------------------#
# Bayesian Clustering using Dirichlet Process Multivariate Normal Distribution (DPMVN)
#-------------------------------------------------------------------------------------#
# Create Dirichlet Process object with adjusted concentration parameter
set.seed(23)
dp <- DirichletProcessMvnormal(Ecomor_scaled_data, alphaPriors =  c(2,50))
tic("DPMVN Runtime")
dp <- Fit(dp, its = 1000)
DPMVN_time_Ecomor <- toc()

# Extract clusters 
dpMVN_clusters_Ecomor <- as.numeric(dp$clusterLabels)
table(true_lab_Ecomor, dpMVN_clusters_Ecomor)

#old way
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

#--------------------------------------------------------------------------------------#
#                 All Classical Except Genie and one Bayesian (DPMVN) Results Together
#--------------------------------------------------------------------------------------#
calculate_clustering_metricsEcomor_2 <- function(dataEcomor_2, true_clusters_Ecomor, estimated_clusters_Ecomor_list, times_list) {
  
  results_list <- list()

  for (method_name in names(estimated_clusters_Ecomor_list)) {
    clusters <- estimated_clusters_Ecomor_list[[method_name]]

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

  return(results_list)
}

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

#-------------------#
# Call the function
#-------------------#
clustering_metricsEcomor_2 <- calculate_clustering_metricsEcomor_2(dataEcomor_2 = Ecomor_scaled_data, true_clusters_Ecomor = true_lab_Ecomor, estimated_clusters_Ecomor_list = cluster_result_Ecomor, times_list = times_list_Pu)

#-------------------#
#     Results 
#-------------------#
for (method_name in names(clustering_metricsEcomor_2)) {
  cat("\n", method_name, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(clustering_metricsEcomor_2[[method_name]]$Table)
  cat("Kappa:", clustering_metricsEcomor_2[[method_name]]$Kappa, "\n")
  cat("Adjusted Rand Index:", clustering_metricsEcomor_2[[method_name]]$ARI, "\n")
  cat("Normalized Mutual Information:", clustering_metricsEcomor_2[[method_name]]$NMI, "\n")
  cat("", clustering_metricsEcomor_2[[method_name]]$CPU_RUN_TIME, "\n")
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

k <- length(unique(true_lab_Ecomor))

for (i in seq_along(gini_values)) {
  gini <- gini_values[i]
  cat("\nTrying gini threshold:", gini, "\n")
  
  tic("Genie Runtime")
  res <- try(gclust(Ecomor_scaled_data, gini_threshold = gini), silent = TRUE)
  run_time_Ecomor<- toc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    cat("Error or invalid clustering result\n")
    kappa_scores[i] <- NA
  } else {
    cluster_labels <- cutree(res, k = k)

    relabeled <- relabel_clusters_optimal(cluster_labels, true_lab_Ecomor)
    aligned_labels <- relabeled$remapped_labels
    
    cat("  Confusion matrix (after relabeling):\n")
    print(relabeled$confusion)

    kappa_result <- kappa2(data.frame(rater1 = true_lab_Ecomor, rater2 = aligned_labels))
    kappa_scores[i] <- kappa_result$value
    cat("Kappa:", kappa_result$value, "\n")
  }
}

# final Kappa scores
cat("\nFinal Kappa scores by Gini threshold:\n")
print(kappa_scores)

#-------------------------------------------------------#
#           Running Genie with best Gini Coefficient
#-------------------------------------------------------#
tic("Genie Runtime")
res <- gclust(Ecomor_scaled_data, gini_threshold = 0.901)
run_time_Ecomor_Genie <- toc()

k <- length(unique(true_lab_Ecomor))

cluster_labels <- cutree(res, k = k)

# relabeling function
relabeled_Genie <- relabel_clusters_optimal(cluster_labels, true_lab_Ecomor)
cluster_labels <- relabeled_Genie$remapped_labels

table( true_lab_Ecomor , aligned_labels)
kappa2(data.frame(rater1 = true_lab_Ecomor, rater2 = aligned_labels))

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
    CPU_RUN_TIME = run_time_Ecomor_Genie$callback_msg  # Store the runtime
  )
  return(results_list)
}

#------------------#
# Call the function
#------------------#
Genie_metricsEcomor_1 <- calculate_Genie_clustering_metrics(data_1 = Ecomor_scaled_data, true_clusters = true_lab_Ecomor, cluster_labels = cluster_labels)

#------------------#
#     Results 
#------------------#
for (Genie in names(Genie_metricsEcomor_1)) {
  cat("\n", Genie, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(Genie_metricsEcomor_1[[Genie]]$Table)
  cat("Kappa:", Genie_metricsEcomor_1[[Genie]]$Kappa, "\n")
  cat("Adjusted Rand Index:", Genie_metricsEcomor_1[[Genie]]$ARI, "\n")
  cat("NMI:", Genie_metricsEcomor_1[[Genie]]$NMI, "\n")
  cat(" CPU RUN TIME:", Genie_metricsEcomor_1[[Genie]]$CPU_RUN_TIME, "\n")
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
summary(Ecomor_data)

#-------------------------------------------------------#
#           2. Check for missing values
#-------------------------------------------------------#
sum(is.na(Ecomor_data))

#-------------------------------------------------------#
#           3. Distribution of each feature
#-------------------------------------------------------#
par(mfrow=c(2, 4)) 

for (i in 1:ncol(Ecomor_data)) {
  hist(Ecomor_data[, i], main=colnames(Ecomor_data)[i], xlab=colnames(Ecomor_data)[i], col="lightblue", border="black")
}

#-------------------------------------------------------#
#           4. Pairwise scatter plots
#-------------------------------------------------------#
pairs(Ecomor_data, main="ECOMOR Data - Pairwise Scatter Plots", pch=21, bg=c("red", "green3")[true_lab_Ecomor])

#-------------------------------------------------------#
#           5. Boxplots for each feature 
#-------------------------------------------------------#
# Normalized data to the original data frame
pid_normalized <- as.data.frame(Ecomor_data)
pid_normalized$diabetes <- ecomor$taxo$Ordre

par(mfrow=c(2, 4)) 

for (i in 1:ncol(Ecomor_data)) {
  boxplot(Ecomor_data[, i] ~ ecomor$taxo$Ordre, main=colnames(Ecomor_data)[i], xlab="Order", ylab=colnames(Ecomor_data)[i], col="lightblue")
}

#-------------------------------------------------------#
#           6. Correlation matrix
#-------------------------------------------------------#
cor_matrix <- cor(Ecomor_data)
print(cor_matrix)

#-------------------------------------------------------#
#           7. Visualizing correlations with heatmap
#-------------------------------------------------------#
par(mfrow=c(1, 1))
corrplot(cor_matrix, method="color", addCoef.col = "black", tl.cex=0.8, number.cex=0.7)

#-------------------------------------------------------#
#          Normalize the Dataset
#-------------------------------------------------------#
Ecomor_scaled_data <- scale(Ecomor_data)

#--------------------------------------------------------------#
#           8. skewness and kurtosis for each feature
#--------------------------------------------------------------#
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

#----------------------------------------------------------#
#             Combined skewness and kurtosis
#----------------------------------------------------------#
combined_data <- as.vector(as.matrix(Ecomor_scaled_data))
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
for (i in 1:ncol(Ecomor_scaled_data)) {
  outliers1[[i]]<- boxplot.stats(Ecomor_scaled_data[,i])$out
  print(paste("Feature:", colnames(Ecomor_scaled_data)[i]))
  print(paste("Outliers:", length(outliers1[[i]])))
  total_outliers1 <- total_outliers1 + length(outliers1[[i]])
  outlier_rows1 <- unique(c(outlier_rows1, outliers1[[i]]))  
}

# Total number of outliers 
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
outliers_matrix <- apply(Ecomor_scaled_data, 2, detect_outlier)
outliers_matrix

# Number of rows that have at least one TRUE (outlier)
rows_with_outliers <- rowSums(outliers_matrix) > 0
num_rows_with_outliers <- sum(rows_with_outliers)

# Result
print(paste("Number of rows with at least one outlier:", num_rows_with_outliers))
