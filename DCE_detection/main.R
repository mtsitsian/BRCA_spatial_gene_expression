########Main script of the DCE (COD) detection pipeline########
###IMPORTANT!!! Start the analysis by installing and loading all the packages in the "install_packages.R" script
###Continue by downloading and preparing the data via the "TCGA_Data_Mining.R", "Initial filtering.R", and the "data_preparation.R" scripts, they are all located in the DCE_detection folder in the github repo

# load required libraries
library(parallel) #for parallel computation
library(doParallel) #for parallel computation
library(features) #for detecting local minima and maxima to detect DCE 
library(feather) #to save files in a lightweight format


# source all the rest required functions and prepare the data
filenames <- c('scripts/calculate_correlations.R',
               'scripts/permutation_test_for_correlations.R',
               'scripts/compute_binsignal.R',
               'scripts/detect_DCEs.R',
               'scripts/merge_DCEs.R')
sapply(filenames, source, echo = T)

# Calculate correlation matrix
chrs <- unique(bins$chromosome)

#cor_final_Healthy <- calc_cor_final(bin_counts = bin_counts_Healthy)
cor_final_NP <- calc_cor_final(bin_counts = bin_counts_NP)
cor_final_CIN <- calc_cor_final(bin_counts = bin_counts_CIN)
cor_final_EBV <- calc_cor_final(bin_counts = bin_counts_EBV)
cor_final_MSI <- calc_cor_final(bin_counts = bin_counts_MSI)
cor_final_GS <- calc_cor_final(bin_counts = bin_counts_GS)

############Extract the correlation matrices to your working directory to use them for the plotgardener heatmaps#################

# write_feather(cor_final_NP, "cor_final_NP.feather")
# read_feather("cor_final_NP.feather")

dir.create("correlation_matrices")

write_feather(cor_final_NP, "cor_final_NP.feather")
write_feather(cor_final_EBV, "cor_final_EBV.feather")
write_feather(cor_final_CIN, "cor_final_CIN.feather")
write_feather(cor_final_MSI, "cor_final_MSI.feather")
write_feather(cor_final_GS, "cor_final_GS.feather")



# Calculate permutation scores
#WARNING!!!! This script is made for windows OS!!! 
#For parallel computing in linux OS, use the parameter "mc.cores = 4" to make the analysis quicker, like the example below

#EXAMPLE FOR LINUX OS
#system.time(
#  permutation_scores_Healthy <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor,
#                                         bin_counts = bin_counts_Healthy, cor_final = cor_final_Healthy, mc.cores = 4)
#)

#EXAMPLE FOR WINDOWS OS, THE PERMUTATION NUMBER OUTPUT DOES NOT GET PRINTED BUT IT IS SIGNIFICANTLY FASTER!!!
# Set the number of cores to be used
num_cores <- 4

# Initialize the parallel backend
cl <- makeCluster(num_cores)

# Register the parallel backend
registerDoParallel(cl)

# Run the parallel computation
# Run the parallel computation
system.time(
  permutation_scores_NP <- foreach(chr = chrs) %dopar% {
    permutations(chr, n_of_permutations = 1000, correlation = cor, bin_counts = bin_counts_NP, cor_final = cor_final_NP)
  }
)


system.time(
  permutation_scores_CIN <- foreach(chr = chrs) %dopar% {
    permutations(chr, n_of_permutations = 1000, correlation = cor, bin_counts = bin_counts_CIN, cor_final = cor_final_CIN)
  }
)

system.time(
  permutation_scores_EBV <- foreach(chr = chrs) %dopar% {
    permutations(chr, n_of_permutations = 1000, correlation = cor, bin_counts = bin_counts_EBV, cor_final = cor_final_EBV)
  }
)

system.time(
  permutation_scores_MSI <- foreach(chr = chrs) %dopar% {
    permutations(chr, n_of_permutations = 1000, correlation = cor, bin_counts = bin_counts_MSI, cor_final = cor_final_MSI)
  }
)


system.time(
  permutation_scores_GS <- foreach(chr = chrs) %dopar% {
    permutations(chr, n_of_permutations = 1000, correlation = cor, bin_counts = bin_counts_GS, cor_final = cor_final_GS)
  }
)


# Stop the parallel backend
stopCluster(cl)




# compute intra cod coexpression
compute_intra_cod_coexp <- function(cods, chr, thres = 0.05, group){
  bin_counts <- get(x = paste0("bin_counts_", group))
  map <- calc_cor_matrix(bins = bins, bin_counts = bin_counts, chr = chr, cor_method = "spearman", correlation = cor)
  rnames <- rownames(map)
  cnames <- colnames(map)
  if (thres < 1){
    permutation_scores <- get(x = paste0("permutation_scores_", group))
    x <- get(x = paste0("cor_final_", group))
    x <- x[which(x$chrom == chr),]
    x$cor[which(permutation_scores[[chr]] > thres)] <- 0
    map <- diag(1, nrow = nrow(map), ncol = ncol(map))
    map[upper.tri(map)] <- x$cor
    map <- t(map)
    map[upper.tri(map)] <- x$cor
    rownames(map) <- rnames
    colnames(map) <- cnames
  }
  average_intra_cod <- numeric()
  for (i in 1:nrow(cods)){
    cod <- map[cods[i,"START"]:cods[i,"END"], cods[i,"START"]:cods[i,"END"]]
    average_intra_cod[i] = mean(cod[upper.tri(cod)])
  }
  return(average_intra_cod)
}

# check intra cod coexpression using different window sizes to define cods
test_for_window <- function(windows, chr, group){
  intra_cod <- list()
  for (w in windows){
    print(w)
    cods <- codfinder(chr = chr, group = group, 
                      window.size = w, 
                      bin_signal_thres = 0.06798727, #Bin signal threshhold is the average binsignal of the control group! Don't forget to change that value! Read below!
                      perm_pvalue_thres = 0.05)
    cods <- cods$CODs
    cods <- cod_merge(cods, interfere = 2, group = group, chr = chr, permut_thres = 0.05)
    intra_cod[[as.character(w)]] <- compute_intra_cod_coexp(cods, chr, thres = 0.05, group)
  }
  return(intra_cod)
}

check_w <- test_for_window(3:19, "chr1", "pretreat")
str(check_w)
sapply(check_w, function(x) length(x))
sapply(check_w, function(x) summary(x))


# Find CODs for all chrs and all groups
#run the the make_cod_list function for the control group to find the average bin_signal and use it as the bin signal threshold (RUN THE FUNCTION AGAIN)

bin_signals <- list()
CODs <- list()
make_COD_list <- function(group, window.size, chrs){
  print(group)
  cods <- list()
  bin_signals <- get('bin_signals', envir = .GlobalEnv)
  bin_signals[[group]] <- list()
  for (chr in chrs){
    print(chr)
    c <- codfinder(chr, group, window.size, 
                   bin_signal_thres = 	0.06798727, #Bin signal threshhold is the average binsignal of the control group! Don't forget to change that value! Read below!
                   perm_pvalue_thres = 0.05)
    cods[[chr]] <- c$CODs
    bin_signals[[group]][[chr]] <- c$bin_signal
  }
  assign('bin_signals', value = bin_signals, envir = .GlobalEnv)
  cods
}
#CODs[['Healthy']] <- make_COD_list('Healthy', 3, chrs)

## Merge cods with a gap of at most two bins between them
merge_cods <- function(cod_list, group){
  for (i in 1:length(cod_list)){
    print(names(cod_list)[[i]])
    cod_list[[i]] <- cod_merge(cod_list[[i]], 2, group, names(cod_list)[[i]], 0.05)
  }
  cod_list
}
#CODs$Healthy <- merge_cods(CODs$Healthy, 'Healthy')


#Run the make_COD_list script for the rest of the subtypes in the analysis
CODs[['NP']] <- make_COD_list('NP', 3, chrs)
CODs$NP <- merge_cods(CODs$NP, 'NP')

CODs[['CIN']] <- make_COD_list('CIN', 3, chrs)
CODs$CIN <- merge_cods(CODs$CIN, 'CIN')

CODs[['EBV']] <- make_COD_list('EBV', 3, chrs)
CODs$EBV <- merge_cods(CODs$EBV, 'EBV')

CODs[['MSI']] <- make_COD_list('MSI', 3, chrs)
CODs$MSI <- merge_cods(CODs$MSI, 'MSI')

CODs[['GS']] <- make_COD_list('GS', 3, chrs)
CODs$GS <- merge_cods(CODs$GS, 'GS')

########Average Binsignal########
#calculate the average binsignal for the control group with the script below (in the functions above, the control group is called "Healthy", just for demonstration purposes, as a comment)
#Then, run the previous scripts again (test_for_window and make_COD_list), this time by changing the "bin_signal_thres" argument to the average binsignal value for the control

average_binsignal <- 
  apply(X = sapply(bin_signals, function(x) sapply(x, function(y) mean(y[,1]))), 
        MARGIN = 2, 
        FUN = mean)

average_binsignal %<>% 
  as.data.frame %>% 
  rownames_to_column('group')
colnames(average_binsignal)[2] <- "Average_binsignal"
average_binsignal$group <- factor(average_binsignal$group, levels = c("NP", "EBV", "CIN", "MSI", "GS"))


## Find the cod coordinates
coCODs <- function(cods){
  c_CODs <- list()
  for (group in names(cods)){
    c_CODs[[group]] <- list()
    for (chr in names(cods[[1]])){
      c_CODs[[group]][[chr]] <- matrix(0, nrow(cods[[group]][[chr]]), ncol(cods[[group]][[chr]]), dimnames = list(NULL, c("START", "END")))
      bin_starts <- bin_starts <- as.numeric(rownames(bin_signals[[group]][[chr]]))
      
      c_CODs[[group]][[chr]][,"START"] <- bin_starts[cods[[group]][[chr]][,"START"]]
      c_CODs[[group]][[chr]][,"END"] <- bin_starts[cods[[group]][[chr]][,"END"]]+9999
    }
  }
  return(c_CODs)
} #find the cod coordinates
c_CODs <- coCODs(CODs)
