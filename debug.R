## do 50 times
simulation_times = 50
set.seed(123)
simulationDF = data.frame()
time = matrix(NA, nrow = simulation_times, ncol = 6)
for (i in 1:simulation_times){
  timeL = rep(NA, 6) # DCATS w/wto simularity matrix, fisher, speckle, diffcyt, scDC
  ## simulation
  simulation = simulator_fastMNN(totals1 = runif(rep1, 1500, 2000), totals2 = runif(rep2, 1500, 2500), probC1, probC2, setresolu)
  if (is.na(simulation)[1] == FALSE){
    ## data convert
    numb_cond1 = simulation$dfRes %>% 
      filter(condition == "Cond1") %>% 
      group_by(clusterRes, batch) %>% 
      summarise(n = n()) %>% 
      pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
      column_to_rownames(var = "batch")
    numb_cond2 = simulation$dfRes %>% 
      filter(condition == "Cond2") %>% 
      group_by(clusterRes, batch) %>% 
      summarise(n = n()) %>% 
      pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
      column_to_rownames(var = "batch")
    conf.mat<-table(Idents(simulation$integratedSamples), simulation$trueLabels)
    ## test
    ## DCATS---betabinLRT
    t1start = Sys.time()
    betabin_noBC = betabinLRT_rw(as.matrix(numb_cond1), as.matrix(numb_cond2), bias_corr = FALSE)
    timeL[1] = Sys.time() - t1start
    ##  DCATS---betabinLRT with Bias correction from similarity matrix
    simil_matI = get_similarity_mat(cluster_num, confuse_rate=0)
    simil_matU = get_similarity_mat(cluster_num, confuse_rate=0.1)
    betabin_wBC_I = betabinLRT_rw(as.matrix(numb_cond1), as.matrix(numb_cond2), bias_corr = TRUE, similarity_mat = simil_matI)
    betabin_wBC_U = betabinLRT_rw(as.matrix(numb_cond1), as.matrix(numb_cond2), bias_corr = TRUE, similarity_mat = simil_matU)
    ## DCATS---betabin with similarity matrix (backward-forward)
    t2start = Sys.time()
    betabin_wMI = dcats_betabin(as.matrix(numb_cond1), as.matrix(numb_cond2), similarity_mat = simil_matI, n_samples = 100, binom_only = FALSE)
    timeL[2] = Sys.time() - t2start
    betabin_wMU = dcats_betabin(as.matrix(numb_cond1), as.matrix(numb_cond2), similarity_mat = simil_matU, n_samples = 100, binom_only = FALSE)
    bin_wMI = dcats_betabin(as.matrix(numb_cond1), as.matrix(numb_cond2), similarity_mat = simil_matI, n_samples = 100)
    bin_wMU = dcats_betabin(as.matrix(numb_cond1), as.matrix(numb_cond2), similarity_mat = simil_matU, n_samples = 100)
    ## Fisher's exact test
    t3start = Sys.time()
    fisher_pvals = getFisher(as.matrix(numb_cond1), as.matrix(numb_cond2))
    timeL[3] = Sys.time() - t3start
    ## speckle
    t4start = Sys.time()
    speckleRes = propeller(clusters = simulation$dfRes$clusterRes, sample = simulation$dfRes$batch, group = simulation$dfRes$condition) %>% 
      dplyr::rename(cluster = BaselineProp.clusters,
                    speckle_pvals = FDR) %>%
      dplyr::select(cluster, speckle_pvals)
    timeL[4] = Sys.time() - t4start
    ## diffcyt
    t5start = Sys.time()
    diffcytP = getDiffcyt(numb_cond1, numb_cond2, simulation$integratedSamples)
    timeL[5] = Sys.time() - t5start
    ## scDC
    t6start = Sys.time()
    res_scDC <- scDC_noClustering(cellTypes = simulation$dfRes$clusterRes, simulation$dfRes$batch, calCI = TRUE, calCI_method = c("percentile", "BCa", "multinom"),nboot = 1000, verbose = FALSE)
    res_GLM <- fitGLM(res_scDC, c(rep("cond1",rep1*cluster_num),rep("cond2",rep2*cluster_num)), pairwise = FALSE, fixed_only = TRUE, verbose = FALSE)
    timeL[6] = Sys.time() - t6start
    scDCRes_temp = summary(res_GLM$pool_res_fixed)
    scDCRes = scDCRes_temp[c(cluster_num+1,(dim(scDCRes_temp)[1]-cluster_num+2):dim(scDCRes_temp)[1]),]
    ## results
    cluster_map = apply(conf.mat, 2, which.max) %>%
      plyr::mapvalues(from = c(1:cluster_num), to = c("A", "B", "C", "D", "E", "F", "G ","H", "I", "J", "K", "L")[1:cluster_num]) %>% 
      as_tibble() %>% 
      mutate(truth = c("N", "N", "P", "P", "N", "P", "N", "P", "N", "N", "P", "P")) %>%
      dplyr::rename(cluster = value) %>% 
      arrange(cluster)
    all_res = cluster_map %>% 
      mutate(betabin_noBC_pvals = betabin_noBC$pvals,
             betabin_wBCI_pvals = betabin_wBC_I$pvals,
             betabin_wBCU_pvals = betabin_wBC_U$pvals,
             betabin_wMI_pvals = betabin_wMI$pvals,
             betabin_wMU_pvals = betabin_wMU$pvals,
             bin_wMI_pvals = bin_wMI$pvals,
             bin_wMU_pvals = bin_wMU$pvals,
             fisher_pvals = fisher_pvals,
             scDC_pvals = scDCRes$p.value) %>%
      merge(speckleRes, by = "cluster") %>% 
      merge(diffcytP, by = "cluster")
    simulationDF = rbind(simulationDF, all_res)
    time[i,] = timeL
  }
}

multinom_EM <- function(X, simMM, min_iter=10, max_iter=1000,
                        logLik_threshold=1e-2, verbose=TRUE) {
  # Be very careful on the shape of simMM; rowSums(simMM) = 1
  K = ncol(simMM)
  
  # initialization
  mu = sample(K)
  mu = mu / sum(mu)
  Z = matrix(NA, K, K)
  logLik_old <- logLik_new <- log(mu %*% simMM) %*% X
  
  if (verbose) {
    print(paste("Iteration 0 logLik:", round(logLik_new, 3)))
  }
  
  for (it in seq_len(max_iter)) {
    ## E step: expectation of count from each component
    # Z = (simMM * mu)
    # Z = t(t(Z) / colSums(Z))
    
    for (i in seq(K)) {
      for (j in seq(K)){
        Z[i, j] = simMM[i, j] * mu[i] / sum(mu * simMM[, j])
      }
    }
    
    ## M step: maximizing likelihood
    # mu = c(X %*% Z) # this is wrong
    
    ## v2
    mu = c(Z %*% X)
    mu = mu / sum(mu)
    
    ## Check convergence
    logLik_new <- log(mu %*% simMM) %*% X
    # sum(X * log(mu %*% t(simMM)))
    if (it > min_iter && logLik_new - logLik_old < logLik_threshold) {
      break
    } else {
      logLik_old <- logLik_new
    }
    if (verbose) {
      print(paste("Iteration", it, "logLik:", round(logLik_new, 3)))
    }
    
  }
  
  ## return values
  list("mu" = mu, "logLik" = logLik_new,
       "simMM" = simMM, "X" = X, "X_prop" = X / sum(X),
       "predict_X_prop" = mu %*% simMM)
}
