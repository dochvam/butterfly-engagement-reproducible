#### Helper GAM functions####

process_spec_diff <- function(this_species, this_region, hex_data, target_species, write_inds = T, 
                              verbose = F, overwrite = F, timeout = Inf, maxK = 20,
                              output_path = "intermediate/ind_results/",
                              warn_path = "intermediate/warnings/") {
  
  withTimeout({
    this_name_clean <- gsub(" ", "_", this_species)
    
    flog.appender(appender.file(paste0(warn_path, this_name_clean, "_warning.log")))
    this_output_file <- paste0(output_path, "res_", this_name_clean, "_", this_region, ".csv")
    if (file.exists(this_output_file) && !overwrite) {
      return(NULL)
    }
    
    target_df <- data.frame(
      species = this_species,
      common_name = this_name_clean,
      difference = NA,
      diff_SE = NA,
      k1 = NA, k2 = NA, k3 = NA,
      inat_medPred = NA,
      inat_SE = NA,
      inat_medResp = NA,
      inat_GAM_pval = NA,
      ebud_medPred = NA,
      ebud_SE = NA,
      ebud_medResp = NA,
      ebud_GAM_pval = NA
    )
    
    if (verbose) {
      cat("Processing", this_species, "-", this_region, "\n")
    }
    
    this_dat_ebud <- hex_data %>% 
      filter(species == this_species, dataset == "eButterfly") %>% 
      mutate(
        unscaled_x = x, unscaled_y = y,
        x = scale(x), y = scale(y)
      ) %>% 
      rename(total = ncelltotal, success = n)
    
    this_dat_inat <- hex_data %>% 
      filter(species == this_species, dataset == "iNat") %>% 
      mutate(
        unscaled_x = x, unscaled_y = y,
        x = scale(x), y = scale(y)
      ) %>% 
      rename(total = ncelltotal, success = n)
    
    nyear <- length(unique(this_dat_ebud$year))
    ngridcell <- length(unique(this_dat_ebud$gridcell))
    
    spacetime_dim <- c(sqrt(ngridcell), sqrt(ngridcell), nyear)
    # density_scalar <- (nrow(this_dat_ebud) / prod(spacetime_dim))^(1/3)
    k_spatial <- floor(sqrt(nrow(this_dat_ebud) / (nyear/2)))
    if (any(k_spatial < 3)) {
      k_spatial[k_spatial < 3] <- 3
      k_time <- floor(nrow(this_dat_ebud) / k_spatial^2)
    } else {
      k_time <- floor(nyear/2)
    }
    
    kvec_bounded <- c(k_spatial, k_spatial, k_time)
    kvec <- kvec_bounded
    kvec[kvec > maxK] <- maxK
    
    target_df$kvec_at_max <- all(kvec == kvec_bounded)
    
    cat(paste0("Data summary: ", ngridcell, " hexes over ", nyear, " years...\n"))
    
    cat("... Fitting eButterfly GAM ...\n")
    this_fit_ebud <- gam(cbind(success, total - success) ~ 
                           ti(x, y, k = kvec[1:2]) + 
                           ti(year, k = kvec[3]) + 
                           ti(x, y, year, k = kvec),
                         family = quasibinomial(link = "logit"), 
                         optimizer = c("outer", "bfgs"),
                         control = list(maxit = 100000),
                         data = this_dat_ebud)
    
    cat("... Fitting iNaturalist GAM ...\n")
    this_fit_inat <- gam(cbind(success, total - success) ~ 
                           ti(x, y, k = kvec[1:2]) + 
                           ti(year, k = kvec[3]) + 
                           ti(x, y, year, k = kvec), 
                         family = quasibinomial(link = "logit"),
                         optimizer = c("outer", "bfgs"),
                         control = list(maxit = 100000),
                         data = this_dat_inat)
    
    # Save a data frame of results.
    # Each row: one dataset/hex-year, w/ a prediction & uncertainty
    locs <- this_dat_ebud[, c("x", "y", "longitude", "latitude", "year")]
    
    ebud_pred_oi <- predict(this_fit_ebud, newdata = locs, se.fit = TRUE)
    inat_pred_oi <- predict(this_fit_inat, newdata = locs, se.fit = TRUE)
    ebud_results <- locs %>% 
      mutate(dataset = "eButterfly", 
             index = as.numeric(ebud_pred_oi$fit),
             se = as.numeric(ebud_pred_oi$se.fit))
    inat_results <- locs %>% 
      mutate(dataset = "iNaturalist", 
             index = as.numeric(inat_pred_oi$fit),
             se = as.numeric(inat_pred_oi$se.fit))
    
    all_pred_inds <- bind_rows(ebud_results, inat_results) %>% 
      select(-x, -y)
    write_csv(all_pred_inds,
              paste0(output_path, "predicted_indices_", this_species, "_", this_region, ".csv")
    )
    
    
    check_ebud = k.check(this_fit_ebud)
    check_inat = k.check(this_fit_inat)
    
    target_df$k1 = kvec[1]
    target_df$k2 = kvec[2]
    target_df$k3 = kvec[3]
    target_df$ebud_GAM_pval = min(check_ebud[, 4])
    target_df$inat_GAM_pval = min(check_inat[, 4])
    
    #FYI check_ebud[2] is edf, check_ebud[3] is kindex
    
    
    rmvn <- function(n,mu,sig) { ## MVN random deviates
      L <- mroot(sig);m <- ncol(L);
      t(mu + L%*%matrix(rnorm(m*n),m,n)) 
    }
    
    Xp_ebud <- predict(this_fit_ebud ,type="lpmatrix") 
    Xp_inat  <- predict(this_fit_inat ,type="lpmatrix") 
    
    br_ebud <- rmvn(10000, coef(this_fit_ebud),this_fit_ebud$Vp) ## 1000 replicate param. vectors
    br_inat <- rmvn(10000, coef(this_fit_inat),this_fit_inat$Vp) ## 1000 replicate param. vectors
    
    res1 <- rep(0,10000)
    res2 <- rep(0,10000)
    res3 <- rep(0,10000)
    for (i in 1:10000) { 
      pr_ebd  <- Xp_ebud %*% br_ebud[i,] ## replicate predictions
      pr_inat <- Xp_inat %*% br_inat[i,] ## replicate predictions
      res1[i] <- median(pr_ebd) ## median eButterfly prediction
      res2[i] <- median(pr_inat) ## median iNat prediction
      res3[i] <- median(pr_inat - pr_ebd) ## median difference
    }
    
    
    saveRDS(list(
      ebud = Xp_ebud, inat = Xp_inat, 
      ebud_coef = coef(this_fit_ebud),
      inat_coef = coef(this_fit_inat),
      species = this_species, name_clean = this_name_clean
    ), paste0("intermediate/lpmtx/lp_matrix_", this_name_clean, ".RDS"))
    
    
    target_df$ebud_medPred <- median(res1)
    target_df$inat_medPred <- median(res2)
    target_df$difference <- median(res3)
    
    target_df$ebud_medResp <- nimble::expit(median(res1))
    target_df$inat_medResp <- nimble::expit(median(res2))
    
    target_df$ebud_SE <- sd(res1)
    target_df$inat_SE <- sd(res2)
    target_df$diff_SE <- sd(res3)
    target_df$region <- this_region
    
    if (write_inds) write_csv(target_df, this_output_file)
    target_df
    
  }, timeout = timeout)
}


#### Helper phylogenetic functions ####

dropFamily <- function(x){
  
  
  sep_val <- unlist(strsplit(x, " "))
  
  if(length(sep_val) == 3){
    toReturn <- paste(sep_val[2]," ", sep_val[3], sep="")
    
  }else if(length(sep_val) == 4){
    
    toReturn <- paste(sep_val[3]," ", sep_val[4], sep="")
    
    
  }else if(length(sep_val) == 5){
    toReturn <- paste(sep_val[4]," ", sep_val[5], sep="")
    
  }else{
    toReturn = NA
  }
  
  
  return(toReturn)
  
  
}

getFamily <- function(x){
  sep_val <- unlist(strsplit(x, " "))
  
  
  
  toReturn <- sep_val[1]
  
  return(toReturn)
}


trimTree_A <- function(tree, toUse){
  
  tree1Names = gsub("_", " ", tree$tip.label)
  
  tree1Names_b <- lapply(tree1Names, dropFamily) %>% unlist()
  tree1Names_fam <- lapply(tree1Names, getFamily) %>% unlist()
  
  keep = which(tree1Names_b %in% intersect(tree1Names_b, toUse$species) )
  remove = setdiff(1:length(tree1Names_b), keep)
  toDrop = tree$tip.label[remove]
  
  trimmed_tree <- drop.tip(tree, toDrop)
  
  phylo_cov <- vcv(trimmed_tree)
  
  return(phylo_cov)
  
}


getReadyMeta_A <- function(phylocov, toUse){
  
  dimName_og <- gsub("_", " ", attr(phylocov, which = "dimnames")[[1]]) 
  dimName <- lapply(dimName_og, dropFamily) %>% unlist()
  family <- lapply(dimName_og, getFamily) %>% unlist()
  
  toM  = cbind.data.frame(id = 1:nrow(phylocov), dimName = dimName, family = family) ## 84
  
  
  testM <- merge(toM, toUse, by.x = "dimName", by.y = "species", all.x = T, all.y = F)
  ## inner join
  
  
  testM2 <- testM %>% group_by(id) %>% slice_sample(n=1) %>% arrange(id)
  
  attr(phylocov, which = "dimnames")[[1]] <- testM2$id
  attr(phylocov, which = "dimnames")[[2]] <- testM2$id
  
  
  testM2$target_species2 <- as.factor(testM2$id)
  
  return(list(data = testM2, matrix = phylocov)) 
  
}
