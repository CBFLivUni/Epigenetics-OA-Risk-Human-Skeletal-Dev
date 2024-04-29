

GetProxySnps <- function(snps_df, lds, genome_build, measure, pop="EUR", tkn="d355f7f8ae2f") {
  
  # Dependencies
  require(LDlinkR)
  
  # Loop over chromosomes
  snp_ld <- lapply(unique(snps_df$chr), function(chrm) {
    
    # Record which chromosomes
    print(paste0("## Currently on :- ", chrm))
    
    # Get SNPs on chromosome
    chrm_snps <- unique(snps_df[snps_df$chr == chrm,]$Variant)
    
    # Get LDs between SNPs
    ld_m <- LDmatrix(c(chrm_snps, unique(lds[lds$chr == chrm,]$PROXY)), 
                     pop = pop, 
                     r2d = measure, 
                     token = tkn, 
                     file = FALSE,
                     genome_build = genome_build
    )
    
    # Return
    return(ld_m)
  }); names(snp_ld) <- unique(snps_df$chr)
  
  # Return
  return(snp_ld)
  
}


AssignProxiesByMeasure <- function(df, seen=c(), measure="R2") {
  
  # Order by LD
  df <- df[df$FILTER == measure,]
  df <- df[order(df[[measure]], decreasing=TRUE),]
  
  # Initialise output
  out_df <- data.frame()
  missing_snps <- c()
  
  # Loop over SNPS
  for (snp in unique(df$OA_SNP)) { 
    
    # Subset to get top SNP-wise proxies
    sub_df <- df[df$OA_SNP == snp,]
    
    # Remove previously seen proxies
    sub_df <- sub_df[!sub_df$PROXY %in% seen,]
    
    # If there is still an available proxy SNP, add this
    if (nrow(sub_df) > 0) {
      
      # Keep top and add to output
      out_df <- rbind(out_df, sub_df[1,])
      
      # Update seen proxies
      seen <- c(seen, sub_df[1,]$PROXY)
      
    }
  }
  
  # Return
  return(out_df)
  
}

KeepBestProxy <- function(snps, lds, common, nonmulti, multi_prox, r2, dp, measure, ld_filter=0, dprime_filter=1) {
  
  # Initialise final output and missing SNPs output of this code
  final_out.df <- data.frame()
  collapsed_snps.df <- data.frame()
  
  # Loop over multi proxies
  for (proxy in names(multi_proxies)) { 
    
    # Subset to get proxty-wise SNPs and SNPs that share the proxy (alongside their respective possible proxies)
    print(paste0("## ON PROXY :- ",proxy))
    proxy_df <- lds[lds$PROXY == proxy,]
    snp_df <- data.frame(lds[which(lds$OA_SNP %in% proxy_df$OA_SNP),])
    
    # Order by R2 then D'
    # snp_df <- snp_df
    
    # Loop over measurements and filter (R2 first, then D prime)
    measure_df <- data.frame()
    for (measure in sort(unique(snp_df$FILTER), decreasing=TRUE)) {
      
      # Assign top proxies
      out <- AssignProxiesByMeasure(snp_df, measure=measure)
      measure_df <- rbind(out, measure_df)
      
    }
    
    # Remove duplicate (ie if both R2 and D prime, remove the D prime)
    measure_df <- measure_df[!duplicated(measure_df$OA_SNP),]
    final_out.df <- rbind(final_out.df, measure_df)
    
    # Get missing SNPs
    missing_snps <- unique(snp_df$OA_SNP)[!unique(snp_df$OA_SNP) %in% final_out.df$OA_SNP]
    
    # If there are any missing SNPs (ie those dropped due to sharing a proxy with another SNP that had greater LD for said proxy)
    if (length(missing_snps) > 0) {
      
      # Loop over missing SNPs
      for (snp in missing_snps) {
        
        # Get missing SNP chromosome
        missing_snp_chr <- unique(snps[snps$Variant == snp,]$chr)
        
        # Get chromosome-wise SNPs
        rownames(r2[[missing_snp_chr]]) <- r2[[missing_snp_chr]]$RS_number
        rownames(dp[[missing_snp_chr]]) <- dp[[missing_snp_chr]]$RS_number
        chrm_snp.r2 <- r2[[missing_snp_chr]][,-1]
        chrm_snp.dp <- dp[[missing_snp_chr]][,-1]
        
        # Get possible SNPs for the current missing SNP
        snp_r2s <- chrm_snp.r2[!rownames(chrm_snp.r2) %in% snp, snp, drop=F]
        snp_dps <- chrm_snp.dp[!rownames(chrm_snp.dp) %in% snp, snp, drop=F]
        
        # Remove SNPs not on the array and get max SNP by R2
        snp_r2s <- snp_r2s[drop=F,which(rownames(snp_r2s) %in% unique(lds$PROXY)),]
        snp_dps <- snp_dps[drop=F,which(rownames(snp_dps) %in% unique(lds$PROXY)),]
        
        next_snp.r2 <- rownames(snp_r2s)[which.max(snp_r2s[[snp]])]
        next_snp.dp <- rownames(snp_dps)[which.max(snp_dps[[snp]])]
        
        # If the R^2 between SNPs is above threshold
        if (snp_r2s[which.max(snp_r2s[[snp]]),] > ld_filter) {
          
          # Go with the next SNP by R^2
          next_snp <- next_snp.r2
          
          # Print out collapsed SNPs
          print(paste0("### Collapsing ",snp, " to ",next_snp, " R2 = ",snp_r2s[which.max(snp_r2s[[snp]]),]))
          
          # Add to collapsed SNPs df
          collapsed_snps.df <- rbind(collapsed_snps.df, 
                                     data.frame(SNP=snp, 
                                                COLLAPSED_TO=next_snp, 
                                                CHRM=missing_snp_chr, 
                                                MEASURE="R2", 
                                                R2=snp_r2s[which.max(snp_r2s[[snp]]),], 
                                                D=snp_dps[which.max(snp_dps[[snp]]),]))
          
          # Else
        } else {
          
          # Go with next SNP by D'
          next_snp <- next_snp.dp
          
          # Print out collapsed SNPs
          print(paste0("### Collapsing ",snp, " to ",next_snp, " D' = ",snp_dps[which.max(snp_dps[[snp]]),]))
          
          # Add to collapsed SNPs df
          collapsed_snps.df <- rbind(collapsed_snps.df, 
                                     data.frame(SNP=snp, 
                                                COLLAPSED_TO=next_snp, 
                                                CHRM=missing_snp_chr, 
                                                MEASURE="Dprime", 
                                                R2=snp_r2s[which.max(snp_r2s[[snp]]),], 
                                                D=snp_dps[which.max(snp_dps[[snp]]),]))
          
        }
      }
    }
  }
  
  # Combine non-duplicated proxies to processed proxies
  final_out.df <- rbind(nonmulti, final_out.df[,-ncol(final_out.df)])
  
  # Add SNPs on microarray
  final_out.df <- rbind(common, final_out.df)
  
  # Deduplicate (this preferentially retains SNPs on the microarray)
  final_out.df <- final_out.df[!duplicated(final_out.df$OA_SNP),]
  
  # Sanity check LD and D prime are above thresholds
  final_out.df[final_out.df$R2 > ld_filter | final_out.df$Dprime > dprime_filter, ]
  
  # Filter collapsed by LD
  collapsed_snps.df[collapsed_snps.df$R2 > ld_filter | collapsed_snps.df$D > dprime_filter,]
  
  # Add sig
  final_out.df["SNP_TYPE"] <- "Proxy required"
  final_out.df[final_out.df$OA_SNP == final_out.df$PROXY,]["SNP_TYPE"] <- "On Array"
  
  # Return
  return(final_out.df)
  
}

