

#function to return the nearest gene to a genomic region (hg19)
annotateGenes <- function(df){
  
  # Dependencies
  require(annotatr)
  require(EnsDb.Hsapiens.v75)
  
  # Get biomart annotation
  ens <- genes(EnsDb.Hsapiens.v75, filter=GeneBiotypeFilter("protein_coding"))
  seqlevels(ens) <- paste0("chr",seqlevels(ens))
  
  # Make granges
  gr_query <- GRanges(df)
  gr_subj <- GRanges(ens)
  
  # Get distance to nearest genes
  nearestDist <- distanceToNearest(gr_query,gr_subj)
  nearestDist <- cbind(data.frame(gr_query[queryHits(nearestDist),])[,c("seqnames","start","end")], 
                       SYMBOL=gr_subj[subjectHits(nearestDist),]$symbol,
                       DISTANCE=data.frame(nearestDist)$distance)
  colnames(nearestDist)[[1]] <- "chr"
  
  # Return
  return(nearestDist)
}

AllContrastDMRFF <- function(fit, m, all_contrasts, annotation, max_gapbase=500, se_stat_origin="limma", est="coefficients", stderr=NA, pval="p.value", extra="", dir="./", reGenerateFile=TRUE) {
  
  # Dependencies
  require(rlist)
  require(dmrff)
  
  # Define save filename
  filename <- paste0(dir,"dmrff_results.maxgap",max_gapbase,extra,".tsv")
  
  # If file doesn't exist, run
  if (!file.exists(filename) | reGenerateFile) {
    
    # Initialise output list
    out_dmrs <- list()
    
    # Loop over contrasts
    print("###### RUNNING DE NOVO........")
    for (contrast in all_contrasts) {
      
      # Print
      print(paste0("##### Working with contrast :- ",contrast))
      
      # Define SE (mostly if not limma)
      if (se_stat_origin == "limma") { 
        
        se_vals <- sqrt(fit$s2.post) * fit$stdev.unscaled[,contrast]
        
        # Else if not using limma as output
      } else {
        
        # Use custom-provided SEs
        se_vals <- se_stat_origin
      } 
      
      # Subset to get contrast-specific results
      stats <- data.frame(estimate=fit[[est]][,contrast],
                          se=se_vals,
                          p.value=fit[[pval]][,contrast])
      
      # Get common DE results, methylation matrix and annotation data
      common <- intersect(rownames(m), rownames(annotation))
      annotation <- annotation[match(common, rownames(annotation)),]
      stats <- stats[match(common, rownames(stats)),]
      m <- m[match(common, rownames(m)),]
      
      # Annotate DE results
      stats <- cbind(stats, annotation)
      
      # Apply Dmrff
      dmrs <- dmrff(estimate=stats$estimate,
                    se=stats$se,
                    p.value=stats$p.value,
                    methylation=m,
                    chr=stats$chr,
                    pos=stats$pos,
                    maxgap=max_gapbase, # Distance between CpGs in bp
                    verbose=T)
      
      # Add contrast and code for the site
      dmrs["CONTRAST"] <- contrast
      dmrs["SITE"] <- paste0(dmrs$chr, ":", dmrs$start, "-", dmrs$end)
      
      # Get closest gene and distance
      geneDist <- annotateGenes(dmrs)
      dmrs <- merge(dmrs, geneDist, by=c("chr","start","end"))
      
      # Append
      out_dmrs <- rlist::list.append(out_dmrs, dmrs)
      names(out_dmrs)[[length(out_dmrs)]] <- contrast
    }
    
    # Combine
    out_dmrs.df <- do.call("rbind",out_dmrs)
    
    # Save
    write.table(out_dmrs.df, filename, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Else Read in
  } else {
    
    # read in
    print("###### READING IN FROM FILE........")
    out_dmrs.df <- read.table(filename, sep="\t", header=TRUE)
    
  }
  
  # Return
  return(out_dmrs.df)
  
}

# Function to loop over max gap parameters
DmrffTuneMaxGap <- function(max_gap_vals=c(500,1000,2000), ...) {
  
  # Loop over max gap
  out_df <- do.call("rbind",lapply(max_gap_vals, function(max_gap) {
    
    # Run DMR
    print(paste0("##### RUNNING DMRFF WITH MAX GAP OF = ",max_gap))
    dmrff_df <- data.frame(AllContrastDMRFF(max_gapbase=max_gap, ...), MAXGAP=max_gap)
    
    # Return
    return(dmrff_df)
    
  }))
  
  # Return
  return(out_df)
  
}


AddDMRCpGDetails <- function(all_dmr_df, sig_cpgs, annots) {
  
  # Loop over DMR
  all_dmr_df["CpG"] <- ""
  all_dmr_df["nSigCpg"] <- ""
  all_dmr_df["nNonSigCpg"] <- ""
  all_dmr_df["pSigCpg"] <- ""
  all_dmr_df["n"] <- NA
  
  # Order
  all_dmr_df <- all_dmr_df[order(as.numeric(gsub("chr","",all_dmr_df$chr)), all_dmr_df$start),]
  
  # Set rownames
  rownames(all_dmr_df) <- paste(all_dmr_df$chr, all_dmr_df$start, all_dmr_df$end, sep="_")
  
  # Loop over chrms (to speed things up)
  for (chr in unique(all_dmr_df$chr)) {
    
    # Subset by chrm
    print(chr)
    chr_annots <- annots[annots$chr == chr,]
    chr_dmr_df <- all_dmr_df[all_dmr_df$chr == chr,]
    
    # Loop over DMRs
    dmrs <- rownames(chr_dmr_df)
    for (dmrid in dmrs) {
      
      # Subset to get DMR coordinates
      dmr_df <- chr_dmr_df[dmrid,]
      
      # Get relevant CpGs
      cpgs <- chr_annots[chr_annots$pos >= dmr_df$start &
                         chr_annots$pos <= dmr_df$end,]$Name
      
      # Check CpGs match - NOTE THAT N ONLY REFERS TO CPGS IN THE REGION, NOT
      # THE ACTUAL CPGS IN THE DMR
      # if (length(cpgs) != dmr_df$n) { stop("We're missing some CpGs!") }
      
      # Enumerate significant Cpgs
      sig_dmr_cpgs <- intersect(cpgs, sig_cpgs$Name)
      nonsig_dmr_cpgs <- setdiff(cpgs, sig_cpgs$Name)
      
      # Add as columns
      all_dmr_df[dmrid,]["nSigCpg"] <- length(sig_dmr_cpgs)
      all_dmr_df[dmrid,]["nNonSigCpg"] <- length(nonsig_dmr_cpgs)
      all_dmr_df[dmrid,]["pSigCpg"] <- length(sig_dmr_cpgs) / length(cpgs)
      
      # Add CpGs
      all_dmr_df[dmrid,]["CpG"] <- paste(cpgs, collapse=", ")
      
      # Update n as wrong
      all_dmr_df[dmrid,]["n"] <- length(cpgs)
      
    }
  
  }
    
  # Force numeric
  all_dmr_df$nSigCpg <- as.numeric(all_dmr_df$nSigCpg)
  all_dmr_df$nNonSigCpg <- as.numeric(all_dmr_df$nNonSigCpg)
  all_dmr_df$pSigCpg <- as.numeric(all_dmr_df$pSigCpg)
  
  # Return
  return(all_dmr_df)
  
}



GetTall_DMR_DF <- function(df) {
  
  # Sort
  df <- df[order(as.numeric(gsub(".*_","",df$DMR_ID)), decreasing=FALSE),]
  
  # Initialise input
  tallDMR_df <- data.frame()
  
  # Loop over DMRs
  for (dmrid in df$DMR_ID) {
    
    n_dmr <- which(df$DMR_ID %in% dmrid)
    if (n_dmr %% 500 == 0) { print(dmrid) }
    
    # Subset
    sub_df <- df[df$DMR_ID == dmrid,]
    
    # Get CpGs
    cpgs <- strsplit(sub_df$CpG, split=", ")[[1]]
    
    # Expand CpG out and form output
    ## Get warning about rownames, can ignore
    cpgwise_df <- suppressWarnings(cbind(sub_df[,!grepl("^CpG$",colnames(sub_df))], CpG=cpgs))
    
    # Combine
    tallDMR_df <- rbind(tallDMR_df, cpgwise_df)
  }
  
  # Return
  return(tallDMR_df)
  
}


GetMeanCpgPerDMR <- function(all_dmr_df, mvals) {
  
  # Initialise
  agg_mvals <- data.frame()
  
  # Loop over DMRs
  for (dmrid in unique(all_dmr_df$DMR_ID)) {
    
    # Subset to get DMR coordinates
    dmr_df <- all_dmr_df[all_dmr_df$DMR_ID == dmrid,]
    
    # Get CpGs
    cpgs <- strsplit(dmr_df$CpG, split=", ")[[1]]
    
    # Get mean methylation in DMR
    mean_mvals <- aggregate(. ~ ID, data=data.frame(ID=dmrid, mvals[cpgs,]), FUN=mean)
    agg_mvals <- rbind(mean_mvals, agg_mvals)
    
  }
  
  # Return
  return(agg_mvals)
  
}


RunMethylEnrichment <- function(df, sig_cpgs, groupCol, cpg_df=NULL, cpg_level=TRUE, groupIsBackGroundSubset=TRUE, pval_col="Bonf_padj", ont_db="GO", dir="./", extra="", CpGCol="Name", reGenerateFile=TRUE, ...) {
  
  # Dependencies
  require(missMethyl)
  
  # Define file
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  filename <- paste0(dir, ont_db, ".methyl_enrich", extra, ".tsv")
  
  # If file doesn't exist, run
  if (!file.exists(filename) | reGenerateFile) {
    
    print("##### RUNNING !!!!")
    
    # Loop over contrasts and fit
    out_df <- do.call("rbind",lapply(names(sig_cpgs), function(contr) {
      
      print(contr)
      
      # Subset to get background genes for said contrast
      contr_df <- df[df[[groupCol]] == contr,]
      
      # If we have no currently set background dataframe
      if (is.null(cpg_df)) { 
        
        # Just use all unique features in the test set. NOTE this may not be
        # appropriate if this is already significant filtered
        background <- unqiue(contr_df[[CpGCol]])
      
      # Else if we have a background dataframe supplied and the levels of contrast
      # are the same as the grouping factor for groupCol, then filter this way
      } else if (groupIsBackGroundSubset & !is.null(cpg_df)) { 
        
        # Filter background by grouping column
        background <- cpg_df[cpg_df[[groupCol]] == contr,][[CpGCol]] 
      
      # Else if we have a background dataframe but the levels for contrast for it 
      # aren't the same as the grouping factor for groupCol, take background as-is,
      # deduplicating the CpGs.
      } else {
        
        # Fitler background by grouping column
        background <- cpg_df[[CpGCol]]
        
      } 
      
      # Keep only unique
      background <- unique(background)
      
      # Print length of background set
      print(paste0("##### Testing with number of background CpGs = ", length(background)))
      
      # Toggle whether CpG...
      if (cpg_level) {
        
        # Run enrichment analysis
        gst <- gometh(sig.cpg=contr_df[contr_df[[pval_col]] < sig_cpgs[[contr]],][[CpGCol]], 
                      all.cpg=background, 
                      plot.bias=TRUE, ...)
        
        # ... or DMR
      } else {
        
        # Filter sig
        contr_df <- contr_df[contr_df[[pval_col]] < sig_cpgs[[contr]],]
        
        # Get as GRanges object
        gr <- makeGRangesFromDataFrame(contr_df, seqnames.field="chr", start.field="start", end.field="end")
        
        # Save plot
        png(paste0(dir, gsub("[:]","-",contr), ".", ont_db, ".methyl_enrich", extra, ".bias_plt.png"), width=1500, height=1500, units="px")
        
        # Run enrichment analysis
        gst <- goregion(regions=gr, 
                        all.cpg=background, 
                        plot.bias=TRUE, 
                        collection=ont_db,
                        sig.genes = TRUE,
                        ...)
        
        # Close device
        dev.off()
        
      }
      
      # Add contrast
      gst_df <- cbind(GROUP=contr, gst)
      
    }))
    
    # Get fold enrichment
    out_df["FOLD_ENRICHMENT"] <- out_df$DE / out_df$N
    
    # Save
    print("### Saving.........")
    print(filename)
    write.table(out_df, filename, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Else read in
  } else {
    
    # Read in
    print("##### READING IN FROM FILE!!!!")
    out_df <- read.table(filename, header=TRUE, sep="\t", quote="\"")
    
  }
  
  # Return
  return(out_df)
  
}


PlotTopEnrichTerms <- function(df, sortCol, topN=10, outDir="./", extraDetails="_go", ...) {
  
  # If any sig terms
  if (nrow(df) > 0) {
    
    # Order
    df <- df[order(df[[sortCol]], decreasing=ifelse(sortCol == "FDR", FALSE, TRUE)),]
    
    # Subset, accounting for when less than topN sig terms
    if (nrow(df) >= topN) { df <- df[1:topN,] }
    
    # Plot
    plt <- ggplot(df, aes(x=reorder(TERM,-log10(FDR)), y=-log10(FDR), fill=FOLD_ENRICHMENT, size=log2(N))) + 
      geom_point(stat="identity", shape=21) + 
      scale_fill_continuous(high="red",low="blue") +
      xlab("") +
      theme_bw(base_size=16) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                                     legend.position="top",
                                     axis.title.y=element_text(size=24),
                                     legend.text=element_text(size=6),
                                     legend.box="vertical", legend.margin=margin())
    
    # Save
    dir.create(outDir, recursive=TRUE)
    ggsave(paste0(outDir, "dotPlot_top",topN,"Terms_by",sortCol,extraDetails,"Results.pdf"), plt, units="px", device="pdf", dpi=300, ...)
    ggsave(paste0(outDir, "dotPlot_top",topN,"Terms_by",sortCol,extraDetails,"Results.png"), plt, units="px", device="png", dpi=300, ...)
    
    # View plot
    return(plt)
    
  }
}


