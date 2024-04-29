

NiceVolcanoPlots <- function(res, comparison, dir, pval_colname="padj", logfc_colname="log2FoldChange", feature_colname="ENSG", 
                             padj_thr=0.01, topn_genes=15, label_sig_thr=1e-45, n.max_ovlp=25, logfc_thr=0.25,
                             fills=c("lightgrey", "#BB5566", "#004488"), colours=c("white","white", "white"), dot_transparency=0.75, dot_shape=1,
                             extra_filename="") {
  
  # Packages
  require(ggplot2)
  require(ggrepel)
  
  # Loop over results df
  plts <- lapply(unique(res[[comparison]]), function(contrast) {
    
    # Print current comparisons
    print(contrast)
    
    # Subset by contrast
    res_df <- res[res[[comparison]] == contrast,]
    
    # Remove NA
    res_df <- res_df[!is.na(res_df[[pval_colname]]),]
    
    # Get sig
    sig.res_df <- res_df[res_df[[pval_colname]] < padj_thr & !is.na(res_df[[pval_colname]]),]
    
    # If any significant results
    n_sig <- nrow(res_df[which(res_df[[feature_colname]] %in% sig.res_df[[feature_colname]]),])
    
    # Add sig key
    res_df["Sig"] <- "Non-Sig"
    if (n_sig > 0) { 
      
      res_df[which(res_df[[feature_colname]] %in% sig.res_df[[feature_colname]]),]["Sig"] <- paste0("FDR<",padj_thr)
      
    }
    
    # Add lfc key
    res_df["Direction"] <- "Up"
    res_df[res_df[[logfc_colname]] < 0,]["Direction"] <- "Down"
    
    # Add colour key
    res_df["ColorKey"] <- paste(res_df$Direction, res_df$Sig, sep=" ")
    res_df[grepl("Non-Sig",res_df$ColorKey),]["ColorKey"] <- "Non-Sig"
    
    # Order plotting
    res_df$ColorKey <- factor(res_df$ColorKey , levels=c("Non-Sig", paste0("Down FDR<",padj_thr), paste0("Up FDR<",padj_thr)))
    
    # Order by LFCs (for labelling)
    res_df <- res_df[order(res_df[[logfc_colname]], decreasing=TRUE),]
    
    # Define genes to label
    sig.res_df <- res_df[res_df[[pval_colname]] < padj_thr & !is.na(res_df[[pval_colname]]),]
    
    
    # Define colours
    names(colours) <- c("Non-Sig", paste0("Up FDR<",padj_thr), paste0("Down FDR<",padj_thr))
    names(fills) <- c("Non-Sig", paste0("Up FDR<",padj_thr), paste0("Down FDR<",padj_thr))
    
    # Plot
    plt <- ggplot(res_df, aes(x=.data[[logfc_colname]], y=-log10(.data[[pval_colname]]), fill=ColorKey, color=ColorKey)) + 
      geom_vline(xintercept=0, linetype="dashed", color="darkgrey") +
      geom_vline(xintercept=logfc_thr, linetype="dashed", color="darkgrey") +
      geom_vline(xintercept=-logfc_thr, linetype="dashed", color="darkgrey") +
      geom_point(alpha=dot_transparency, shape=dot_shape, size=3) + 
      scale_fill_manual(values=fills) +
      scale_color_manual(values=colours) +
      geom_hline(yintercept=-log10(padj_thr), linetype="dashed", color="darkgrey") +
      ylab("-log10 FDR") + xlab("Log2 FC") + ggtitle(paste0(comparison, ": ", mgsub::mgsub(contrast, c("-Day_","_"), c(": "," ")))) +
      theme_bw(base_size=16) + theme(legend.title=element_blank(), legend.position="top", legend.direction="horizontal", 
                                     axis.title.x = element_text(size=32, face="bold"), axis.title.y = element_text(size=32, face="bold"),
                                     plot.title = element_text(size=40), legend.text = element_text(size=20)) +
      guides(color = guide_legend(override.aes = list(size=10), nrow=1, byrow=TRUE), text=FALSE)
    
    # If any labels
    res_df["SHAPE"] <- "Circle"
    if (topn_genes > 0) {
      
      # Get labels of interest
      label_df <- rbind(head(sig.res_df, topn_genes), tail(sig.res_df, topn_genes), sig.res_df[sig.res_df[[pval_colname]] < label_sig_thr,])
      label_df <- label_df[!duplicated(label_df),]
      
      # Add size key
      res_df[which(res_df[[feature_colname]] %in% label_df[[feature_colname]]),]["SHAPE"] <- "Triangle"
      
      # Add gene annotation
      plt <- plt + geom_text_repel(data=label_df, mapping=aes(x=.data[[logfc_colname]], y=-log10(.data[[pval_colname]]), label=.data[[feature_colname]]), size=3, show.legend = FALSE, max.overlaps=n.max_ovlp)
      
    }
    
    # I
    
    # Save
    # setwd(dir)
    ggsave(gsub(" ","_",paste0(dir,comparison,"-",contrast, ".volcano_plot.padj_thr",padj_thr,".",extra_filename,".png")), plt, units="px", width=4000, height=3500)
  })
  
  # Return
  return(plts)
}