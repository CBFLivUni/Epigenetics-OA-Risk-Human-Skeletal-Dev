
RunPCA <- function(m)
{
  # Run PCA
  res <- prcomp(t(m), center=TRUE, scale=TRUE) # Scale and center values
}


# Function for plotting a legend for chosen variable
PlotLegend <- function(data, variable, color_vals=NA, save_dir="./", extra_details="", ncol=2, width=1000, height=250, chosen_unit="px") {
  
  if (all(is.na(color_vals))) { 
    
    # Plot legends
    plt <- ggplot(data, aes(.data[[variable]], fill = .data[[variable]])) + 
      geom_bar() + guides(fill=guide_legend(ncol=ncol))
    
  } else {
    
    # Plot legends
    plt <- ggplot(data, aes(.data[[variable]], fill = .data[[variable]])) + 
      geom_bar() + guides(fill=guide_legend(ncol=ncol)) + scale_fill_manual(values=color_vals)
    
  }
  
  # Make legend
  legend <- cowplot::get_legend(plt)
  
  # Plot
  png(paste0(save_dir,variable,".",extra_details,".legend.png"), width=width, height=height, units=chosen_unit)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
}



Plot2DPCA <- function(meta, var, pcs=c(1,2), pch_size=3, color_by=NA, shape_by=NA, label_id="SAMPLE_ID", default_color="black", reorder_colors=NA)
{
  # Set cell lines and axis labels
  xlabel <- paste("PC",pcs[[1]]," ~ ",round(var$variance.percent[pcs[[1]]],2),"%",sep="")
  ylabel <- paste("PC",pcs[[2]]," ~ ",round(var$variance.percent[pcs[[2]]],2),"%",sep="")
  
  # Set colour factor levels
  if (all(!is.na(reorder_colors))) { meta[[color_by]] <- factor(meta[[color_by]], levels=reorder_colors)} 
  
  # If adding a separate shape
  if (is.na(shape_by))
  {
  
    # Plot PCA plot, coloured by color-by
    pca_plt <- ggplot(meta,
                      aes(x=.data[[paste0("PC",pcs[[1]])]],
                          y=.data[[paste0("PC",pcs[[2]])]], 
                          fill=.data[[color_by]],
                          label=.data[[label_id]],
                          order=.data[[color_by]],
                          color=.data[[color_by]])) + 
      geom_point(size=pch_size, alpha=0.75, color=default_color, shape=21) + xlab(xlabel) + ylab(ylabel) + 
      geom_vline(xintercept=0,linetype="dashed",color="grey") + 
      geom_hline(yintercept=0,linetype="dashed",color="grey") +
      theme_bw(base_size=16) + theme(legend.position = "none")
  
  } else {
    
    # Plot PCA plot, coloured by color-by
    pca_plt <- ggplot(meta,
                      aes(x=.data[[paste0("PC",pcs[[1]])]],
                          y=.data[[paste0("PC",pcs[[2]])]], 
                          shape=.data[[shape_by]],
                          group=.data[[shape_by]],
                          label=.data[[label_id]],
                          fill=.data[[color_by]],
                          order=.data[[color_by]],
                          color=.data[[color_by]])) + 
      geom_point(size=pch_size, alpha=0.75, color=default_color) + xlab(xlabel) + ylab(ylabel) + 
      geom_vline(xintercept=0,linetype="dashed",color="grey") + 
      geom_hline(yintercept=0,linetype="dashed",color="grey") +
      theme_bw(base_size=16) + theme(legend.position = "none") 
    
  }
    
  # Return
  return(pca_plt)
}

CorrPCs <- function(meta, t)
{
  # Correlate PCs 1-4
  cor_df <- data.frame(
    TRAIT=unname(t),
    TRAIT_CLASS=unlist(lapply(names(t), function(name) ifelse(is.null(name), "", name))),
    PC1=unlist(lapply(unname(t), function(col) { print(col); cor(meta$PC1, meta[[col]], method="spearman", use="p") })),
    PC2=unlist(lapply(unname(t), function(col) cor(meta$PC2, meta[[col]], method="spearman", use="p"))),
    PC3=unlist(lapply(unname(t), function(col) cor(meta$PC3, meta[[col]], method="spearman", use="p"))),
    PC4=unlist(lapply(unname(t), function(col) cor(meta$PC4, meta[[col]], method="spearman", use="p")))
  )
  
  # Make tall
  tall.cor_df <- reshape2::melt(cor_df, id.vars=c("TRAIT", "TRAIT_CLASS"))
  
  # Return
  return(tall.cor_df)
}


SubsetRotations <- function(pca, n="ENSG")
{
  # Subset by
  ## If its a character
  if (is.character(n))
  {
    # Subset with grepl
    pca <- as.matrix(pca[!grepl(n, rownames(pca)),])
    
  ## Else subset by indices
  } else {
    # Subset with grepl
    pca <- as.matrix(pca[1:n,])
    
  }
  
  # Return
  return(pca)
}


PlotContribBarPlot <- function(data, annot, colorcol, colors, cols=c("PC1","PC2"), reorder_col=FALSE)
{
  # If multiple cols of interest (need to use rownames)
  if (length(cols) == 1) {
    
    # Get tall format
    tall.data <- reshape2::melt(data$rotation[,cols])
    tall.data <- merge(annot, tall.data, by.x="GENE", by.y="row.names")
    
  } else {
    
    # Get tall format
    tall.data <- reshape2::melt(data$rotation[,cols])
    tall.data <- merge(annot, tall.data, by.x="GENE", by.y="Var1")
    
  }
  
  # Toggle if reorder values by stat
  if (reorder_col == FALSE) {
    
    # Plot
    plt <- ggplot(tall.data, aes(x=GENE, y=value, fill=.data[[colorcol]])) + 
      geom_bar(stat="identity") + 
      scale_fill_manual(values=colors) +
      ylab("Contribution") + xlab("") +
      theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1), legend.title = element_blank(), legend.position = "top") + 
      guides(fill = guide_legend(nrow=2, byrow=TRUE))
    
  } else { 
    
    # Plot
    plt <- ggplot(tall.data, aes(x=reorder(GENE, value), y=value, fill=.data[[colorcol]])) + 
      geom_bar(stat="identity") + 
      scale_fill_manual(values=colors) +
      ylab(paste0(cols," Contribution")) + xlab("") +
      theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1), legend.title = element_blank(), legend.position = "top") + 
      guides(fill = guide_legend(nrow=2, byrow=TRUE))
    
  }
  
  # Return
  return(plt)
}


# Function 
PlotPCABiplot <- function(data, var, annotation, meta, colors, colorcol, point_color_col, point_cols, pcs=c(1,2))
{
  # Create fector rotation df
  rot_df <- data.frame(data$rotation)
  pc_df <- data.frame(data$x)
  
  # Scale
  rot_df <- data.frame(scale(rot_df, scale=TRUE, center=TRUE))
  pc_df <- data.frame(scale(pc_df, scale=TRUE, center=TRUE))
  
  # Add metadata
  pc_df <- merge(meta[,c("SRR",point_color_col),], pc_df, by.x="SRR", by.y="row.names")
  
  # Join annotation with rotation
  rot_df <- merge(annotation, rot_df, by.x="GENE", by.y="row.names")
  
  # Set cell lines and axis labels
  xlabel <- paste("PC",pcs[[1]]," ~ ",round(var$variance.percent[pcs[[1]]],2),"%",sep="")
  ylabel <- paste("PC",pcs[[2]]," ~ ",round(var$variance.percent[pcs[[2]]],2),"%",sep="")
  
  plt <- ggplot() + 
    geom_vline(xintercept=0, linetype="dashed", color="darkgrey") +
    geom_hline(yintercept=0, linetype="dashed", color="darkgrey") +
    geom_circle(rot_df, mapping=aes(x0=0, y0=0, r=max(c(PC1*1.25, PC2*1.25))), linetype="dashed", color="darkgrey") +
    geom_segment(rot_df, mapping=aes(x=0, y=0, xend=PC1, yend=PC2, color=.data[[colorcol]]),
                 arrow = arrow(length = unit(0.25, "cm"))) +
    geom_text(rot_df[abs(rot_df$PC1) > mean(rot_df$PC1) | abs(rot_df$PC2) > mean(rot_df$PC1),], 
              mapping=aes(label=GENE, x=PC1, y=PC2, color=.data[[colorcol]]), size=3, nudge_y=0.00025) +
    geom_point(pc_df, mapping=aes(x=PC1, y=PC2, fill=.data[[point_color_col]]), pch=21, alpha=0.65, size=5) +
    scale_color_manual(values=colors) + 
    scale_fill_manual(values=point_cols) +
    xlab(xlabel) + ylab(ylabel) +
    theme_bw() + theme(legend.title = element_blank(),
                       legend.position = "top") + 
    guides(color = guide_legend(override.aes = list(size = 0.5, shape=3), nrow=2, byrow=TRUE))
  
  # Return
  return(plt)
}

