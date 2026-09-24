#' Single Gene, All SNVs Plot
#'
#' This function generates individual SNV plots for VAF, N_VAR, and N_REF. It uses processed SNV data
#' and a Seurat object to create plots with options for saving and including slingshot trajectories.
#'
#' @importFrom Seurat Embeddings as.SingleCellExperiment
#' @importFrom slingshot slingshot slingCurves
#' @importFrom SingleCellExperiment reducedDims reducedDims<-
#' @importFrom plotly plotly_json
#'
#' @param seurat_object Processed Seurat object.
#' @param processed_snv Data frame of processed SNV information, typically output from the preprocess_snv_data function.
#' @param gene_of_choice Selected gene
#' @param output_dir Directory to save plots and HTML (if save_each_plot is TRUE).
#' @param slingshot Logical; whether to include slingshot trajectories. Default: TRUE.
#' @param dimensionality_reduction Dimensionality reduction method ('UMAP', 'PCA', 'tSNE'). Default: "UMAP".
#' @param dynamic_cell_size Logical; whether to scale cell size dynamically based on SNV and reference read counts. Default: FALSE.
#' @param save_each_plot Logical; whether to save each plot individually. Default: FALSE.
#' @param gridlines Logical; whether to include gridlines in the plots. Default: TRUE.
#' @return A list containing JSON content for VAF plots.
#' @details
#' This function generates an individual SNV plot using processed SNV data (processed_snv) and the dimensionality
#' reduction embeddings from a Seurat object, and the SNV of interest in the following format, 1:155169447:C:T. The plots visualize key metrics such as mean VAF (Variant Allele Fraction).
#'
#' The plot can be saved individually in the specified output_dir if save_each_plot is set to TRUE.
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' single_gene_plot <- single_gene_plot(
#'        seurat_object=processed_data$SeuratObject,
#'        processed_snv=processed_data$ProcessedSNV,
#'        gene_of_choice=gene,
#'        output_dir=paste0('output/',gene),
#'        slingshot=TRUE,
#'        dimensionality_reduction='UMAP',
#'        dynamic_cell_size=FALSE,
#'        save_each_plot=TRUE
#'        )
#'        }
#'
#' @export
#'
single_gene_plot <- function(seurat_object, processed_snv, gene_of_choice, output_dir = NULL, slingshot = T,
                                 dimensionality_reduction = "UMAP", dynamic_cell_size = F, save_each_plot = F,
                                 gridlines = T) {

  cat("\nGenerating individual gene SNVs plot...\n")
  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("Invalid dimensionality_reduction method. Please use one of: 'umap', 'pca', 'tsne'.")
  }
  dim.title <- switch(dimensionality_reduction, umap = "UMAP",
                      pca = "PCA", tsne = "tSNE")
  pal <- c("#EBEBEB", "#85C1E9", "#E74C3C", "#B03A2E", "#641E16")
  df.dim <- as.data.frame(Embeddings(seurat_object, reduction = dimensionality_reduction))
  colnames(df.dim) <- c("x", "y", "z")
  df.snv <- processed_snv
  df.snv <- df.snv[c("CHROM", "POS", "REF", "ALT", "ReadGroup", "VAF", 'GENE')]
  gene_options <- paste(df.snv$GENE)
  if (gene_of_choice %in% gene_options) {
  gene_options <- gene_of_choice
  } else {
        stop("gene not present")
    }



  individual_gene_html <- NULL
  curves <- NULL

  if (slingshot) {
    sce <- as.SingleCellExperiment(seurat_object, assay='SCT')
    if (!dimensionality_reduction %in% names(reducedDims(sce))) {
      reducedDims(sce)[[dimensionality_reduction]] <- Embeddings(seurat_object,
                                                                 reduction = tolower(dimensionality_reduction))
    }
    sce <- slingshot(sce, clusterLabels = "seurat_clusters",
                     reducedDim = dimensionality_reduction)
    curves <- slingCurves(sce, as.df = T)
    colnames(curves)[1:3] <- c(paste0(dim.title, "_1"),
                               paste0(dim.title, "_2"),
                               paste0(dim.title, "_3"))
  }


  generate_gene_plots <- function(selected_gene, title_color = "blue", dynamic_cell_size = F) {
    df_subset <- df.snv[df.snv$GENE == selected_gene, ]

    vaf <- df_subset$VAF[match(colnames(seurat_object), df_subset$ReadGroup)]
    readgroups <- df_subset$ReadGroup[match(colnames(seurat_object), df_subset$ReadGroup)]
 
    y <- data.frame(x = df.dim[, 1], y = df.dim[, 2], z = df.dim[, 3],
                    vaf = vaf, ReadGroup = readgroups)
    mean_vaf = aggregate(vaf ~ ReadGroup, data=y, FUN=mean, na.rm=TRUE)
    y <- y[!duplicated(y$ReadGroup),]
    y$vaf <- NULL
    y = merge(mean_vaf, y, by='ReadGroup')
    plots <- list()

    y$vaf_label <- "Undetected"
    y$vaf_label[y$vaf == 0] <- "0 VAF, N_REF Only"
    y$vaf_label[0 < y$vaf & y$vaf <= 0.25] <- "0<VAF<=0.25"
    y$vaf_label[0.25 < y$vaf & y$vaf <= 0.75] <- "0.25<VAF<=0.75"
    y$vaf_label[0.75 < y$vaf & y$vaf <= 1] <- "0.75<VAF<=1.00"
    y$vaf_label <- factor(y$vaf_label, levels = c("Undetected", "0 VAF, N_REF Only",
                                                  "0<VAF<=0.25", "0.25<VAF<=0.75",
                                                  "0.75<VAF<=1.00"))

    # VAF plots
    f_vaf <- plot_ly(type = "scatter3d", mode = "markers+lines")
    for (j in 1:5) {
      cur_label <- levels(y$vaf_label)[j]
      if (dynamic_cell_size) {
        f_vaf <- f_vaf %>%
          add_trace(
            data = subset(y, vaf_label == cur_label), x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "markers", size = 0.05,
            marker = list(color = pal[j], line = list(width = 0)), name = cur_label
            )
      }
      else {
        f_vaf <- f_vaf %>%
          add_trace(
            data = subset(y, vaf_label == cur_label), x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "markers", size = 0.05,
            marker = list(color = pal[j], line = list(width = 0)), name = cur_label
            )
      }
    }
    if (!is.null(curves)) {
      f_vaf <- f_vaf %>% add_trace(data = curves,
                                   x = ~get(paste0(dim.title, "_1")),
                                   y = ~get(paste0(dim.title, "_2")),
                                   z = ~get(paste0(dim.title, "_3")),
                                   split = ~Lineage, mode = "lines", line = list(width = 2))
    }
    f_vaf <- f_vaf %>% layout(title = list(text = "VAF_RNA",
                                           font = list(color = title_color)),
                              scene = list(xaxis = list(title = paste0(dim.title, "_1")),
                                           yaxis = list(title = paste0(dim.title, "_2")),
                                           zaxis = list(title = paste0(dim.title, "_3"))))

    plots[['VAF']] <- f_vaf
    return(plots)
  }


  # function to save plots' json
  plots_json <- lapply(gene_options, function(gene) {
    plots <- generate_gene_plots(gene, title_color = 'blue')
    list(
      VAF = list(
        id = paste0("plot_VAF_", gene),
        json = plotly::plotly_json(plots[['VAF']], jsonedit = F)
      )
    )}
  )

  plots_json <- unlist(plots_json, recursive = F)
  if (save_each_plot && !is.null(output_dir)) {
    save_snv_plot <- function(plot, selected_gene, plot_type) {
      file_path <- file.path(
        output_dir,
        plot_type, paste0(plot_type, "_", gene, ".html")
      )
      dir.create(dirname(file_path), showWarnings = F, recursive = T)

      gene_title_format <- selected_gene
      plot <- plot %>% layout(
        title = list(text = paste(plot_type, "<br>", gene_title_format),
                     font = list(color = "black")), margin = list(t = 50)
      )

      suppressWarnings(saveWidget(as_widget(plot), file = file_path, selfcontained = F, libdir = "lib"))
    }
    for (gene in gene_options) {
      plots <- generate_gene_plots(selected_gene = gene, title_color = "black")
      save_snv_plot(plots[["VAF"]], gene, "VAF")
    }
  }
  ind_snv_out <- list()
  ind_snv_out[["plots_json"]] <- plots_json
  ind_snv_out[["gene_options"]] <- gene_options

  cat("\nIndividual gene SNVs plots saved.\n")

  if (!gridlines){
    
  no_axis <- list(
  showgrid = FALSE,
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showspikes = FALSE,
  title = "",
  backgroundcolor = "rgba(0,0,0,0)",
  showbackground = FALSE)

  for (i in seq_along(ind_snv_out[["plots_json"]])){
  parsed <- jsonlite::fromJSON(ind_snv_out[["plots_json"]][[i]]$json, simplifyVector = FALSE)
  parsed$layout$scene$xaxis <- no_axis
  parsed$layout$scene$yaxis <- no_axis
  parsed$layout$scene$zaxis <- no_axis
  ind_snv_out[["plots_json"]][[i]]$json <- jsonlite::toJSON(parsed, auto_unbox = TRUE)
    }
  }

  return(ind_snv_out)

}