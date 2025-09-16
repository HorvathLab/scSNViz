pckages <- c('dplyr', 'ggplot', 'htmlwidgets', 'Matrix', 'parallel', 'plotly', 'Rtsne', 'randomcoloR', 'saveWidget', 'Seurat', 'slingshot')
sapply(pckages, require, character=TRUE)

test_that(desc = 'plot_snv_data test', code = {
    snv_file <- '../../input/sample1_SNVs.tsv'
    srt_obj_file <- '../../input/sample1_Seurat_object.rds'
    output_dir = "output_singlesnvplot"

    snv_file <- '../../input/sample1_SNVs.tsv'
    srt_obj_file <- '../../input/sample1_Seurat_object.rds'

    sample1 <- readRDS(srt_obj_file)
    sample1@project.name = 'Sample1'
    sample1$orig.ident = 'Sample1'

    sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")
    sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot

    sample1 <- SCTransform(object = sample1, vst.flavor = "v2", method = "glmGamPoi",
           vars.to.regress = "percent.mt", verbose = F)
    sample1 <- RunPCA(sample1)
    sample1 <- FindNeighbors(sample1, dims = 1:10)
    sample1 <- FindClusters(sample1, resolution = 0.5)

    processed_data <- preprocess_snv_data(rds_obj = sample1,
                                        snv_file = snv_file,
                                        dimensionality_reduction = "UMAP",
                                        th_vars = 0,
                                        th_reads = 0,
                                        enable_sctype = TRUE,
                                        tissue_type = "Immunesystem",
                                        generate_statistics = TRUE,
                                        output_dir = output_dir)
    plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
                        processed_data$ProcessedSNV,
                        processed_data$AggregatedSNV,
                        processed_data$PlotData,
                        output_dir = output_dir,
                        include_histograms = TRUE,
                        dimensionality_reduction = "umap",
                        include_cell_types = TRUE,
                        include_copykat = FALSE,
                        include_snv_dim_red = FALSE,
                        slingshot = TRUE,
                        color_scale = "YlOrRd",
                        cell_border = 0,
                        save_each_plot = TRUE)
    ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                          processed_snv = processed_data$ProcessedSNV,
                                          sig_snvs = processed_data$SigSNV,
                                          output_dir = output_dir,
                                          slingshot = TRUE,
                                          save_each_plot = TRUE,
                                          dimensionality_reduction = "UMAP",
                                          dynamic_cell_size = FALSE)
    one_snv_plot <- single_snv_plot(
       seurat_object = processed_data$SeuratObject,
       processed_snv = processed_data$ProcessedSNV,
       snv_of_choice = "1:155169447:C:T",
       output_dir = "output/1_155169447_C_T",
       slingshot = TRUE,
       dimensionality_reduction = "UMAP",
       dynamic_cell_size = FALSE,
       save_each_plot = TRUE
     )
    expect_false(is.null(one_snv_plot$plots_json))
    expect_false(is.null(one_snv_plot$snv_options))
})
