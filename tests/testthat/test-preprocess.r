pckages <- c('HGNChelper', 'Seurat')
sapply(pckages, require, character=TRUE)
test_that(desc = 'preprocess_snv_data test', code = {
    here::i_am('sample1_Seurat_object.rds')
    snv_file <- here::here('sample1_SNVs.tsv')
    srt_obj_file <- here::here('sample1_Seurat_object.rds')
    output_dir = 'output_preprocess'

    sample1 <- readRDS(srt_obj_file)
    sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")
    sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA < 50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot

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

    expect_false(is.null(processed_data$ProcessedSNV))
    expect_false(is.null(processed_data$AggregatedSNV))
    expect_false(is.null(processed_data$PlotData))
    expect_false(is.null(processed_data$SeuratObject))
    expect_true(any(grepl('SNVCount',colnames(processed_data$AggregatedSNV))))
    expect_true(any(grepl('ReadGroup',colnames(processed_data$ProcessedSNV))))
    expect_true(any(grepl('HasSNV',colnames(processed_data$PlotData))))
    expect_true(length(list.files('output')) > 0)
    expect_true(class(processed_data$SeuratObject)[1]=="Seurat")
})
