# scSNViz: visualization and analysis of Cell-Specific Expressed SNVs
scSNViz is a specialized tool for the visualization and analysis of single-cell expressed SNVs from cell-barcoded scRNA-seq data. scSNViz enables quantitative assessment of SNV expression, 2D and 3D visualization of individual or user-defined groups of variants, expression-based clustering of SNVs, and cross-sample comparisons. Beyond visualization, scSNViz enables the estimation, summarization, and graphical representation of SNV expression metrics, facilitating the study of allelic dynamics across somatic, germline, and RNA-originating SNVs. To support integrative transcriptomic analyses, scSNViz interoperates with established frameworks including Seurat for clustering, Slingshot for trajectory inference, scType for cell type annotation, and CopyKat for copy number profiling.

 
<img src='https://github.com/HorvathLab/scSNViz/blob/5b9c6326effd5bf75d08b17ef54abb2c393bf8aa/docs/scSNViz_PanelA.png' width=50% height=50%>


## Installation

Ubuntu, MacOS and Windows are currently supported.

If using Debian OS, you may need to install the following libraries:
sudo apt install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libcurl4-openssl-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev libcairo2

Note: If using macOS or a Linux-based distribution, the plot_snv_data() and individual_snv_plots() functions will prompt the user to select the number of cores for parallel processing. This can drastically improve the runtime for large samples. The number of cores used is at the user’s discretion. For certain versions of macOS, this is limited to two cores, if you encounter issues using more than two, this may be the cause. If you are unsure how many cores are available, the R parallel package (loaded with scsnviz) provides the detectCores() function to check.

```
# Enter commands in R (or R studio)
install.packages('devtools')
library(devtools)

install_github('navinlabcode/copykat')

if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('glmGamPoi')
BiocManager::install('SingleCellExperiment')
BiocManager::install('slingshot')

install_github('HorvathLab/scSNViz')
library(scsnviz)
```
If the above fails due to download timeout, try to increase the global options timeout. E.g. options(timeout=3600)

If the above fails due to rate limits, try generating a GitHub Personal Access Token (PAT), add it into your environment and then run again.

If the above gives errors related to randomcoloR, it may be due to the V8 package dependency. Please refer to the V8 README https://github.com/jeroen/V8.

## Load Libraries

```
packages <- c('copykat', 'HGNChelper', 'parallel', 'readr', 'Seurat', 'SingleCellExperiment', 'slingshot')
lapply(packages, library, character.only=TRUE)
```

## Workflow for a single sample

#### Define paths to input files, and define the output directory.

The input files are located in the input folder on github. The snv file is an output from SCReadCounts. The user may provide a .tsv file that is not from SCReadCounts as long as it is also a .tsv and contains the following columns: CHROM, POS, REF, ALT, ReadGroup, SNVcount, and RefCount.
```
snv_file <- 'input/sample1_SNVs.tsv'
srt_obj_file <- 'input/sample1_Seurat_object.rds'
output_dir = 'output'    # or output directory of your choice
```

#### Read in the counts matrix (from either an .RDS file of an existing Seurat object or a counts matrix)
```
#gene.matrix <- Read10X(data.dir=countsmatrix_file) # for reading in a countsmatrix, the data.dir may also be the directory for that contains barcodes.tsv, genes.tsv and matrix.mtx, such as: /user/filtered_gene_bc_matrices/hg19/
#sample1 <- CreateSeuratObject(counts=gene.matrix, min.cells=3, min.features=200, project='Sample1')
sample1 <- readRDS(srt_obj_file)
sample1@project.name = 'Sample1' # set the project name
sample1$orig.ident = 'Sample1' # set the project name
```

#### Quality Control: Filter data and perform scaling and normalization
```
# define the percentage of counts per cell that originate from mitochondrial genes 
sample1[['percent.mt']] <- PercentageFeatureSet(sample1, pattern='^MT-') # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'

# plot the number of features (or genes) per cell, the number of counts per cell and the percent.mt per cell
VlnPlot(sample1, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3) # As described in Seurat introductory Vignettes

# filter the Seurat object based on the violin plot
sample1 <- subset(sample1, subset= nFeature_RNA>1000 & nFeature_RNA<7500 & nCount_RNA<50000 & percent.mt<15) # Modify numbers appropriate to your violin plot
```


#### Scale and normalize the data. Then, run a PCA.
```
sample1 <- SCTransform(object=sample1, vst.flavor='v2', method='glmGamPoi',
           vars.to.regress='percent.mt', verbose = F)
sample1 <- RunPCA(sample1)
sample1 <- FindNeighbors(sample1, dims=1:10)
sample1 <- FindClusters(sample1, resolution=0.5)

```

#### Preprocess the SNV data and incorporate the Seurat object into the workflow
```
# preprocessing the SNV data
processed_data <- preprocess_snv_data(rds_obj=sample1,
                                      snv_file=snv_file,
                                      dimensionality_reduction='UMAP',
                                      th_vars=0,
                                      th_reads=0,
                                      enable_sctype = TRUE, #to classify cell types using sctype
                                      tissue_type='Immunesystem', # other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics=TRUE,
                                      output_dir=output_dir)
```

#### Generate 3d dimensionality reduction plots
```
plots <- plot_snv_data(seurat_object=processed_data$SeuratObject,
                       processed_data$ProcessedSNV,
                       processed_data$AggregatedSNV,
                       processed_data$PlotData,
                       output_dir=output_dir,
                       include_histograms=TRUE,  
                       dimensionality_reduction='UMAP',
                       include_cell_types=TRUE,
                       include_copykat=FALSE, # CNV metrics produced by copykat; this may significantly increase processing time depending on the size of gene counts matrix provided
                       include_snv_dim_red=FALSE, # If set to TRUE, this function transposes the SNVxBarcode matrix and generates a dimensionality reduction plot to view similarity between SNVs.
                       slingshot=TRUE,
                       color_scale='YlOrRd',
                       cell_border=0,
                       save_each_plot=TRUE)
```

<img src='https://github.com/HorvathLab/scSNViz/blob/0fc797784b905b29c744968114a5e2192de213b2/docs/sample_outputs.png'>

#### Generate individual SNV plots
```
#Individual SNV's plottable capped at 50 unique.
ind_snv_plots <- individual_snv_plots(seurat_object=processed_data$SeuratObject,
                                      processed_snv=processed_data$ProcessedSNV,
                                      sig_snvs=processed_data$SigSNV,
                                      output_dir=output_dir,
                                      slingshot=TRUE,
                                      save_each_plot=TRUE,
                                      dimensionality_reduction='UMAP',
                                      dynamic_cell_size=FALSE)
```

<img src='https://github.com/HorvathLab/scSNViz/blob/0fc797784b905b29c744968114a5e2192de213b2/docs/individual_snv_plots.png'>

#### Plot individual SNV
This function allows you to pass either a single SNV of interest as a text string. Or a list of SNVs of interest. Lists should contain SNV of interest in double quotes.
```
one_snv_plot <- single_snv_plot(
       seurat_object=processed_data$SeuratObject,
       processed_snv=processed_data$ProcessedSNV,
       snv_of_choice='1:155169447:C:T',
       snv_list_of_choice = NULL,
       output_dir='output/1_155169447_C_T',
       slingshot=TRUE,
       dimensionality_reduction='UMAP',
       dynamic_cell_size=FALSE,
       save_each_plot=TRUE
     )

```

#### Plot all SNVs for a single gene
If using labeled data with gene symbol annotations in the 'GENE' column of the SNV (e.g. sample1_SNVs_large.tsv has gene annotations in the column 'GENE'), the following function may be run to plot the N_REF, N_VAR, and VAF for all SNVs contained within a single cell.
```
one_gene_plot <- single_gene_plot(
    seurat_object=processed_data$SeuratObject,
    processed_snv=processed_data$ProcessedSNV,
    gene_of_choice=gene,
    output_dir=paste0('output/', 'S1PR1'),
    slingshot=TRUE,
    dimensionality_reduction='UMAP',
    dynamic_cell_size=FALSE,
    save_each_plot=TRUE
    )
```

#### Generate exploratory combined plots report
```
generate_report(plot_object=plots,
                ind_snv_object=ind_snv_plots,
                hide_ind_plots=TRUE, # Set this to FALSE in order to see plots for each individual SNV.
                output_dir=output_dir)
```


<img src='https://github.com/HorvathLab/scSNViz/blob/dev/docs/Exploratory_combined_plots.png'>


#### Generate exploratory combined plot for single SNV of interest
```
generate_report(plot_object=plots,
                ind_snv_object=one_snv_plot,
                hide_ind_plots=FALSE,
                output_dir=output_dir)
```
#### Generate Transposed SNV plot (with or without labels)

In order to generate a transposed SNV plot, you must have a minimum of 100 unique SNVs in your SNV file. To run the tutorial with sample data, please re-run the above workflow for an individual sample, but using the 'input/sample1_SNVs_large.tsv' file. Once you have generated the processed_data variable, plot_snv_data may be run with the include_snv_dim_red variable set to TRUE.


<img src='https://github.com/HorvathLab/scSNViz/blob/d4bf67d0d2a33897c9e317e84c6e3301dc1a6c03/docs/transposed_snv_plot.png'>

## Workflow for Multiple Samples
The following is a workflow that calculates and overlays basic SNV metrics on top of a dimensionality reduction integrated from multiple samples.

#### Prepare integrated data
The below workflow is shown for unprocessed Seurat objects. Ideally this workflow is started with Read10X()/CreateSeuratObject() or a Seurat object that has not been processed. If the Seurat object has been processed, already, this may result in errors.

```
output_dir = "output_integrated_samples"    # or output directory of your choice

sample1 <- readRDS('input/sample1_Seurat_object.rds')
sample2 <- readRDS('input/sample2_Seurat_object.rds')
sample1$orig.ident = 'sample1'
sample2$orig.ident = 'sample2'

## QC ##
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes
sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot


sample2[["percent.mt"]] <- PercentageFeatureSet(sample2, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(sample2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sample2 <- subset(sample2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot
########

srt_merged <- merge(x=sample1, y=c(sample2), add.cell.ids=c('sample1', 'sample2'))
srt_merged <- SCTransform(srt_merged, vars.to.regress='percent.mt', verbose=F)
srt_merged <- RunPCA(srt_merged)
srt_integrated <- IntegrateLayers(object=srt_merged, method=CCAIntegration, normalization.method='SCT', new.reduction='integrated', verbose=F) #any of the suggested integration methods in Seurat may be applied here for the methods parameter
srt_integrated <- FindNeighbors(srt_integrated, reduction='integrated', dims=1:30)
srt_integrated <- FindClusters(srt_integrated, resolution=0.6)
```

#### Prepare SNV data
The sample labelling must be consistent with the orig.ident set in the step that prepares the integrated data.

```
snv_sample1 <- read.table('input/sample1_SNVs.tsv', sep='\t', header=T)
snv_sample2 <- read.table('input/sample2_SNVs.tsv', sep='\t', header=T)

snv_sample1$ReadGroup = paste0('sample1_', snv_sample1$ReadGroup) # these IDs must match the added cell IDs from above
snv_sample2$ReadGroup = paste0('sample2_', snv_sample2$ReadGroup) # these IDs must match the added cell IDs from above

snv_file <- rbind(snv_sample1,snv_sample2)
write_tsv(snv_file,'snv_file_integrated.tsv')
```


#### Preprocess the SNV data and incorporate the integrated Seurat object into the workflow. Generate plots.
```
processed_data <- preprocess_snv_data(rds_obj=srt_integrated,
                                      snv_file='snv_file_integrated.tsv',
                                      dimensionality_reduction='UMAP', #you may only generate UMAP plots with integrated samples
                                      th_vars=0,
                                      th_reads=0,
                                      enable_integrated=TRUE,
                                      integrated_reduction_name='integrated',
                                      enable_sctype=TRUE, # to classify cell types using scType
                                      tissue_type='Immunesystem', #other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics=TRUE,
                                      output_dir=output_dir)
```

#### Generate 3D dimensionality reduction plots
```
plots <- plot_snv_data(seurat_object=processed_data$SeuratObject,
                       processed_data$ProcessedSNV,
                       processed_data$AggregatedSNV,
                       processed_data$PlotData,
                       output_dir=output_dir,
                       include_histograms=TRUE,  
                       dimensionality_reduction='umap',
                       include_cell_types=TRUE,
                       include_copykat=FALSE, # this is not currently an option for integrated objects
                       slingshot=FALSE, # this is not currently an option for integrated objects
                       color_scale='YlOrRd',
                       cell_border=0,
                       enable_integrated=TRUE,
                       save_each_plot=TRUE)
```

<img src='https://github.com/HorvathLab/scSNViz/blob/0fc797784b905b29c744968114a5e2192de213b2/docs/integrated_plot.png'>


#### Generate individual SNV plots
```
ind_snv_plots <- individual_snv_plots(seurat_object=processed_data$SeuratObject,
                                      processed_snv=processed_data$ProcessedSNV,
                                      sig_snvs=processed_data$SigSNV,
                                      output_dir=output_dir,
                                      slingshot=FALSE,
                                      save_each_plot=TRUE,
                                      dimensionality_reduction='UMAP',
                                      dynamic_cell_size=FALSE,
                                      enable_integrated=TRUE)
```

#### Generate exploratory combined plots report
```
generate_report(plot_object=plots,
                ind_snv_object=ind_snv_plots,
                hide_ind_plots=TRUE, # individual plots for each SNV are hidden
                output_dir=output_dir)
```

<img src='https://github.com/HorvathLab/scSNViz/blob/dev/docs/integrated_output_example.png'>

#### Contact Information
Please contact Siera Martinez (siera.martinez@gwu.edu) or Luke Johnson (luke.johnson@gwu.edu) with any questions.

#### Copyright and licensing
Code copyright 2024 scSNViz https://github.com/HorvathLab/scSNViz/blob/main/LICENSE.md
