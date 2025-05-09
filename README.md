# README

This is a collection of command line R scripts for analyzing spatial
transcriptomics data.  It is based on Seurat 3.2 workflows with a focus on
multi-sample analyses (technical replicates and treatment/control pairs).

The main purpose of this effort is to implement best practices that can be
launched in automated pipelines. It further provides converter functions that
make it easier to use methods implemented in different workspaces.

## Installation

These scripts require Seurat 3.2 or later with Visium support. SpaceRanger 2.0
or later output requires Seurat 4.2 or later. SpaceRanger 3.0 currently
requires the `visium-hd` branch from seurat and seurat-object GitHub
repositories.

### Dependencies

If you are a conda user, you can find all dependencies available as
conda packages in the conda_environment.yml file.

Otherwise install Seurat directly from CRAN:

```
install.packages("Seurat")
```

A few optional additional packages extending functionality:

```
remotes::install_github("satijalab/seurat-wrappers")
remotes::install_github("navinlabcode/CellTrek")
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

# following packages not necessary with conda_environment.yml
BiocManager::install(c("batchelor",
    "harmony",
    "NMF",
    "corrplot",
    "optparse",
    "pheatmap",
    "patchwork"))
```
 
For the [scvi-tools](https://github.com/scverse/scvi-tools) wrapper, we
recommend using our conda environment and additionally installing the following
packages (using 'mamba' instead of 'conda' should be a major speedup):

```
conda install pytorch torchvision torchaudio "pytorch-cuda>=11.7" -c pytorch -c nvidia
conda install "jaxlib=*=*cuda*" jax cuda-nvcc -c conda-forge -c nvidia
conda install scvi-tools -c conda-forge
```
For [cell2location](https://github.com/BayraktarLab/cell2location), install scvi-tools as above
and then install it via pip in the activate conda environment:

```
pip install cell2location[tutorials]
```

For [Giotto](https://github.com/drieslab/Giotto), optionally if you use conda,
install a few missing dependencies:

```
conda install r-terra r-checkmate r-pak r-rfast leidenalg python-louvain -c conda-forge
```

Then install it via pak in R:

```
pak::pkg_install("drieslab/Giotto")
```

### sttkit

Common functionality in this toolkit are provided in an R package called
`sttkit`.  Install it from GitHub:

```
remotes::install_github('lima1/sttkit')
```

Start R and enter the following to get the path to the command line
scripts:

```
system.file("extdata", package = "sttkit")
``` 

Exit R and store this path in an environment variable, for example in 
BASH:

```
export STTKIT="/path/to/sttkit/extdata"
Rscript $STTKIT/st_normalize.R --help
Usage: /path/to/sttkit/inst/extdata/st_normalize.R [options] ...
```

## Tools

### st_normalize.R

A simple script that takes data from the ST pipeline and uses several Seurat
features to normalize counts and to generate QC plots.

10X Visium SpaceRanger example:
```
#  spaceranger_dir is the path to filtered_feature_bc_matrix.h5
Rscript $STTKIT/st_normalize.R --spaceranger_dir $SAMPLE/outs \
     --sampleid $SAMPLE \ 
     --outprefix OUTDIR/${PIPELINE}/$SAMPLE/normalize/$SAMPLE 
```

We can also additional gene annotation if needed:

```
#  spaceranger_dir is the path to filtered_feature_bc_matrix.h5
Rscript $STTKIT/st_normalize.R --spaceranger_dir $SAMPLE/outs \
     --sampleid $SAMPLE \ 
     --outprefix OUTDIR/${PIPELINE}/$SAMPLE/normalize/$SAMPLE \
     --gtf $REFERENCES/cellranger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
     --spaceranger_probe_set $REFERENCES/spaceranger-2.1.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
```

The GTF extracts a stable gene id such as ENSEMBL that can be used instead of
the gene name. For FFPE Visium, we can provide probe information. Both can be
accessed in R via the `[[]]` operator on the `Spatial` assay. In case
SpaceRanger was run without probe filtering, features without `included` flag
are ignored in the `SCT` assay.

```
head(ndata$Spatial[[]],2)
               gene_id symbol
OR4F5  ENSG00000186092  OR4F5
SAMD11 ENSG00000187634 SAMD11
                                                                                                                                                     probe_seqs
OR4F5  AATGGAGAAAGCCAATTCCCCATGTGACAGCCATAATGCCGACACATGCG
SAMD11 TGTCACCATCGCTGGCAGAGAAGCTGGGAGTTCGCTCCTTCTTCAGGTTC|[...]
         all.included included                       regions
OR4F5           FALSE    FALSE                     unspliced
SAMD11 TRUE|TRUE|TRUE     TRUE unspliced|unspliced|unspliced
```

SpatialTranscriptomics example:
```
PIPELINE="standard"
Rscript $STTKIT/st_snormalize.R --infile $SAMPLE/${PIPELINE}_pipeline/${SAMPLE}_${PIPELINE}_ensembl_adjusted.tsv \
     --sampleid $SAMPLE \ 
     --hejpeg {SAMPLE}_HE_bw_scaled.jpg \
     --outprefix OUTDIR/${PIPELINE}/$SAMPLE/normalize/$SAMPLE \
```


This script will generate a few files with `--outprefix` as filename prefix.
Note that `--hejpeg` is ignored for Visium data. 

![ex_sagittal_1_he_counts](https://user-images.githubusercontent.com/364466/75631970-879b2c80-5bc5-11ea-8e45-87f488357210.png)

### st_cluster.R

Cluster the SpatialTranscriptomics data generated by `st_normalize.R`

Simple single-sample example:
```
Rscript $STTKIT/st_cluster.R \
     --infile $OUTDIR/$SAMPLE/normalize/serialize/${SAMPLE}_scaled.rds \
     --outprefix OUTDIR/$SAMPLE/cluster/$SAMPLE 
```
![ex_sagittal_1_he_cluster](https://user-images.githubusercontent.com/364466/75631996-be714280-5bc5-11ea-88d8-2028430f7316.png)

Advanced multi-sample example with NMF clustering:
```
#! /bin/bash
#$ -S /bin/bash
#$ -pe orte 20 # number of parallel jobs
#$ -N Wtester
#$ -cwd
#$ -j y
#$ -o  ${OUTDIR}/${NORMALIZATION_METHOD}/${SAMPLE}/cluster
#$ -l h_rt=345600      #this need to adapted to your needs
#$ -l m_mem_free=4G    #this need to adapted to your needs
mpirun --mca mpi_warn_on_fork 0 -v -np \$NSLOTS  R --slave \
    -f $STTKIT/st_cluster.R --args \
    --infile lists/${SAMPLE}_${NORMALIZATION_METHOD}_spatial.list \
    --labels lists/${SAMPLE}_labels.list \
    --outprefix $OUTDIR/${NORMALIZATION_METHOD}/$SAMPLE/cluster/$SAMPLE \
    --gmt ../../signatures/all.gmt \
    --extra_gmt ../../signatures/pathways_kegg.gmt \
    --min_features $MIN_FEATURES \
    --nmf --nmf_rank 4:16 --nmf_nruns \$NSLOTS $NMF_RANDOMIZE --nmf_method nsNMF \
    --verbose --mpi
```

Since NMF clustering is slow, we may need to use the doMPI package to run in in
parallel (provide the `--mpi` flag).  This example shows that for multi-sample,
we provide `--infile` a text file with suffix .list containing multiple input
files.  For the non-default `seurat2` or `scran` normalizations, use the unscaled
`${SAMPLE}_unscaled.rds` files to normalize all samples jointly before clustering.

We provide a `--gmt` file with signatures of interest to make sure that the
corresponding genes are not filtered out for lower variance than other genes.
Gene signatures in `--extra_gmt` are not forced to be included and instead
broadly tested against NMF cluster markers. This is useful for providing large
pathway databases such as KEGG or REACTOME.

We use the nsNMF algorithm instead of the default to get a more sparse solution
at the cost of a significantly longer runtime.

Here the results of NMF clustering on the 10X mouse brain example data:

![all_he_nmf_cluster_9_ex_sagittal_1_small](https://user-images.githubusercontent.com/364466/75268883-1410ae00-57c6-11ea-9adf-00bef6b05fef.png)
![all_he_nmf_cluster_9_ex_sagittal_a1_small](https://user-images.githubusercontent.com/364466/75268894-17a43500-57c6-11ea-9dcd-ea44cb1a85ad.png)


Note that the cluster ids are consistent across sections.
  
### st_score.R

Takes the output of `st_cluster.R` and gene signatures in GMT format as input
and plots signature scores (when a clustered RDS was provided, additional
signature per cluster plots will be generated)

Example:

```
Rscript $STTKIT/st_score.R \
    --infile $OUTDIR/$SAMPLE/normalize/serialize/${SAMPLE}_scaled.rds \
    --gmt mm10_io_sigs.gmt \
    --outprefix OUTDIR/$SAMPLE/signatures/$SAMPLE  
```
Advanced feature: `--infile` can be again a list of input files (see
`st_cluster.R`).  In this case violin plots are generated to compare the
signatures across samples.

### st_benchmark.R

Compare Spatial data with bulk RNA-seq.

Example:

```
Rscript $STTKIT/st_benchmark.R \
    --infile $OUTDIR/$SAMPLE/cluster/serialize/${SAMPLE}.rds \
    --htseq ${SAMPLE_BULK}.gene_counts.cts \
    --outprefix OUTDIR/$SAMPLE/benchmark/$SAMPLE  

```
Advanced feature: both `--infile` and `--htseq` can be again a list of input 
files (see `st_cluster.R`).

### st_integrate.R

Integrates SpatialTranscriptomics with a (matched) scRNA reference. The default
is simply following the Seurat best practices as outlined in their [Spatial Vignette](https://satijalab.org/seurat/articles/spatial_vignette.html):

```
Rscript $STTKIT/st_integrate.R \
    --infile $OUTDIR/$SAMPLE/cluster/serialize/${SAMPLE}.rds \
    --outprefix $OUTDIR/$SAMPLE/integrate/$SAMPLE \
    --singlecell allen_cortex.rds \
    --labels_singlecell allen_cortex \
    --refdata subclass
```

Here, the reference scRNA-seq dataset is expected to be normalized by
`sctransform` and contains cell type annotation in a `type` meta data column
(the column can be changed with `--refdata` as in this example).  Again,
`--singlecell` can be a list of reference datasets. Specify
`--integration_method rctd` to use [RCDT](https://github.com/dmcable/spacexr),
`--integration_method scvi_destvi` to use [DestVI](https://github.com/scverse/scvi-tools),
`--integration_method scvi_cell2location` to use [cell2location](https://github.com/BayraktarLab/cell2location),
or `--integration_method giotto` for [SpatialDWLS from
Giotto](https://github.com/drieslab/Giotto) instead. Output files and plots are
equivalent.


![ex_sagittal_1_he_labels_allen_cortex_1_small](https://user-images.githubusercontent.com/364466/75380489-21e93080-58a5-11ea-8d1a-75950b0dd104.png)
![ex_sagittal_a1_he_labels_allen_cortex_1_small](https://user-images.githubusercontent.com/364466/75380495-2281c700-58a5-11ea-97d7-efa00e79914e.png)

All celltype predictions can be easily loaded in Seurat and compared:

```
ls $OUTDIR/$SAMPLE/integrate/serialize/*transfer*
LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_giotto_transfer_predictions.rds
LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_rctd_multi_transfer_predictions.rds
LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_scvi_destvi_transfer_predictions.rds
LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_scvi_cell2location_transfer_predictions.rds
LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_seurat_transfer_predictions.rds
```

In R:

```
x <- readRDS("cluster/serialize/LIB-021633rd1.rds")
x$predictions <- readRDS("integrate/serialize/LIB-021633rd1_742abcb4d6052d8416d7d7a47d0f6749_rctd_transfer_predictions.rds")
```

This can now be used following the Seurat best practices.

We also provide a convenient way of averaging prediction in a simple consensus
method:

```
files <- dir("integrate/serialize", pattern = "transfer_predictions.rds",
    full.names = TRUE)
tp_consensus <- find_assayobject_consensus(lapply(files, function(x)
    readRDS(x)[[1]]), labels = labels, plot_correlations = run_plots,
        plot_cor_method = "kendall")
x$predictions <- tp_consensus
```

We also support the [CellTrek](https://github.com/navinlabcode/CellTrek)
package that performs coembedding of the single-cell and spatial data to
generate the training model. The single cells are then charted on to their
spatial locations using non-linear interpolation to augment the ST spots.
This method works especially well when matched single cell and spatial data
are available.

We have adapted the same to work using command line inside of sttkit, and also
splitting the various cell-types on to separate panels as shown below

```
 Rscript $STTKIT/st_integrate.R \
    --infile $OUTDIR/$SAMPLE/cluster/serialize/${SAMPLE}.rds \
    --outprefix $OUTDIR/$SAMPLE/celltrek/$SAMPLE \
    --singlecell allen_cortex.rds \
    --labels_singlecell allen_cortex \
    --refdata subclass --png --serialize \
    --integration_method celltrek
```
![brain_he_celltrek_dots_labels_brain](https://user-images.githubusercontent.com/405520/181096917-260e847f-074d-41ee-a225-54dab2121184.png)

![ex_sagittal_1_he_celltrek_dots_labels_allen_cortex](https://user-images.githubusercontent.com/405520/181086643-de237893-006d-4a27-b41b-2d5c7c72ba64.png)

The `CellTrek::celltrek_vis` function uses RShiny to visualize all cell-types
in the mouse brain sample.  We can easily load the `celltrek` object in R:

```
cd $OUTDIR/$SAMPLE/celltrek
R
```
```
library(CellTrek)
library(dplyr)
options("browser" = "google-chrome")
# The serialized RDS object is a list for cases when multiple
# single cell references were provided
brain_celltrek <- readRDS("serialize/ex_sagittal_a2_celltrek.rds")[[1]]
brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=sort(unique(brain_celltrek$cell_type)))
CellTrek::celltrek_vis(brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
                       brain_celltrek@images$ex_sagittal_a2@image, brain_celltrek@images$ex_sagittal_a2@scale.factors$lowres)
```

Now choose `cell_type` under "Color" and then click "Plot".

![celltrek_example_output](https://user-images.githubusercontent.com/364466/181936155-51c0f0c1-878b-4d3e-97d9-53b61f487f6d.png)

### st_enhance.R

Imputes data from neighboring spots. Currently only
[BayesSpace](https://github.com/edward130603/BayesSpace) supported.

```
Rscript $STTKIT/st_enhance.R \
   --infile $OUTDIR/$SAMPLE/cluster/serialize/${SAMPLE}.rds \
   --outprefix $OUTDIR/$SAMPLE/enhance/$SAMPLE \
```

When the provided `--infile` contains SCTransform normalized data, it will
use those log counts. Otherwise BayesSpace's own normalization is used.

![mouse_10x_bayesspace](https://user-images.githubusercontent.com/364466/122808326-c19cb000-d29a-11eb-95b4-da14b0f6d4f0.png)

### st_hejpeg.R

Some standard edits to H&E jpegs (obsolete with Visium).

Example:
```
Rscript $STTKIT/st_hejpeg.R  --infile LP_L10012_S085_TGFB_EX2_LIB-026528rd1.jpg \
    --outfile LIB-026528rd1_HE_bw_scaled.jpg --dither
```



## Example workflow

In the following we run `sttkit` on 5 technical replicates of a breast cancer
sample.

```
PROJECT="/mnt/tmplabdata/ngdx/projects/dev/spatialTranscriptomics/NGDX-P00273"
PIPELINE="standard"
NORMALIZATION="sctransform"
OUTDIR="../../data/$NORMALIZATION"
MIN_FEATURES=400   # exclude spots with fewer than 400 detected genes
NUM_FEATURES=3000  # aim for including ~3000 genes
MIN_SPOTS=1        # include genes detected in a single spot

mkdir -p $OUTDIR/$PIPELINE

SAMPLES=("LIB-021633rd1" "LIB-021634rd1" "LIB-021635rd1" "LIB-021636rd1" "LIB-021637rd1" )
for SAMPLE in "${SAMPLES[@]}"
do
    rm -rf $OUTDIR/$PIPELINE/$SAMPLE
    SLIDE="ST_LP_L4_009_02JUN2018_Breast_EX2"

    Rscript ~/git/CancerGenetics/ncgs-in-spatial_tools/sttkit/inst/extdata/st_normalize.R \
        --infile $PROJECT/$SLIDE/$SAMPLE/${PIPELINE}_pipeline/${SAMPLE}_ensembl_adjusted.tsv \
        --outprefix $OUTDIR/${PIPELINE}/$SAMPLE/normalize/$SAMPLE \
        --sampleid $SAMPLE \
        --hejpeg $PROJECT/$SLIDE/Images/${SAMPLE}_HE_bw_scaled.jpg \
        --min_features $MIN_FEATURES \
        --min_spots $MIN_SPOTS \
        --num_features $NUM_FEATURES \
        --normalization_method $NORMALIZATION

    Rscript ~/git/CancerGenetics/ncgs-in-spatial_tools/sttkit/inst/extdata/st_cluster.R \
        --infile $OUTDIR/$PIPELINE/$SAMPLE/normalize/serialize/${SAMPLE}_scaled.rds \
        --outprefix $OUTDIR/$PIPELINE/$SAMPLE/cluster/$SAMPLE 

done

# We can specify groups of samples and cluster them together.
#
# In this example, we use all high quality samples and name the group
# "all_good" (I'm very good at naming things...)
#
SAMPLES=("all_good")
for SAMPLE in "${SAMPLES[@]}"
do
    rm -rf $OUTDIR/$PIPELINE/$SAMPLE

echo "#! /bin/bash
#$ -S /bin/bash
#$ -pe orte 20 # number or parallel jobs
#$ -N Wtester
#$ -cwd
#$ -j y
#$ -o  $OUTDIR/$PIPELINE/${SAMPLE}/cluster
#$ -l h_rt=345600      
#$ -l m_mem_free=8G    
mpirun --mca mpi_warn_on_fork 0 -v -np \$NSLOTS  R --slave \
    -f ~/git/CancerGenetics/ncgs-in-spatial_tools/sttkit/inst/extdata/st_cluster.R \
    --args --infile lists/${SAMPLE}_${NORMALIZATION}_spatial.list \
    --outprefix $OUTDIR/$PIPELINE/${SAMPLE}/cluster/${SAMPLE} \
    --nmf --nmf_ranks 2:12 --nmf_randomize --nmf_method nsNMF \
    --png --force --mpi --verbose

" > ${OUTDIR}/$PIPELINE/$SAMPLE/${SAMPLE}_cluster.sh

qsub ${OUTDIR}/$PIPELINE/$SAMPLE/${SAMPLE}_cluster.sh

done
```

The .list file simply list input files line by line:

```
cat all_good_sctransform_spatial.list
../../data/sctransform/standard/LIB-021633rd1/normalize/serialize/LIB-021633rd1_scaled.rds
../../data/sctransform/standard/LIB-021634rd1/normalize/serialize/LIB-021634rd1_scaled.rds
../../data/sctransform/standard/LIB-021635rd1/normalize/serialize/LIB-021635rd1_scaled.rds
../../data/sctransform/standard/LIB-021636rd1/normalize/serialize/LIB-021636rd1_scaled.rds
../../data/sctransform/standard/LIB-021637rd1/normalize/serialize/LIB-021637rd1_scaled.rds
```
