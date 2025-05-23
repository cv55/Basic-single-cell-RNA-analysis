suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
suppressWarnings(library(cowplot))
suppressWarnings(library("optparse"))
library(glmGamPoi)
library(sctransform)
})




option_list = list(
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, after QC", metavar="character"),
make_option(c("--norm-method"), type="character", default="vst",dest='norm',
              help="Normalization method, either SCT or vst. Default is LogNorm (vst)", metavar="character"),
make_option(c("--regress-out"), type="character", dest='regress', default=NA,
              help="Values to regress out", metavar="vector"),
make_option(c("--exclude"), type="character", dest='exclude', default=NA,
              help="Multiomes to exclude", metavar="vector"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")
);

parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

if(!is.na(opt$out)){
  outpath=(opt$out)
  if(str_sub(outpath,-1,-1) =='/'){
  outpath=substr(outpath,1,nchar(outpath)-1)
  }}


output=paste0('/', opt$out)
ifelse(!dir.exists(output), dir.create(output), FALSE)
ifelse(!dir.exists(paste0(output, '/objects')), dir.create(paste0(output, '/objects')), FALSE)

gex_norm_list=c()
gex_qc_list=readRDS(opt$rds)


#Set regression values
if(!is.na(opt$regress)){regress= as.vector(strsplit(opt$regress, ',')[[1]])}
if(is.na(opt$regress)){
        opt$regress=c()
        regress=opt$regress}
print(paste0("vars to regress: ", regress))


#Set exclude values
if(!is.na(opt$exclude)){exclude= as.vector(strsplit(opt$exclude, ',')[[1]])}
if(is.na(opt$exclude)){
        opt$exclude=c()
        exclude=opt$exclude}

print(paste0("Experiments to exclude: ", exclude))


seurat_norm<- function(sobject, method=NULL, regress.out=c()){
  DefaultAssay(sobject) <- "RNA"

  if(method=="SCT"){
    print("SCT is applied")
    sobject=SCTransform(sobject, vars.to.regress = regress.out, verbose=FALSE )
    return(sobject)
  }

  else if(method=="vst"){
    sobject <- NormalizeData(sobject)
    sobject <- FindVariableFeatures(sobject, selection.method = "vst",nfeatures = 2000, verbose=FALSE)
    sobject <- ScaleData(sobject, vars.to.regress = regress.out)
    return(sobject)
  }
}

gex_qc_list<- gex_qc_list[which(!names(gex_qc_list) %in% exclude)]

for (i in names(gex_qc_list)) {
    if(i!="merged"){
	    gex_norm_list[[i]] <- seurat_norm(gex_qc_list[[i]],method=opt$norm, regress.out=regress)
    }}

gex_norm_list_clean<-gex_norm_list
gex_norm_list_clean<- gex_norm_list_clean[which(!names(gex_norm_list_clean) %in% exclude)]
saveRDS(gex_norm_list_clean, paste0(output, '/objects/gex_norm_list_clean.rds.gz'), compress = "gzip")

print('Here0')
# Find most variable features across samples to integrate
integ_features_clean <- SelectIntegrationFeatures(object.list = gex_norm_list_clean, nfeatures = 8000) 
print('Here00')

# Merge normalized samples
merged_seurat_norm <- merge(x = gex_norm_list_clean[[1]],
		       y = gex_norm_list_clean[2:length(gex_norm_list_clean)],
		       merge.data = TRUE)

print('Here')

if(opt$norm=="SCT"){DefaultAssay(merged_seurat_norm) <- "SCT"}
print('Here1')

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat_norm) <- integ_features_clean
print('Here2')

# Calculate PCs using manually set variable features
if(opt$norm=="SCT"){merged_seurat_norm <- RunPCA(merged_seurat_norm, assay = "SCT", npcs = 70)}

#saveRDS(gex_norm_list, paste0(cd, output, '/objects/Normalized_list.rds.gz'), compress = "gzip")

saveRDS(merged_seurat_norm, paste0(output, '/objects/Normalized_mergedObject.rds.gz'), compress = "gzip")

#saveRDS(combined.sct, paste0(cd, output, '/objects/Seurat_integration.rds.gz'), compress = "gzip")
