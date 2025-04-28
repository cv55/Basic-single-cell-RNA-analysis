suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(stringr))
suppressWarnings(library(SoupX))
suppressWarnings(library(DropletUtils))
suppressWarnings(library("optparse"))
suppressWarnings(library(scDblFinder))
})

option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices.", metavar="character"),

make_option(c("--raw-matrix"), type="character", default=NA,dest='raw',
              help="Path with input raw matrices.", metavar="character"),

make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")

);


parser=OptionParser(option_list=option_list)
opt=parse_args(parser)


NM='vst'
DR='UMAP'

path=(opt$matrix)
if(str_sub(path,-1,-1) =='/'){path=substr(path,1,nchar(path)-1)}

rawpath=(opt$raw)
if(str_sub(rawpath,-1,-1) =='/'){rawpath=substr(rawpath,1,nchar(rawpath)-1)}


#output=paste0('/', opt$out)
output=opt$out
cd=getwd()
dir.create(paste0(output), showWarnings = FALSE)
dir.create(paste0(output, '/objects'), showWarnings = FALSE)

gex_list=list()

runSoupx_h5<- function(filt, rawh5, resolution, name, norm, DR){

    filt.matrix=Read10X_h5(filt,use.names = T)
    raw.matrix=Read10X_h5(rawh5,use.names = T)
    soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
    srat=CreateSeuratObject(counts = filt.matrix, min.cells=50)

        if(norm=='SCT'){srat<-SCTransform(srat, verbose = F)}

        else{
                srat <- NormalizeData(srat)
                srat <- FindVariableFeatures(srat, selection.method = "vst",nfeatures = 2000, verbose=FALSE)
                srat <- ScaleData(srat)
            }

            srat    <- RunPCA(srat, verbose = F)
            if(DR=='UMAP'){srat    <- RunUMAP(srat, dims = 1:30, verbose = F)}
            if(DR=='TSNE'){srat    <- RunTSNE(srat, dims = 1:30, verbose = F,check_duplicates = FALSE)}
            srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
            srat    <- FindClusters(srat, verbose = T, resolution = resolution)

    soup.channel  <- setClusters(soup.channel, setNames(srat@meta.data$seurat_clusters, rownames(srat@meta.data)))
    if(DR=='UMAP'){soup.channel  <- setDR(soup.channel, srat@reductions$umap@cell.embeddings)}
    if(DR=='TSNE'){soup.channel  <- setDR(soup.channel, srat@reductions$tsne@cell.embeddings)}
    soup.channel  <- autoEstCont(soup.channel,forceAccept=TRUE)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

    return(adj.matrix)

    }

 cd=getwd()

print('REMOVE AMBIENT')
  dir.create(paste0(output,'/SoupX'), showWarnings = TRUE)

   n=0
   for(file in list.files(path, pattern='h5')){
     n=n+1
     name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
#	name=strsplit(file, '\\.')[[1]][1]
     print(paste0("Removing ambient RNA from ", name))

     rawp=paste0(rawpath,'/',name,'_raw.h5')
     filteredp=paste0(path,'/',name,'_filtered.h5')

	
     adj.matrix <-runSoupx_h5(filteredp, rawp, 0.5, name,NM,DR)
#      DropletUtils:::write10xCounts(paste0(output,'/SoupX/',name), adj.matrix)
      fobj=CreateSeuratObject(counts=adj.matrix, assay="RNA", min.cells=10)
      gex_list[[name]]=fobj
    }
dev.off()

saveRDS(gex_list, paste0(output, '/objects/preprocessed_list_SoupX.rds.gz'), compress = "gzip")


