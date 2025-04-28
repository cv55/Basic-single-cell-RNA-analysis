suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(stringr))
suppressWarnings(library("optparse"))
suppressWarnings(library(scDblFinder))
})

option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices.", metavar="character"),
make_option(c("--rds"), type="character", default=NA,dest='rds',
              help="Path with input rds object file (eg after SoupX).", metavar="character"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")
);


parser=OptionParser(option_list=option_list)
opt=parse_args(parser)


print("HERE")
if(!is.na(opt$matrix)){
path=(opt$matrix)
if(str_sub(path,-1,-1) =='/'){
path=substr(path,1,nchar(path)-1)
}}

if(!is.na(opt$out)){
  outpath=(opt$out)
  if(str_sub(outpath,-1,-1) =='/'){
  outpath=substr(outpath,1,nchar(outpath)-1)
  }}


output=paste0('/', opt$out)
ifelse(!dir.exists(paste0(output)), dir.create(paste0(output)), FALSE)
ifelse(!dir.exists(paste0(output, '/objects')), dir.create(paste0(output, '/objects')), FALSE)

gex_list=list()
nodbl_list=list()
if(is.na(opt$matrix)){gex_list=readRDS(opt$rds)}


print('Find DOUBLETS')
min_umi=100

if(!is.na(opt$matrix) & is.na(opt$rds)){
   for(file in list.files(path, pattern='h5')){
     name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
	print(name)
     filteredp=paste0(path,'/',name,'_filtered.h5')
	print(filteredp)
	fobj<- Read10X_h5(filteredp)
     fobj=CreateSeuratObject(counts=fobj, assay="RNA", min.cells=10)
     gex_list[[name]]=fobj
}
}

for(name in names(gex_list)){
	current=gex_list[[name]]
	print(current)
	print(paste0("Removing cells with <0 UMIs from ", name))
	print(paste0("Removing doublets from ", name))

	current=subset(current, subset = nCount_RNA>min_umi)
	print(current)

	sce <- scDblFinder(GetAssayData(current, slot="counts"), clusters=FALSE)
	current$scDblFinder.score <- sce$scDblFinder.score
	current$scDblFinder.class <- sce$scDblFinder.class
	nodbl_list[[name]]<-current
	
}

saveRDS(nodbl_list, paste0(output, '/objects/preprocessed_list_DblFinder.rds.gz'), compress = "gzip")
