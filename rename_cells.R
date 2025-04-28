suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(stringr))
suppressWarnings(library(tidyverse))
suppressWarnings(library("optparse"))
})



 
option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices.", metavar="character"),
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, eg after removal of ambient RNA", metavar="character"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")
); 

parser=OptionParser(option_list=option_list) 
opt=parse_args(parser)

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

 
cd=getwd()
output=paste0('/', opt$out)
ifelse(!dir.exists(paste0(output)), dir.create(paste0(output)), FALSE)
ifelse(!dir.exists(paste0(output, '/objects')), dir.create(paste0(output, '/objects')), FALSE)



gex_list=list()

if(is.na(opt$matrix)){gex_list=readRDS(opt$rds)}

if(!is.na(opt$matrix)){
for(file in list.files(path)) {
        print(file)
        name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
        gex_list[[name]]<- Read10X_h5(paste0(path,'/',file))
				}
        gex_list= lapply(X=gex_list, FUN=function(x){
        x=CreateSeuratObject(counts=x, assay="RNA", min.cells=10) })
			}



##RENAME CELLS
rename_list=list()
for(i in names(gex_list)){
  current=gex_list[[i]]
	print(i)
  current@meta.data['FullInfo']<- i
  current@meta.data['Sample']<- strsplit(i,'_')[[1]][2]
  current@meta.data['Treatment']<- strsplit(i,'_')[[1]][1]
 
  if(strsplit(i,'_')[[1]][1] =='pro' ){
      current@meta.data['Months']<- 'zero'}

  if(strsplit(i,'_')[[1]][1] =='post' ){
	if(strsplit(i,'_')[[1]][3] == '6m'){
        current@meta.data['Months']<- 'six'}
  
        if(strsplit(i,'_')[[1]][3] == '12m'){
        current@meta.data['Months']<- 'twelve'}
}


#  rownames(current@meta.data)<- paste(i,rownames(current@meta.data), sep="_")
#  current<- RenameCells(current, new.names=rownames(current@meta.data))
current<- RenameCells(current,new.names=paste(i,rownames(current@meta.data), sep="_"))
	

#print(head(current@meta.data))
  rename_list[[i]]<- current
			  }

if(length(names(rename_list)) > 1){
rename_list[['merged']]<- merge(rename_list[[1]], y=c(rename_list[-1]))

}


saveRDS(rename_list, paste0(output, '/objects/renamed_objects.rds.gz'), compress = "gzip")
