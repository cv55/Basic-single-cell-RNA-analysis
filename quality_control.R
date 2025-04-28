
suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
suppressWarnings(library(cowplot))
suppressWarnings(library("optparse"))
suppressWarnings(library(dplyr))
})

option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices.", metavar="character"),
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, eg after removal of ambient RNA", metavar="character"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character"),
make_option(c("--use-quantiles"), type="character", default="FALSE",dest='qlogical',
              help="Whether to perform per sample QC using the quantile defined", metavar="character"),
make_option(c("--quantiles"), type="numeric", default=0.9,dest='qthreshold',
              help="When quantiles are used, specify the threshold", metavar="character"),
make_option(c("--performQC"), type="character", default="FALSE",
              help="To either perform QC or not", metavar="character")

);
parser=OptionParser(option_list=option_list)
opt=parse_args(parser)


nFeature_RNA_min=300
nFeature_RNA_max=40000
nCount_RNA_min=600
nCount_RNA_max=15000
percent_mito=15
percent_rb=40
QC=c(nFeature_RNA_min,nFeature_RNA_max,nCount_RNA_min,nCount_RNA_max,percent_mito,percent_rb)

sample_order=c('pro_903_0m','post_903_6m','pro_2553_0m','post_2553_6m',
	       'pro_3377_0m','post_3377_6m','pro_3631_0m','post_3631_6m',
	       'pro_3639_0m','post_3639_6m',
	       'pro_3738_0m','post_3738_6m','post_3738_12m',
	       'pro_2732_0m','post_2732_6m',
	       'pro_3916_0m','post_3916_6m','pro_4054_0m','post_4054_6m',
	       'pro_4156_0m','post_4156_6m','pro_805_0m','post_805_6m',
	       'pro_3835_0m','post_3835_6m','post_3835_12m', 'pro_4153_0m', 'post_4153_6m',
	       'pro_4283_0m','post_4283_6m','pro_3888_0m','post_3888_6m','pro_4375_0m','post_4375_6m',
	       'pro_3890_0m','post_3890_6m','pro_3916BCR_0m','post_3916BCR_6m','pro_4443_0m','post_4443_6m',
	       'pro_4493_0m','post_4493_6m','pro_4549_0m','post_4549_6m','pro_4651_0m','post_4651_6m',
	       'pro_778_0m','post_778_6m','post_778_12m','post_778_R', 'pro_4329_0m', 'post_4329_6m',
		'healthy_5319_0m', 'healthy_5320_0m', 'pro_4961_0m', 'post_4961_6m', 'pro_5070_0m', 'post_5070_6m',
		'pro_1093_0m','post_1093_6m', 'pro_3156_0m', 'post_3156_6m')

print(length(sample_order))



#sample_order=c('pro_778_0m','post_778_6m','post_778_12m','post_778_R')
		

seurat_qc<- function(sobject,name,folder,nFeature_RNA_min,nFeature_RNA_max,
                     nCount_RNA_min, nCount_RNA_max,percent_mito,percent_ribo){

#  sobject[["percent.mt"]] <- PercentageFeatureSet(sobject, pattern = "^MT-")
#  sobject[["percent.rb"]] <- PercentageFeatureSet(sobject, pattern = "^RP[SL]")

  counts_df<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='OG')
  colnames(counts_df)<- c('Read_Counts','Metric','Filter')
  genes_df<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='OG')
  colnames(genes_df)<- c('Gene_Counts','Metric','Filter')
  mito<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='OG')


  sobject <- subset(sobject, subset = nFeature_RNA>nFeature_RNA_min & nFeature_RNA<nFeature_RNA_max & percent.mt<percent_mito & percent.rb< percent_ribo &
                      nCount_RNA < nCount_RNA_max & nCount_RNA > nCount_RNA_min )

  counts_df_qc<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='QC')
  colnames(counts_df_qc)<- c('Read_Counts','Metric','Filter')
  genes_df_qc<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='QC')
  colnames(genes_df_qc)<- c('Gene_Counts','Metric','Filter')
  mito_qc<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='QC')

  counts_df<-rbind(counts_df,counts_df_qc)
  genes_df<-rbind(genes_df,genes_df_qc)
  mito<- rbind(mito,mito_qc)

  p1<-ggplot(counts_df, aes(x=Metric, y=Read_Counts)) + geom_violin(fill='firebrick4') +
    theme_bw() + facet_grid(cols  = vars(Metric), rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)

  p2<-ggplot(genes_df, aes(x=Metric, y=Gene_Counts)) + geom_violin(fill='turquoise4') +
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 )

  p3<-ggplot(mito, aes(x=Metric, y=percent.mt)) + geom_violin(fill='forestgreen') +
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)

#  setEPS()
  pdf(paste0(folder,name,".pdf"))
  plot(p1+p2+p3)
  dev.off()
  return(sobject)
}



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


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",'royalblue', 'grey','red','forestgreen','firebrick'
)




rename_list=list()

if(is.na(opt$matrix)){rename_list=readRDS(opt$rds)}

if(!is.na(opt$matrix)){
for(file in list.files(path)) {
        print(file)
        name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
        gex_list[[name]]<- Read10X_h5(paste0(path,'/',file))
                                }
        rename_list= lapply(X=gex_list, FUN=function(x){
        x=CreateSeuratObject(counts=x, assay="RNA", min.cells=10) })
                        }


infos=c()
cells=c()
samples=c()
for(n in names(rename_list)){
    if(n!="merged"){
    cell_number=nrow(rename_list[[n]]@meta.data)
    infos=c(infos,n)
    cells=c(cells,cell_number)
    samples=c(samples,strsplit(n, split="_")[[1]][2])    
    }}

df=data.frame(Samples=samples, Cell_Number=cells, Info=infos)
df$Info<- factor(df$Info, levels=sample_order)

p<-ggplot(df, aes(x=Info, y=Cell_Number, fill=Samples)) + 
geom_bar(stat = "identity") + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=13)) +
scale_fill_manual(values = c25)

infos=c()
cells=c()
samples=c()

for(n in names(rename_list)){
    if(n!="merged"){
    doublet_number=table(rename_list[[n]]@meta.data$scDblFinder.class)
        if(length(doublet_number)>1){
	    infos=c(infos,c(n,n))
	    samples=c( samples, c(strsplit(n, split="_")[[1]][2],strsplit(n, split="_")[[1]][2]) )
	    cells=c(cells,doublet_number) }
        else{
            infos=c(infos,n)
            samples=c(samples,strsplit(n, split="_")[[1]][2])
            cells=c(cells,doublet_number) }        
    }
}
print(cells)
df=data.frame(Samples=samples, Cell_Number=cells, Cell_Type=names(cells), Info=infos)
print(df)

df$Info<- factor(df$Info, levels=sample_order)


p1<-ggplot(df, aes(x=Info, y=Cell_Number, fill=Cell_Type)) + 
geom_bar(position="fill",stat = "identity") + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13), axis.text.y = element_text(size=13)) +
scale_fill_manual(values = c('steelblue', 'gold')) + theme(legend.text = element_text(size=14))
#p2<-ggplot(df, aes(x=Pools, y=Cell_Number, fill=Cell_Type)) + geom_bar(position="stack",stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
#scale_fill_manual(values = c('steelblue', 'gold')) + theme(legend.text = element_text(size=14))


metadata=rename_list$merged@meta.data
metadata$FullInfo<- factor(metadata$FullInfo, levels=sample_order)

p2<-metadata %>% 
  ggplot(aes(x=FullInfo, y=nFeature_RNA, fill=Sample)) + 
  geom_boxplot() + 
  theme_bw() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13), axis.text.y = element_text(size=14)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGenes") +
  scale_fill_manual(values = c25)


rename_list$merged[["percent.mt"]] <- PercentageFeatureSet(rename_list$merged, pattern = "^MT-")
metadata=rename_list$merged@meta.data
metadata$FullInfo<- factor(metadata$FullInfo,levels=sample_order) 

p3<-metadata %>%
  ggplot(aes(x=FullInfo, y=percent.mt, fill=Sample)) +
  geom_boxplot() +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13), axis.text.y = element_text(size=14)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("%MT") +
  scale_fill_manual(values = c25)

metadata<-  rename_list$merged@meta.data
metadata['Complexity']<- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)

p4<- metadata %>% 
  ggplot(aes(x=FullInfo, y=Complexity, fill=Sample)) + 
  geom_boxplot() + 
  theme_bw() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13), axis.text.y = element_text(size=14)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity") +
  scale_fill_manual(values = c25)


sobject <- subset(rename_list$merged, subset = nFeature_RNA>nFeature_RNA_min & nFeature_RNA<nFeature_RNA_max & percent.mt<percent_mito  &
                      nCount_RNA < nCount_RNA_max & nCount_RNA > nCount_RNA_min )

qc<-ifelse(rownames(metadata) %in% rownames(sobject@meta.data), "PASS", "FAIL")
metadata['QC']<- qc
doublets=rownames(metadata[metadata$scDblFinder.class=="doublet",])
sobject <- sobject[,!colnames(sobject) %in% doublets]
metadata <- metadata[!rownames(metadata) %in% doublets,]

df=data.frame()

for(n in unique(metadata$FullInfo)){
    
    tmp<-  as.data.frame(metadata %>% subset(FullInfo==n) %>% .$QC %>% table())
    tmp['FullInfo']<- n
    df =rbind(df,tmp)
}

colnames(df)<- c('QC', 'Cell_Number', 'FullInfo')
df$FullInfo<- factor(df$FullInfo, levels=sample_order)

p5<-ggplot(df, aes(x=FullInfo, y=Cell_Number, fill=QC)) + 
geom_bar(position="fill",stat = "identity") + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13), axis.text.y = element_text(size=14)) +
scale_fill_manual(values = c('steelblue', 'gold')) + theme(legend.text = element_text(size=14))


nam=c()
lg=c()
hg=c()
lr=c()
hr=c()
hm=c()
hrib=c()

for(i in names(rename_list)){
         current<-rename_list[[i]]
          current[["percent.mt"]] <- PercentageFeatureSet(current, pattern = "^MT-")
          current[["percent.rb"]] <- PercentageFeatureSet(current, pattern = "^RP[SL]")


        nam=c(nam,i)
        low_genes=length(which(current@meta.data$nFeature_RNA<QC[1]))/nrow(current@meta.data)
        lg=c(lg,low_genes)
        high_genes=length(which(current@meta.data$nFeature_RNA>QC[2]))/nrow(current@meta.data)
        hg=c(hg,high_genes)
        low_rna=length(which(current@meta.data$nCount_RNA<QC[3]))/nrow(current@meta.data)
        lr=c(lr,low_rna)
        high_rna=length(which(current@meta.data$nCount_RNA>QC[4]))/nrow(current@meta.data)
        hr=c(hr,high_rna)
        high_mito=length(which(current@meta.data$percent.mt>QC[5]))/nrow(current@meta.data)
        hm=c(hm,high_mito)
        high_ribo=length(which(current@meta.data$percent.rb>QC[6]))/nrow(current@meta.data)
        hrib=c(hrib,high_ribo)

}
qc_data=data.frame('Experiment'= nam, 'Low_Gene_Count'= lg, 'High_Gene_Count'=hg, 'Low_RNA_Count'=lr, 'High_RNA_Count'=hr, 'High_MT%'=hm, 'High_RB%'=hrib)
write.table(qc_data, paste0(output,'/QC_overview.tsv'), sep='\t', row.names = F, col.names = T)


p6<- VlnPlot(rename_list$merged,'nFeature_RNA', group.by='FullInfo', pt.size =0) + 
     VlnPlot(rename_list$merged,'nCount_RNA', group.by='FullInfo',pt.size =0) +
     VlnPlot(rename_list$merged,'percent.mt', group.by='FullInfo',pt.size =0)


if(opt$performQC=="TRUE"){

	ifelse(!dir.exists(paste0(output, '/perSample')), dir.create(paste0(output, '/perSample')), FALSE)
	rown=c()
	cbqc=c()
	caqc=c()
	qc_list=list()
	for(i in names(rename_list)){
	    if(i!="merged" ){
	        print(paste0("Sample to QC:", i))
	current_nd=rename_list[[i]]
	doublets=rownames(current_nd@meta.data[current_nd@meta.data$scDblFinder.class=="doublet",])
	current_nd <- current_nd[,!colnames(current_nd) %in% doublets]
	current_nd[["percent.mt"]] <- PercentageFeatureSet(current_nd, pattern = "^MT-")
	current_nd[["percent.rb"]] <- PercentageFeatureSet(current_nd, pattern = "^RP[SL]")
   


	if(opt$qlogical=="TRUE"){

		nFeature_RNA_min=quantile(current_nd@meta.data$nFeature_RNA, probs=1-as.numeric(opt$qthreshold))
		nFeature_RNA_max=quantile(current_nd@meta.data$nFeature_RNA, probs=as.numeric(opt$qthreshold))
		nCount_RNA_min=quantile(current_nd@meta.data$nCount_RNA, probs=1-as.numeric(opt$qthreshold))
		nCount_RNA_max=quantile(current_nd@meta.data$nCount_RNA, probs=as.numeric(opt$qthreshold))
		percent_mito=quantile(current_nd@meta.data$percent.mt, probs=as.numeric(opt$qthreshold))
		percent_rb=quantile(current_nd@meta.data$percent.rb, probs=as.numeric(opt$qthreshold))


#		nFeature_RNA_min=500
#		nFeature_RNA_max=4000
#		nCount_RNA_min=500
#		nCount_RNA_max=15000
#		percent_mito=15
#		percent_rb=25
		QC=c(nFeature_RNA_min,nFeature_RNA_max,nCount_RNA_min,nCount_RNA_max,percent_mito,percent_rb)
		print(QC)
				}

             current<-seurat_qc(current_nd,name=i,folder=paste0(output, '/perSample/'), 
				nFeature_RNA_min = QC[1], nFeature_RNA_max = QC[2], 
				nCount_RNA_min = QC[3], nCount_RNA_max = QC[4],
				percent_mito = QC[5],percent_ribo=QC[6])

          rown=c(rown,i)
          cbqc=c(cbqc,ncol(rename_list[[i]]))
          caqc=c(caqc,ncol(current))

           print(paste0(i,' Cells before QC: ', ncol(rename_list[[i]]) ))
           print(paste0(i,' Cells after QC: ', ncol(current) ))
                    
          qc_list[[i]]<- current
        }
 
summary=data.frame('dataset'=rown, 'Cells_Before_QC'=cbqc, 'Cells_after_QC'=caqc)
summary$pct=summary$Cells_after_QC / summary$Cells_Before_QC                        }
write.table(summary, paste0(output,'/QC_before_after.tsv'), sep='\t', row.names = F, col.names = T)

merged_seurat <- merge(x = qc_list[[1]],
                       y = qc_list[2:length(qc_list)],
                       merge.data = TRUE)

doublets=rownames(merged_seurat@meta.data[merged_seurat@meta.data$scDblFinder.class=="doublet",])
merged_seurat <- merged_seurat[,!colnames(merged_seurat) %in% doublets]

qc_list[['merged']]<- merged_seurat

p7<- VlnPlot(qc_list$merged,'nFeature_RNA', group.by='FullInfo', pt.size =0) +
     VlnPlot(qc_list$merged,'nCount_RNA', group.by='FullInfo',pt.size =0) +
     VlnPlot(qc_list$merged,'percent.mt', group.by='FullInfo',pt.size =0)



merged_seurat[["joined"]] <- JoinLayers(merged_seurat[["RNA"]])
merged_seurat <- NormalizeData(merged_seurat,assay = "joined")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_seurat.cc <- CellCycleScoring(merged_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "joined")


#merged_seurat <- NormalizeData(merged_seurat)
#merged_seurat <- CellCycleScoring(object = merged_seurat, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)

alld <- data.frame()
Idents(object = merged_seurat.cc)=merged_seurat.cc@meta.data$FullInfo
clusters=levels(Idents(object = merged_seurat.cc))
for(i in clusters){
    x=as.data.frame(merged_seurat.cc@meta.data %>% subset(FullInfo == i ) %>% .$Phase %>% table()/table(merged_seurat.cc$Phase))
    x$Freq <- x$Freq * 1/(sum(x$Freq))
    x$cluster<- i
    alld <- rbind(alld, x)
}

colnames(alld)<- c('Phase','Counts','SampleID')
alld$SampleID<- factor(alld$SampleID, levels=sample_order)

p8 <- ggplot(alld, aes(y=Counts, x=SampleID, fill=Phase)) +
  geom_bar(stat='identity')+scale_fill_manual(values = c("firebrick",'forestgreen','steelblue2')) +
  scale_y_continuous(expand = c(0,0)) +
  ylab('Freq') +theme(
    axis.text.x = element_text(angle=45, hjust=1, size = 13), axis.title=element_text(size=22), axis.text.y = element_text(angle=45, hjust=1, size = 15)) + 
    theme(legend.text = element_text(size=14), legend.title = element_text(size=15))
  

ggsave(filename=paste0(output, '/QC_post-filtering.pdf'), device='pdf', width=10, height=10, plot=p7)
ggsave(filename=paste0(output, '/CC_post-filtering.pdf'), device='pdf', width=10, height=10, plot=p8)

saveRDS(qc_list, paste0(output, '/objects/QC_list.rds.gz'), compress = "gzip")


}


ggsave(filename=paste0(output, '/cells_per_sample.pdf'), device='pdf', width=10, height=10, plot=p)
ggsave(filename=paste0(output, '/doublets_per_sample.pdf'), device='pdf', width=10, height=10, plot=p1)
ggsave(filename=paste0(output, '/genes_per_sample.pdf'), device='pdf', width=10, height=10, plot=p2)
ggsave(filename=paste0(output, '/MT_per_sample.pdf'), device='pdf', width=10, height=10, plot=p3)
ggsave(filename=paste0(output, '/Complexity_per_sample.pdf'), device='pdf', width=10, height=10, plot=p4)
ggsave(filename=paste0(output, '/QC_per_sample.pdf'), device='pdf', width=10, height=10, plot=p5)
ggsave(filename=paste0(output, '/QC_pre-filtering.pdf'), device='pdf', width=10, height=10, plot=p6)
