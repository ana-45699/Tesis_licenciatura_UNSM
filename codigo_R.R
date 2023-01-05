# Librerias

library("fastqcr")
library('ggplot2')
library(ggpubr)
library('tidyverse')
library(DT)
library(ggVennDiagram)
library("reticulate")
library("stringr")
library(scales)
library(RColorBrewer)

# Exploración de los archivos crudos (extensión bam)
read_bam <- read.delim("file_bam_reads_1.txt",sep = ' ',header = FALSE)
colnames(read_bam)<-c('Samples','Bam reads')
bam_bed<-merge(read_bam,min_bed,by='Samples',all=TRUE)
colnames(bam_bed)[3]<-'Bed'
datatable(bam_bed[order(bam_bed$`Bam reads`,decreasing = T),])
min_bed %>% ggplot(aes(x=Samples,y=`Number of reads`))+geom_bar(stat="identity",fill='#4682b4',color='black')+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle('Numero de reads por muestra')+theme(plot.title = element_text(hjust = 0.5))+ylab('Número de reads')+xlab('Muestras')
## Frecuencia de la cantidad de reads
min_bed <- as.data.frame(min_bed)
ggplot(min_bed,aes(y=`Number of reads`))+geom_boxplot()+ylab('Número de reads')+
  theme_classic() +  theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank())
quantile(min_bed$`Number of reads`)

# Tablas de starFusion
starFusion <- read.csv("all_samples.csv",sep = '\t')
temp_samples <- as.character(levels(starFusion$Sample)) 
temp_samples<-gsub('_1','',temp_samples)
starFusion$Sample <- gsub('_1','',starFusion$Sample)
starFusion$X.FusionName <- as.character(starFusion$X.FusionName)
filter_starFusion <- starFusion %>% filter(FFPM>=0.4045)
temp_samples_filter <- as.character(unique(filter_starFusion$Sample)) 
freq_genes_star <- starFusion  %>% select(X.FusionName,Sample) %>% unique() %>% select(X.FusionName) %>% table() %>% data.frame() %>% arrange(desc(Freq))
freq_genes_star_filter <- filter_starFusion %>% select(X.FusionName,Sample) %>% unique() %>% select(X.FusionName) %>% table() %>% data.frame() %>% arrange(desc(Freq))
starFusion_genes_raw <- read.delim('StarFusion_genes_0.txt')
starFusion_genes<- starFusion_genes_raw[,c(1,2,3)]
gene_A <-starFusion_genes[,c(1,2)]
colnames(gene_A)[2]<-'Gene'
gene_B <-starFusion_genes[,c(1,3)]
colnames(gene_B)[2]<-'Gene'
starFusion_genes_single <- rbind(gene_B,gene_A) 
single_A <- starFusion %>% select(LeftGene,Sample)
colnames(single_A)[1] <- c('Gene')
single_B <- starFusion %>% select(RightGene,Sample)
colnames(single_B)[1] <- c('Gene')
single <- rbind(single_A,single_B)
freq_single <- single %>% unique() %>% select(Gene) %>% table() %>% data.frame() %>% arrange(desc(Freq))
datatable(freq_genes_star,caption = 'STAR-Fusion Gene fusions of all samples ')
datatable(freq_single,caption = 'STAR-Fusion genes of all samples ')
# Tabla de frecuencia de genes compañeros
starFusion_unique_sample <- starFusion %>% distinct(X.FusionName,Sample,.keep_all = TRUE)
starFusion_single_5 <- starFusion_unique_sample[,c(7,9)]
colnames(starFusion_single_5)<-c('Gene',"Fusion_gene")
starFusion_single_3 <- starFusion_unique_sample[,c(9,1)]
colnames(starFusion_single_3)<-c('Gene',"Fusion_gene")
star_partner <- rbind(starFusion_single_5,starFusion_single_3)
star_partner <- as.data.frame(table(star_partner))%>%filter(Freq>0)
star_partner_1 <- as.data.frame(table(Gene=star_partner$Gene))

# Tablas de Fuseq

fuseq_raw <- read.delim('bed_results.csv',sep = '\t')
fuseq_raw$Sample<- gsub('/','',fuseq_raw$Sample)
fuseq_raw$Sample <- gsub('_1','',fuseq_raw$Sample)
fuseq_raw$Gene_Name <- paste(fuseq_raw$symbol5,fuseq_raw$symbol3)
## Pasar de formato fuseq a STARFusion
fuseq = fuseq_raw
fuseq_table = data.frame(FusionName = paste0(fuseq$symbol5,"--",fuseq$symbol3),
                         JunctionReadCount = fuseq$supportRead,
                         SpanningFragCount = fuseq$supportRead,
                         SpliceType = "ONLY_REF_SPLICE",
                         LeftGene = paste0(fuseq$symbol5,"^",fuseq$gene5),
                         LeftBreakpoint = paste(paste0("chr",fuseq$chrom5),fuseq$cds.brpos5.start,fuseq$strand5,sep=":"), 
                         RightGene = paste0(fuseq$symbol3,"^",fuseq$gene3),
                         RightBreakpoint = paste(paste0("chr",fuseq$chrom3),fuseq$cds.brpos3.start,fuseq$strand3,sep=":"),
                         LargeAnchorSupport = "YES_LDAS",
                         FFPM=0.1,
                         LeftBreakDinuc = "GT",
                         LeftBreakEntropy = 2.0,
                         RightBreakDinuc = "AG",
                         RightBreakEntropy = 2.0,
                         annots = "[\"NA\"]",
                         Sample = fuseq$Sample) 
muestra_conteo_bed <-fuseq_raw %>% dplyr::select(Sample) %>% table() %>% sort() %>%data.frame()
fusion_freq_fuseq <- data.frame(table(fuseq_raw$fusionName,fuseq_raw$Gene_Name))
fusion_freq_fuseq <- fusion_freq_fuseq %>% dplyr::filter(Freq>0)
colnames(fusion_freq_fuseq) <- c('Gene_ID','Gene_Name','Count')
fusion_freq_fuseq$Gene_Name <-str_replace(fusion_freq_fuseq$Gene_Name,' ','--')
genes_ID_fuseq <- c(as.character(fuseq_raw$gene5),as.character(fuseq_raw$gene3))
genes_name_fuseq <- c(as.character(fuseq_raw$symbol5),as.character(fuseq_raw$symbol3))
single_genes_fuseq <- data.frame('Gene'=genes_ID_fuseq,'Gene_Name'=genes_name_fuseq)
gene_freq_fuseq <- data.frame(read.delim('genes_count.txt'))
gene_freq_fuseq_table <- merge(single_genes_fuseq,gene_freq_fuseq,by='Gene')
gene_freq_fuseq_table <- gene_freq_fuseq_table %>% arrange(desc(Freq))%>% distinct(Gene,Freq,.keep_all = TRUE)
## Compañeros Fuseq
fuseq_partner <- as.data.frame(cbind(Gene=as.character(single_genes_fuseq$Gene),Names=paste0(fuseq$Gene_Name,"_",fuseq$fusionName)))
fuseq_partner <- as.data.frame(table(fuseq_partner))%>%filter(Freq>0)
fuseq_partner_1 <- as.data.frame(table(Gene=fuseq_partner$Gene))
genes_single_names <- unique(cbind(Gene_Name=as.character(single_genes_fuseq$Gene_Name),Gene=as.character(single_genes_fuseq$Gene)))
fuseq_partner_1 <- merge(fuseq_partner_1,genes_single_names,by="Gene")
## Gráfica de genes compañeros
ggdf = data.frame(fuseq_partner_1[,c(1,2)])
fuseq_graph<-ggplot(ggdf,aes(x=Freq))+
  theme_classic()+
  geom_histogram(aes(y=..density..),color="black",fill="#00008B",binwidth=1)+
  xlab("Número de genes compañeros de cada gen")+
  ylab("Frecuencia del número de compañeros")+
  scale_x_continuous(limits =c(0,10),n.breaks = 9)+
  #scale_y_continuous(limits=c(0,820),n.breaks = 10)+
  geom_density(alpha=.2, fill="#63B8FF")
ggdf = star_partner_1
star_graph<-ggplot(ggdf,aes(x=Freq))+
  theme_classic()+
  geom_histogram(aes(y=..density..),color="black",fill="#008B45",binwidth=1)+
  xlab("Número de genes compañeros de cada gen")+theme(axis.title.y = element_blank())+
  ylab("Frecuencia del número de compañeros")+
  scale_x_continuous(limits =c(0,10),n.breaks = 9)+
  #scale_y_continuous(limits=c(0,873),n.breaks = 10)
  geom_density(alpha=.2, fill="#63B8FF") 
ggarrange(fuseq_graph,star_graph,ncol = 2,nrow = 1)
# Frecuencia de genes que integran los genes fusion
gen_5_fuseq <-fuseq_raw[,c(1,12,19)]
colnames(gen_5_fuseq)<-c('Gen_ID','Gene_Name','Sample')
gen_3_fuseq <-fuseq_raw[,c(6,13,19)]
colnames(gen_3_fuseq)<-c('Gen_ID','Gene_Name','Sample')
single_gene_sample_fuseq <- rbind(gen_5_fuseq,gen_3_fuseq)
single_gene_sample_fuseq$Gene_Name <- as.character(single_gene_sample_fuseq$Gene_Name)
gene_sample_freq_fuseq <- as.data.frame(table(paste(single_gene_sample_fuseq$Gen_ID,single_gene_sample_fuseq$Gene_Name,single_gene_sample_fuseq$Sample,sep = "/")))%>% separate(Var1,c('Gen_ID','Gene_Name','Sample'),sep = '/')
freq_gene_one_fuseq <- gene_sample_freq_fuseq  %>% dplyr::select(Gen_ID, Gene_Name) %>% table() %>% data.frame() %>% dplyr::filter(Freq>0)
fuseq_raw$Gene_Name<- gsub(" ",'-',fuseq_raw$Gene_Name)

## Frecuencia de genes que integran las fusiones génicas Fuseq
ggdf = gene_freq_fuseq_table %>% dplyr::slice_head(n=15)
ggdf$Gene = factor(ggdf$Gene,levels=ggdf$Gene[order(ggdf$Freq)])
ggplot(ggdf,aes(x=Gene,y=Freq,fill=Gene))+geom_bar(stat="identity",color='black')+coord_flip()+
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.25),legend.position = "none")+
  ylab('Frecuencia')+
  xlab('Gene_ID')
datatable(fusion_freq_fuseq[order(fusion_freq_fuseq$Count,decreasing=T),],caption="Fusion genes from Fuseq")
## Frecuencia de fusiones génicas en Fuseq
ggdf = fusion_freq_fuseq %>% arrange(desc(Count))%>% slice_head(n=15)
ggdf$Gene_ID = factor(ggdf$Gene_ID,levels=ggdf$Gene_ID[order(ggdf$Count)])
ggdf$Gene_Name = factor(ggdf$Gene_Name,levels=ggdf$Gene_Name[order(ggdf$Count)]) 
colourCount = length(unique(ggdf$Gene_Name))
getPalette = colorRampPalette(brewer.pal(15, "GnBu"))
fuseq_graph<-ggplot(ggdf,aes(x=Gene_Name,y=Count,fill=Gene_Name))+
  geom_bar(stat="identity",color='black')+ ylim(c(0,90))+ ylab('Frecuencia')+coord_flip()+
  theme_classic() + scale_fill_manual(values = getPalette(colourCount))+
  theme(axis.title.y = element_blank(),plot.title = element_text(hjust = 0.25),legend.position = "none")
ggdf_1 = fuseq_raw  %>% filter(as.numeric(chrom3)==as.numeric(chrom5)) %>% select(chrom3) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
ggdf_2 = fuseq_raw  %>% filter(as.numeric(chrom3)!=as.numeric(chrom5))
ggdf_2_3 <- as.data.frame(ggdf_2[,7])
ggdf_2_5 <- as.data.frame(ggdf_2[,2])
colnames(ggdf_2_3)= colnames(ggdf_2_5) ="chrom"
ggdf_2_all= rbind(ggdf_2_3,ggdf_2_5) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
all_freq_chrom = c(ggdf_1$Freq + ggdf_2_all[-24,]$Freq,ggdf_2_all[24,]$Freq)
table_chrom <- data.frame(Chrom=as.character(ggdf_2_all$.),Freq=as.integer(all_freq_chrom))
table_chrom$Chrom<- factor(table_chrom$Chrom,levels=table_chrom$Chrom)
Fuseq_chrom<-table_chrom %>% ggplot(aes(x=Chrom,y=Freq))+ geom_bar(stat = "identity",color='#00008B',fill="#00008B",alpha=0.4)+theme_classic()+theme(axis.title.x = element_blank())+
  ylab("Frecuencia")

# Frecuencia de fusiones génicas en STARFusion
ggdf= freq_genes_star %>% arrange(desc(Freq))%>% slice_head(n=15)
colnames(ggdf)[1]='Gene_Fusion'
ggdf$Gene_Fusion = factor(ggdf$Gene_Fusion,levels=ggdf$Gene_Fusion[order(ggdf$Freq)])
colourCount = length(unique(ggdf$Gene_Fusion))
getPalette = colorRampPalette(brewer.pal(15, "GnBu"))
ggdf[is.na(ggdf)] = 0

star_graph<-ggplot(ggdf,aes(x=Gene_Fusion,y=Freq,fill=Gene_Fusion))+
  scale_fill_manual(values = getPalette(colourCount))+ ylim(c(0,90))+
  geom_bar(stat="identity",color='black')+ ylab('Frecuencia')+ coord_flip()+
  theme_classic() + 
  theme(axis.title.y = element_blank(),plot.title = element_text(hjust = 0.25),legend.position = "none")+
  #ggtitle('Top 15 Paired Fusion genes-STARFusion')
  theme(axis.title.y = element_blank(),plot.title = element_text(hjust = 0.25),legend.position = "none")
ggarrange(fuseq_graph,star_graph,ncol = 2,nrow=1)

## Gráfica de distribución de fusiones por cromosoma en Fuseq y STARFusion
# Fuseq
ggdf_1 = fuseq_raw  %>% filter(as.numeric(chrom3)==as.numeric(chrom5)) %>% select(chrom3) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
ggdf_2 = fuseq_raw  %>% filter(as.numeric(chrom3)!=as.numeric(chrom5))
ggdf_2_3 <- as.data.frame(ggdf_2[,7])
ggdf_2_5 <- as.data.frame(ggdf_2[,2])
colnames(ggdf_2_3)= colnames(ggdf_2_5) ="chrom"
ggdf_2_all= rbind(ggdf_2_3,ggdf_2_5) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
all_freq_chrom = c(ggdf_1$Freq + ggdf_2_all[-24,]$Freq,ggdf_2_all[24,]$Freq)
table_chrom <- data.frame(Chrom=as.character(ggdf_2_all$.),Freq=as.integer(all_freq_chrom))
table_chrom$Chrom<- factor(table_chrom$Chrom,levels=table_chrom$Chrom)
Fuseq_chrom<-table_chrom %>% ggplot(aes(x=Chrom,y=Freq))+ 
  geom_bar(stat = "identity",color='#00008B',fill="#00008B",alpha=0.4)+
  theme_classic()+theme(axis.title.x = element_blank())+
  ylab("Frecuencia")
#STAR
star_sep <- starFusion_unique_sample %>%
  separate(LeftBreakpoint,c('Chr5','Point_5','Strand_5'),sep = ':')%>%
  separate(RightBreakpoint,c('Chr3','Point_3','Strand_3'),sep = ':')
star_sep$Chr3 <- substr(star_sep$Chr3,4,6)
star_sep$Chr5 <- substr(star_sep$Chr5,4,6)
ggdf_1 = star_sep  %>% filter(Chr5==Chr3) %>% select(Chr3) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
ggdf_2 = star_sep  %>% filter(Chr3!=Chr5)
ggdf_2_3 <- as.data.frame(ggdf_2[,12])
ggdf_2_5 <- as.data.frame(ggdf_2[,8])
colnames(ggdf_2_3)= colnames(ggdf_2_5) ="chrom"
ggdf_2_all= rbind(ggdf_2_3,ggdf_2_5) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
all_freq_chrom = c(ggdf_1$Freq[1:23] + ggdf_2_all$Freq[1:23])
all_freq_chrom = c(all_freq_chrom,ggdf_2_all$Freq[24])
table_chrom <- data.frame(Chrom=as.character(ggdf_2_all$.),Freq=as.integer(all_freq_chrom))
table_chrom$Chrom<- factor(table_chrom$Chrom,levels=as.character(table_chrom$Chrom))
STAR_chrom<-table_chrom %>% ggplot(aes(x=Chrom,y=Freq))+ geom_bar(stat = "identity",color='#008B45',fill="#008B45",alpha=0.4)+theme_classic()+
  ylab("Frecuencia")+xlab("Cromosoma")

ggarrange(Fuseq_chrom,STAR_chrom,ncol = 1,nrow = 2)

# Gráfica de eventos intracromosomicos e inter cromosomicos

intra<-function(X){
  grepl("INTRACHROMOSOMAL",X,fixed = TRUE)
}
inter<-function(X){
  grepl("INTERCHROMOSOMAL",X,fixed = TRUE)
}
fuseq_intra <- fuseq %>% group_by(Sample) %>% summarise(event_per_sample=sum(as.numeric(chrom5)==as.numeric(chrom3),na.rm = TRUE))%>% mutate(type='Intracromosomal')
starFusion_unique_sample <- starFusion %>% distinct(X.FusionName,Sample,.keep_all = TRUE) 
star_intra <- starFusion_unique_sample %>% group_by(Sample) %>% summarise(event_per_sample= sum(sapply(annots,intra),na.rm = TRUE)) %>% mutate(type='Intracromosomal')
##
fuseq_inter <- fuseq %>% group_by(Sample) %>% summarise(event_per_sample=sum(as.numeric(chrom5)!=as.numeric(chrom3),na.rm = TRUE))%>% mutate(type='Intercromosomal')
star_inter <- starFusion_unique_sample %>% group_by(Sample) %>% summarise(event_per_sample= sum(sapply(annots,inter),na.rm = TRUE)) %>% mutate(type='Intercromosomal')
##
table_events_fuseq <- rbind(fuseq_inter,fuseq_intra)
table_events_star <- rbind(star_inter,star_intra)

fuseq_event <- ggplot(table_events_fuseq,aes(x=type,y=event_per_sample)) +
  geom_boxplot(color='#00008B') + theme_classic()+
  theme(axis.title.x = element_blank())+
  xlab("Tipos de rearreglos cromosomales") +
  stat_compare_means(method = "wilcox.test",label="p.format",label.y = 30)+
  ylab("Conteo")
star_event <- ggplot(table_events_star,aes(x=type,y=event_per_sample)) +
  geom_boxplot(color='#008B45') + theme_classic()+theme(axis.title.x = element_blank(),axis.title.y =element_blank())+
  stat_compare_means(method = "wilcox.test",label="p.format",label.y = 30)+
  ylab("Conteo")
ggarrange(fuseq_event,star_event,ncol = 2,nrow = 1)

# Resultados compartidos entre STARFusion y Fuseq

common_fusion_genes <-function(genes,Fuseq_genes,starFusion_genes){
  common_table <- 'NA'
  for (i in genes){
    fuseq_table_sub <- Fuseq_genes %>% dplyr::filter(Gene == i) %>% unique()
    starFusion_sub <- starFusion_genes %>% dplyr::filter(Gene == i) %>% unique()
    samples_vector<- intersect(fuseq_table_sub$Sample,starFusion_sub$Sample)
    if (length(samples_vector)>=1){
      common_table_sub <- data.frame(starFusion_sub[which(starFusion_sub$Sample %in% samples_vector),])
      common_table <- data.frame(rbind(common_table,common_table_sub))
    }}
  rownames(common_table) <- NULL
  compare_table<- data.frame(common_table[-1,])
  return(compare_table)
}
starFusion_genes_1 <- starFusion_genes_raw %>%
  unite(Gene,c('LeftGene','RightGene'))%>% dplyr::group_by(Sample,Gene) %>% dplyr::filter(FFPM==max(FFPM))
colnames(Fuseq_genes)[2]<-'Gene'
colnames(starFusion_genes)[2]<-'Gene'
total_genes_fusion <- unique(as.character(Fuseq_genes$Gene,starFusion_genes$Gene))
compare_fusion_genes <- common_fusion_genes(total_genes_fusion,Fuseq_genes,starFusion_genes_1)
raw_total_compare <- compare_fusion_genes %>% 
  separate(Gene,c('Pos_Gene_5','Pos_Gene_3'),sep = '_') 
datatable(compare_fusion_genes,caption = "Fusion Genes with share samples")

## Gráficos
ggdf= raw_total_compare %>% dplyr::select(Sample,LeftGene_Name,RigthGene_Name)%>%mutate(Fusion_Gene=paste0(LeftGene_Name,'--',RigthGene_Name))%>% dplyr::select(Sample,Fusion_Gene)%>%
  table() %>% as.data.frame() %>% dplyr::filter(Freq>0)
ggdf= as.data.frame(table(ggdf$Fusion_Gene))
colnames(ggdf)[1]<-"Fusion_Gene"
ggdf = ggdf %>% arrange(desc(Freq)) %>% dplyr::slice_head(n=15)
ggdf$Fusion_Gene = factor(ggdf$Fusion_Gene,levels=ggdf$Fusion_Gene[order(ggdf$Freq)])
svgFilename <- paste0('top_freq_common.svg')
svg(svgFilename,width = 7, height = 7)
colourCount = length(unique(ggdf$Fusion_Gene))
getPalette = colorRampPalette(brewer.pal(15, "GnBu"))
ggplot(ggdf,aes(x=Fusion_Gene,y=Freq,fill=Fusion_Gene))+geom_bar(stat="identity",color='black')+ylab('Frecuencia')+scale_fill_manual(values = getPalette(colourCount))+
  coord_flip()+
  theme_classic() + 
  theme(axis.title.y = element_blank(),plot.title = element_text(hjust = 0.25),legend.position = "none")
dev.off()
# Tabla completa de resultados compartidos
star_check<- read.csv("STAR_per_unique_sample.tsv",sep='\t')
star_check<-star_check  %>% mutate(search=paste0(LeftGene_ID,"-",RightGene_ID,"-",Sample))
fuseq_check<-fuseq_raw %>% mutate(search=paste0(gene5,"-",gene3,"-",Sample))
common_all_table <- merge(fuseq_check,star_check,by="search",all=FALSE)
# Chequear los FPM de los genes
FPM_check = common_all_table %>% select(X.FusionName,FFPM) %>% group_by(X.FusionName) %>% summarise(Media=mean(FFPM))

ggdf_1 = common_all_table %>% filter(as.numeric(chrom3)==as.numeric(chrom5)) %>% select(chrom3) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
ggdf_2 = common_all_table  %>% filter(as.numeric(chrom3)!=as.numeric(chrom5))
ggdf_2_3 <- as.data.frame(ggdf_2[,8])
ggdf_2_5 <- as.data.frame(ggdf_2[,3])
colnames(ggdf_2_3)= colnames(ggdf_2_5) ="chrom"
ggdf_2_all= rbind(ggdf_2_3,ggdf_2_5) %>% table() %>% as.data.frame() %>% mutate(order=as.numeric(as.character(.))) %>% arrange(order)
all_freq_chrom = c(ggdf_1$Freq + ggdf_2_all[-24,]$Freq,ggdf_2_all[24,]$Freq)
table_chrom <- data.frame(Chrom=as.character(ggdf_2_all$.),Freq=as.integer(all_freq_chrom))
table_chrom$Chrom<- factor(table_chrom$Chrom,levels=table_chrom$Chrom)
table_chrom %>% ggplot(aes(x=Chrom,y=Freq))+ geom_bar(stat = "identity",color='#3A5FCD',fill="#3A5FCD",alpha=0.4)+theme_classic()+
  ylab("Frecuencia")+xlab("Cromosoma")

common_intra <-common_all_table %>% group_by(Sample.x) %>% summarise(event_per_sample=sum(as.numeric(chrom5)==as.numeric(chrom3),na.rm = TRUE))%>% mutate(type='Intracromosomal')
##
common_inter <- common_all_table %>% group_by(Sample.x) %>% summarise(event_per_sample=sum(as.numeric(chrom5)!=as.numeric(chrom3),na.rm = TRUE))%>% mutate(type='Intercromosomal')
table_events_common <- rbind(common_inter,common_intra)
ggplot(table_events_common,aes(x=type,y=event_per_sample)) +
  geom_boxplot(color='#3A5FCD') + theme_classic()+
  theme(axis.title.x = element_blank())+
  xlab("Tipos de rearreglos cromosomales") +stat_compare_means(method = "t.test",label="p.format",label.y = 15)+
  ylab("Conteo")

## Chequeo en bases de datos de fusiones génicas
fusion_GDB <- read.csv("combinedFGDB2genes_genes_ID_new.txt",
                       sep='\t',header = FALSE)
fusion_cosmo <- read.table("cosmic_data_base.txt",header = TRUE)
fusion_cosmo <- fusion_cosmo %>% separate(Genes,c("Gen_5","Gen_3"),sep="-") %>%
  separate(Gen_5,c("Gen_5","Gen_5_ID"),sep="_") %>% 
  separate(Gen_3,c("Gen_3","Gen_3_ID"),sep="_")
check_GDB <- intersect(unique(c(as.character(fusion_GDB$V3),as.character(fusion_GDB$V5))),
                       unique(c(as.character(common_all_table$symbol3),
                                as.character(common_all_table$Symbol5))))
check_cosmo <- intersect(unique(c(as.character(fusion_cosmo$Gen_5),
                                  as.character(fusion_cosmo$Gen_3))),
                         unique(c(as.character(common_all_table$symbol3),
                                  as.character(common_all_table$Symbol5))))

# Filtro de kinasas de los resultados compartidos

kinase_tab <- read.delim("kinases.txt")
kinase_fuseq <- unique(merge(single_genes_fuseq,kinase_tab,by='Gene_Name'))
single_star_gene_Name <-data.frame(Gene=c(as.character(starFusion_genes_raw$LeftGene),
                                          as.character(starFusion_genes_raw$RightGene)) ,
                                   Gene_Name=c(as.character(starFusion_genes_raw$LeftGene_Name),
                                               as.character(starFusion_genes_raw$RigthGene_Name)))
kinase_star <- unique(merge(single_star_gene_Name,kinase_tab,by='Gene_Name'))

kinase_fuseq$Gene <- as.character(kinase_fuseq$Gene)
compare_genes_table$Gene<- as.character(compare_genes_table$Gene)

kinases_share_5 <- merge(kinase_fuseq,raw_total_compare,by.x="Gene_Name",by.y="LeftGene_Name")
kinases_share_3 <- merge(kinase_fuseq,raw_total_compare,by.x="Gene_Name",by.y="RigthGene_Name")

kinases_share_total <- c(paste0(kinases_share_5$Gene_Name,'--',kinases_share_5$RigthGene_Name),paste0(kinases_share_3$LeftGene_Name,'--',kinases_share_3$Gene_Name))
k_fusion_freq <- as.data.frame(table(kinases_share_total))

ggdf = k_fusion_freq
ggdf = as.data.frame(table(ggdf$Freq)) %>%
  mutate(pcnt = (Freq / sum(Freq))*100) %>% 
  arrange(pcnt)%>% mutate(etiquetas = paste0(as.character(round(pcnt,0)),'%'))
ggdf %>% 
  ggplot(aes(x= "",y=pcnt,fill=Var1))+ labs(fill = paste0("Número de",
                                                          '\n','Compañeros'))+scale_fill_manual(values=c("#00C5CD","#00868B"))+ geom_col(color='black',alpha=0.7) +
  geom_label(aes(label = paste0('Porcentaje: ',
                                etiquetas,'\n','Conteo: ',Freq),alpha=0.3,size=7),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y")+ theme_void()


## Figura de las kinasas con sus respectivos dominios
fusion_esqueme <- function(data){
  # Tamaño del marco
  begin=end=NULL
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0.40,0.65)
  p <- p + ggplot2::xlim(min(data$Start, na.rm=TRUE)-1,max(data$End, na.rm=TRUE)+1)
  p <-  p + ggplot2::labs(x = "Número de aminoácidos")
  p <- p + ggplot2::labs(y = "")
  # Longitud de la cadena
  p<- p +ggplot2::ggtitle(data$Gen_Fusion_label[1])+ 
    ggplot2::geom_rect(data = data[data$Type == "GENE",],
                       mapping=ggplot2::aes(xmin = min(data$Start, na.rm=TRUE),
                                            xmax = max(data$End,
                                                       na.rm=TRUE),ymin= 0.5,ymax=0.55),
                       colour = "black",fill = '#7A8B8B',size = 0.5,alpha =1.0)
  # Gráfico de dominios no kinasa
  domains_1=data[data$Type == "DOMAIN",]
  breakpoint = data[data$Type == "BREAK",]
  domains_np = domains_1[-(grep('KINASE_DOM',domains_1$Domain_label)),]
  colourCount = length(unique(domains_np$Domain_label))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))
  if (dim(domains_np)[1]>0){
    for (i in seq(length(domains_np$End))){
      if (domains_np$End[i]>breakpoint$Start && domains_np$Start[i]<breakpoint$Start){
        domains_np$End[i]=breakpoint$Start
      }
    }
    p <-p + ggplot2::geom_rect(data=domains_np ,
                               mapping=ggplot2::aes(xmin=Start,
                                                    xmax=End,
                                                    ymin=0.47,
                                                    ymax=0.58,fill=Domain_label))+
      scale_fill_manual(values = getPalette(colourCount)) 
  }
  # Gráfico de dominios kinasa
  domain_p <-domains_1[(grep('KINASE_DOM',domains_1$Domain_label)),]
  if (domain_p$End>breakpoint$Start && domain_p$Start< breakpoint$Start){
    domain_p$End=breakpoint$Start
  }
  p <- p + ggplot2::geom_rect(domain_p,
                              mapping=ggplot2::aes(xmin=Start,
                                                   xmax=End,
                                                   ymin=0.45,
                                                   ymax=0.60),fill="#FFEC8B",colour = "black")+geom_text(
                                                     data=domain_p,mapping=aes(x=Start+((End-Start)/2), y=0.52,
                                                                               label="Kinase"),size=4)
  # Breakpoint
  p <- p + geom_vline(data = breakpoint, mapping=aes(xintercept=Start),
                      color="#CD0000",linetype = "longdash")
  # Detalles de forma
  p <- p + theme(
    axis.ticks = element_blank(),panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()) +
    theme(axis.text.y = element_blank()) + theme(panel.border = element_blank(),
                                                 panel.background = element_rect(fill ="white"),legend.position = "bottom")
  p <- p + theme(plot.title.position = "panel",
                 plot.title = element_text(color = "black",face="bold"),
                 axis.title.x = element_text(size = rel(0.7)),
                 legend.text = element_text(size = rel(0.6)),
                 legend.title = element_blank(),
                 legend.key.size = unit(0.5,"line")
  )
  return(p)
}

table_graph <- read.delim('kinase_table_graphics_2.txt',sep = '\t')
fusion_genes <- table_graph %>% dplyr::filter(Type=='GENE')
A <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='BRIX1--MAP2K5',])

B <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='EIF2B3--MAST2',])

C <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='GOLGA4--BRAF',])
D <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='ROCK2--CACNA2D3',])

E <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='NLK--ACACA',])
F <- fusion_esqueme(table_graph[table_graph$Gen_Fusion_label=='NLK--HLF',])

svgFilename <- paste0('kinasas_final_prosite.svg')
svg(svgFilename,width = 10, height = 5)
ggarrange(A,B,C,D,E,F,ncol = 2,nrow = 3)
dev.off()