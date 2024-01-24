####################################################################
## Análisis de los resultados de HISAT2 y STRINGTIE               ##
## usando los paquetes de R ballgown y limma.                     ##  
##                                                                ##
## Autor: Francisco J. Romero-Campero fran@us.es                  ##
####################################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para 
## realizar un análisis de expresión génica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciación realizados con 
## hisat2 y stringtie. 

## Para ejecutar con éxito este script es necesario descargar la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el 
## Working Directory To Source File Location. 

## Instalación y carga de los paquetes necesarios. Sólo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastará cargar los paquetes simplemente. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ballgown")

library(ballgown)


## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra típicamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y demás caracteriśticas de cada muestra. 
pheno.data <- read.csv("01-documentation/pheno_data.csv")
pheno.data

## La función ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentra las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio. 
bg.data <- ballgown(dataDir = "03-data/samples/", samplePattern = "sample", pData=pheno.data)
bg.data
sampleNames(bg.data)

## La función gexpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)
gene.names <- rownames(gene.expression)

## Nombramos las columnas con los nombres de nuestras muestras. 
colnames(gene.expression) <- c("h_alta1","h_alta2","h_alta3","h_alta4","h_baja1","h_baja2","h_baja3","h_baja4")

## Por motivos técnicos sumamos 1 a todos los niveles de expresión. 
## El problema viene provocado por x < 1 --> log2(x) < 0
gene.expression.1 <- gene.expression + 1

## Guardamos los datos de expresión génica sin procesar
write.table(x = gene.expression.1,file = "04-results/gene_expression_unprocessed.tsv",
            quote = F,row.names = F,
            sep = "\t")

## Representación de la distribución global de la expresión génica
color=c(rep("cyan4",4), rep("grey",4))
boxplot(gene.expression, outline=F,col=color,ylab="Gene Expression (FPKM)",
        cex.lab=1.5)

boxplot(log2(gene.expression.1), outline=F,col=color,
        ylab="log2(Gene Expression)",
        cex.lab=1.5)

## En este punto ballgown no realiza ninguna normalización de los datos más alla
## del cálculo de los niveles de expresión por FPKM.
## Utilizamos el paquete de R NormalyzerDE para esta tarea. Para ello es necesario generar
## un fichero con un formato específico.

library(NormalyzerDE)

design <- data.frame(sample=colnames(gene.expression),
                     group=c(rep("h_alta",4),rep("h_baja",4)))

write.table(x = design,file = "01-documentation/normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "AT_Humedad_norm",designPath = "01-documentation/normalyzer_design.tsv",
           dataPath = "04-results/gene_expression_unprocessed.tsv",outputDir = "04-results/"[])
# Tras ver los resultados, parece ser que CycLoess es el mejor método de normalización pues consigue una baja variabilidad intragrupo, un alto coeficiente de correlación intragrupo, agrupa bien las muestras y, además, no elimina diferencias significativas entre condiciones.

normalized.gene.expression <- read.table(file="04-results/AT_Humedad_norm/Quantile-normalized.txt", header=T)
head(normalized.gene.expression)
rownames(normalized.gene.expression) <- gene.names
 
## Representación de la distribución global de la expresión génica tras 
## normalización
boxplot(normalized.gene.expression, outline=F,col=color,
        ylab="log2(FPKM + 1)",cex.lab=1.5)

## Previsualizamos la similitud entre las réplicas


scatter_replicates <- function(obj,comparation_matrix, color=c("cyan4","grey")) {
  par(mfrow=c(ncol(comparation_matrix),c(nrow(comparation_matrix)/ncol(comparation_matrix))))
  for (i in 1:nrow(comparation_matrix)) {
    a=comparation_matrix[i,1]
    b=comparation_matrix[i,2]
    index_color=1
    if (i>c(nrow(comparation_matrix)/ncol(comparation_matrix))) {
      index_color=2
    }
    plot(x = obj[,a],
         y = obj[,b],
         pch=19,col=color[index_color],xlab=colnames(obj)[a],ylab=colnames(obj)[b],cex=0.5)
    text(x=3,y=10,
         labels = paste(c(
           "R2 = ",
           round(100*cor(obj[,a],
                         obj[,b]),
                 digits = 2),
           "%"), collapse=""))
  }
  par(mfrow=c(1,1))
}

comparaciones = matrix(c(1,2,1,3,1,4,2,3,2,4,3,4,5,6,5,7,5,8,6,7,6,8,7,8), nrow = 12, ncol = 2, byrow = T)
scatter_replicates(normalized.gene.expression,comparaciones)


## Realizamos un análisis de componentes principales y clustering 
## jerárquico para continuar con la exploración de los datos.
library(FactoMineR)
library(factoextra)
library(cluster)

## Los dos paquetes anteriores son geniales para análisis estadísticos en
## datos de secuencaición de nueva generación. Tienen una documentación muy 
## completa: http://www.sthda.com/english/articles/tag/factominer/
## Documentación genérica: http://www.sthda.com/english/

## Por ejemplo para PCA y clustering jerárquico:
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/65-pca-in-r-using-factominer-quick-scripts-and-videos/
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/74-hcpc-using-factominer-video/

pca.gene.expression <- data.frame(colnames(normalized.gene.expression),
                                  t(normalized.gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"
head(pca.gene.expression)

## PCA
res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_screeplot(res.pca)
fviz_pca_ind(res.pca, col.ind = c(rep("Humedad Alta",4),rep("Humedad Baja",4)), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

## Cluster jerárquico
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 2)   

fviz_dend(res.hcpc,k=2,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
)

res.agnes <- cluster::agnes(pca.gene.expression, method = "ward")
fviz_dend(res.agnes,k=2,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
)


## Cluster partición
fviz_nbclust(x=pca.gene.expression[,-1], method = "gap_stat",FUNcluster = cluster::pam, k=5)
fviz_nbclust(x=pca.gene.expression[,-1], method = "wss",FUNcluster = cluster::pam, k=5)
fviz_nbclust(x=pca.gene.expression[,-1], method = "silhouette",FUNcluster = cluster::pam, k=5)
res.pam <- cluster::pam(pca.gene.expression, k=2)


## Calculamos la matrix de expresión media. 
h_alta <- apply(normalized.gene.expression[1:4], MARGIN = 1,FUN = mean)
h_baja <- apply(normalized.gene.expression[5:8], MARGIN = 1,FUN = mean)

mean.expression <- matrix(c(h_alta,h_baja),ncol=2)
colnames(mean.expression) <- c("h_alta","h_baja")
rownames(mean.expression) <- rownames(normalized.gene.expression)
head(mean.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
plot(h_baja,h_alta,pch=19,cex=0.7,xlab="Humedad Baja",
     ylab=substitute(italic("Humedad Alta")),
     cex.lab=1.25,
     col="grey")

##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

## Especificamos el diseño experimental

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,1,2,2,2,2)))
colnames(experimental.design) <- c("h_baja","h_alta")

##A continuación, ajustamos la estimación de los niveles de expresión de cada
##gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
##fundamentalmente se calcula la media de las réplicas en cada condición.

# Cambio el orden de las columnas para poner la humedad baja como referencia. 
normalized.gene.expression <- normalized.gene.expression[c(5:8,1:4)]
linear.fit <- lmFit(normalized.gene.expression, experimental.design)

##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:

contrast.matrix <- makeContrasts(h_alta-h_baja,levels=c("h_baja","h_alta"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(normalized.gene.expression)

control.h_alta <- topTable(contrast.results, number=7507,coef=1,sort.by="logFC")
head(control.h_alta)

log.fold.change <- control.h_alta$logFC
q.value <- control.h_alta$adj.P.Val
genes.ids <- rownames(control.h_alta)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes <- genes.ids[log.fold.change > 2 & q.value < 0.05]
repressed.genes <- genes.ids[log.fold.change < - 2 & q.value < 0.05]

length(activated.genes)
length(repressed.genes)

## Volcano plot
log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-12,12),ylim = c(0,4), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5)

points(x = log.fold.change[activated.genes],
       y = log.q.val[activated.genes],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes],
       y = log.q.val[repressed.genes],col="blue",cex=0.8,pch=19)


## Enriquecimiento funcional. 
library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(activated.atha.enrich.go,showCategory = 20, font.size=8)
dotplot(activated.atha.enrich.go,showCategory = 20, font.size=8)
emapplot(pairwise_termsim(activated.atha.enrich.go),showCategory = 20, repel=T, cex_label_category=0.6)
cnetplot(activated.atha.enrich.go,showCategory = 20,repel=T, cex_label_category=0.6,cex_label_gene=0.4, max.overlaps=5)


repressed.atha.enrich.go <- enrichGO(gene          = repressed.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")
barplot(activated.atha.enrich.go,showCategory = 20, font.size=8)
dotplot(activated.atha.enrich.go,showCategory = 20, font.size=8, group=T)
emapplot(pairwise_termsim(activated.atha.enrich.go),showCategory = 20, repel=T, cex_label_category=0.6)
cnetplot(activated.atha.enrich.go,showCategory = 20,repel=T, cex_label_category=0.6,cex_label_gene=0.4, max.overlaps=5)


activated.atha.enrich.kegg <- enrichKEGG(gene  = activated.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)
df.activated.atha.enrich.kegg <- as.data.frame(activated.atha.enrich.kegg)
head(df.activated.atha.enrich.kegg)


repressed.atha.enrich.kegg <- enrichKEGG(gene  = repressed.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)

df.repressed.atha.enrich.kegg <- as.data.frame(repressed.atha.enrich.kegg)
head(df.repressed.atha.enrich.kegg)

## Podemos visualizar las rutas KEGG usando el paquete pathveiw
library(pathview)
rutas=c("ath01100","ath01110","ath00030", "ath00190","ath00510", "ath00513", "ath00905","ath04075", "ath04626", "ath00010", "ath00513", "ath00514")

for (i in rutas) {
  pathview(gene.data = sort(log.fold.change,decreasing = TRUE),
           pathway.id = i,
           species = "ath",
           limit = list(gene=max(abs(log.fold.change)), cpd=1),gene.idtype = "TAIR",kegg.dir="04-results/pathview/")
}


coolmap(normalized.gene.expression)

## Código para desarrollar una función gráfico de barras

gene <- activated.genes[1]

expression.matrix <- 2^normalized.gene.expression - 1
cond.names <- c("Humedad Baja","Humedad Alta")

gene.expression.barplot <- function(gene, expression.matrix,cond.names)
{
  expr.1 <- unlist(c(expression.matrix[gene, 1:4]))
  expr.2 <- unlist(c(expression.matrix[gene, 5:8]))
  
  mean.1 <- mean(expr.1)
  mean.2 <- mean(expr.2)
  
  sd.1 <- sd(expr.1)
  sd.2 <- sd(expr.2)
  
  means <- c(mean.1, mean.2)
  sds <- c(sd.1, sd.2)
  
  arrow.top <- means + sds
  arrow.bottom <- means - sds
  
  
  xpos <- barplot(means,ylim=c(0,1.5*max(arrow.top)),col=rainbow(2),
                  main=gene,names.arg = cond.names,
                  ylab="FPKM")
  arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.05)
}

## Análisis de expresión génica diferencial con DESeq2
library(DESeq2)

pheno.data

gene.count.matrix <- read.table(file = "gene_count_matrix.csv",header = T,sep = ",")
head(gene.count.matrix)

sapply(X = strsplit(x = gene.count.matrix$gene_id,split = "\\|"),FUN = function(x){return(x[1])})

gene.ids <- sapply(X = strsplit(x = gene.count.matrix$gene_id,split = "\\|"),FUN = function(x){return(x[1])})

gene.count.matrix <- gene.count.matrix[,-1]
rownames(gene.count.matrix) <- gene.ids
head(gene.count.matrix)

dds <- DESeqDataSetFromMatrix(countData=gene.count.matrix, colData=pheno.data, design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "col0")

dds <- DESeq(dds) 
res <- results(dds) 
res

log.fold.change <- res$log2FoldChange
q.value <- res$padj
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.deseq2 <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]

repressed.genes.deseq2 <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]

length(activated.genes.deseq2)
length(repressed.genes.deseq2)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(activated.atha.enrich.go,showCategory = 20)
dotplot(activated.atha.enrich.go,showCategory = 20)
emapplot(activated.atha.enrich.go,showCategory = 20)
cnetplot(activated.atha.enrich.go,showCategory = 20)


