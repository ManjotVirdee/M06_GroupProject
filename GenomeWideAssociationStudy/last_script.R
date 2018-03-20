###### GENE LIST MANIPULATIONS
### Libraries to load
library(tidyverse)
library(VennDiagram)

### Import manjot files ####
T18_top20MV <- read.delim("~/Desktop/Results/T18_top20_Genes_per_association_test.txt")
T18d24_top20MV <- read.delim("~/Desktop/Results/T18d24_top20_Genes_per_association_test.txt")

### Annotation ####
# Import the annotation
gff_scaff <- read.table("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/ref/reference_genome/Dmagna_scaffID.gff", sep="\t", quote="")

### extract only mRNA and separate the attributes in column 9
#gff_scaff2 <-  dplyr::filter(gff_scaff, V3 == "mRNA") %>% tidyr::separate(V9, c("ID", "Parent", "Note", "gbkey", "product"), sep = ";")
gff_scaff2 <-  dplyr::filter(gff_scaff, V3 == "mRNA") %>% tidyr::separate(V9, c("V9", "product"), sep = "product=")
gff_scaff$V1 <- as.factor(gff_scaff$V1)

# Intitialise column for two products and two IDs
T18_top20MV$product <- rep(NA, length(T18_top20MV$CHR))
T18d24_top20MV$product <- rep(NA, length(T18d24_top20MV$CHR))



for (i in 1:length(T18_top20MV$CHR)) {
   
   # Annotate T18
   chr <- T18_top20MV[i, "CHR"]  
   gff_temp <- dplyr::filter(gff_scaff2, V1 == as.character(chr))
   gff_temp <- dplyr::filter(gff_temp, V4 <= T18_top20MV[i, "BP"])
   gff_temp <- dplyr::filter(gff_temp, V5 >= T18_top20MV[i, "BP"])
   
   T18_top20MV[i, "product"] <-  gff_temp[1, "product"]

   # Annotate T18
   chr <- T18d24_top20MV[i, "CHR"]  
   gff_temp <- dplyr::filter(gff_scaff2, V1 == as.character(chr))
   gff_temp <- dplyr::filter(gff_temp, V4 <= T18d24_top20MV[i, "BP"])
   gff_temp <- dplyr::filter(gff_temp, V5 >= T18d24_top20MV[i, "BP"])
   
   T18d24_top20MV[i, "product"] <-  gff_temp[1, "product"]   
   
   
   print (length(T18_top20MV$CHR) - i)
}

write.table (T18_top20MV, "~/Desktop/Results/T18_top20_GeneProduct.txt", sep="\t", row.names=FALSE)
write.table (T18d24_top20MV, "~/Desktop/Results/T18d24_top20_GeneProduct.txt", sep="\t", row.names=FALSE)

### Histograms ####
### Plot histogram gene name T18
T18_hist.df <- T18_top20MV %>% dplyr::select(gene_name, label) %>% unique
T18_hist.df <- plyr::count(as.character(T18_hist.df$gene_name))
plyr::count(T18_hist.df$freq)

T18_hist <- ggplot(data = T18_hist.df, aes(x = x, y = freq)) + theme_bw() +
   geom_bar(stat="identity") + 
   xlab("genes") +
   ylab("n. of life history trait") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18_hist)
ggsave("~/Desktop/Results/T18_hist.png", device = NULL, width=10, height=7)



### Plot histogram gene name T18d24
T18d24_hist.df <- T18d24_top20MV %>% dplyr::select(gene_name, label) %>% unique
T18d24_hist.df <- plyr::count(as.character(T18d24_hist.df$gene_name))
plyr::count(T18d24_hist.df$freq)

T18d24_hist <- ggplot(data = T18d24_hist.df, aes(x = x, y = freq)) + theme_bw() +
   geom_bar(stat="identity") + 
   xlab("genes") +
   ylab("n. of life history trait") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24_hist)
ggsave("~/Desktop/Results/T18d24_hist.png", device = NULL, width=10, height=7)






### write list of genes that affect multiple phenotipes ####
which_genes <- T18_hist.df %>% dplyr::filter(freq > 1) %>% dplyr::select(x)
T18_multiPhen <- T18_top20MV %>% dplyr::filter(gene_name %in% which_genes$x) %>% dplyr::select(gene_name, label) %>% unique %>% dplyr::arrange(gene_name)
write.table (T18_multiPhen, "~/Desktop/Results/T18_multiPhen.txt", sep="\t", row.names=FALSE)

which_genes <- T18d24_hist.df %>% dplyr::filter(freq > 1) %>% dplyr::select(x)
T18d24_multiPhen <- T18d24_top20MV %>% dplyr::filter(gene_name %in% which_genes$x) %>% dplyr::select(gene_name, label) %>% unique %>% dplyr::arrange(gene_name)
write.table (T18d24_multiPhen, "~/Desktop/Results/T18_multiPhen.txt", sep="\t", row.names=FALSE)


### Venn Diagram 
T18_genes_num <- length(T18_hist.df$x)
T18d24_genes_num <- length(T18d24_hist.df$x)

overlap_gene_list <- dplyr::intersect(T18_hist.df$x, T18d24_hist.df$x)
overlap_gene_num <- length(overlap_gene_list)

T18_genes_list <- dplyr::setdiff(T18_hist.df$x, T18d24_hist.df$x)
T18d24_genes_list <- dplyr::setdiff(T18d24_hist.df$x, T18_hist.df$x)


grid.newpage()
Venn <- draw.pairwise.venn(T18_genes_num, T18d24_genes_num, overlap_gene_num, category = c("Evolution", "Plasticity"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))



# Create table for gene products
T18_genes_tab <- T18_top20MV %>% dplyr::filter(gene_name %in% T18_genes_list) %>% dplyr::select(CHR, gene_name, product) %>% unique
T18_genes_tab$label <- rep("evolution", length(T18_genes_tab$gene_name))

T18d24_genes_tab <- T18d24_top20MV %>% dplyr::filter(gene_name %in% T18d24_genes_list) %>% dplyr::select(CHR, gene_name, product) %>% unique
T18d24_genes_tab$label <- rep("plasticity", length(T18d24_genes_tab$gene_name))

common_genes_tab <- T18d24_top20MV %>% dplyr::filter(gene_name %in% overlap_gene_list) %>% dplyr::select(CHR, gene_name, product) %>% unique
common_genes_tab$label <- rep("both", length(common_genes_tab$gene_name))

gene_products.df <- dplyr::bind_rows(T18_genes_tab, T18d24_genes_tab, common_genes_tab, .id = NULL)

length(T18_genes_tab$gene_name)
length(T18d24_genes_tab$gene_name)
length(common_genes_tab$gene_name)
length(gene_products.df$gene_name)



write.table (gene_products.df, "~/Desktop/Results/gene_products.txt", sep="\t", row.names=FALSE)


scaffold.df <-  plyr::count(as.character(gene_products.df$CHR))
scaffold.df <- scaffold.df %>% dplyr::arrange(freq)
plyr::count(scaffold.df$freq)
View(scaffold.df)

hist (scaffold_hist$freq)

scaffold_hist <- ggplot(data = scaffold.df, aes(freq)) + theme_bw() +
   geom_histogram (breaks=seq(0.5, 9.5, by=0.5)) + 
   coord_cartesian(ylim = c(0, 20)) +
   xlab("n. of genes") +
   ylab("n. of scaffolds/contigs") +
   scale_x_continuous(breaks=seq(0,9,1)) +
   theme (axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          legend.position = "none")

print(scaffold_hist)
ggsave("~/Desktop/Results/scaffold_hist.png", device = NULL, width=10, height=7)


