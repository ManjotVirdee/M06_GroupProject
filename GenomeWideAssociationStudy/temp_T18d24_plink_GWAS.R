# Author:
# Gian Marco Baranzoni

### Libraries
library(tidyverse)

### Import all the data ####

### Import all the data from T18d24 ####
T18d24.1b_size.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_1b_size.qassoc", na.strings = "NA", sep = "")
T18d24.2b_size.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_2b_size.qassoc", na.strings = "NA", sep = "")
T18d24.age_maturity.df <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_age_maturity.qassoc", na.strings = "NA", sep = "")
T18d24.ctmax.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_ctmax.qassoc", na.strings = "NA", sep = "")
T18d24.day1b.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_day1b.qassoc", na.strings = "NA", sep = "")
T18d24.day2b.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_day2b.qassoc", na.strings = "NA", sep = "")
T18d24.mort_day.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_mort_day.qassoc", na.strings = "NA", sep = "")
T18d24.size_maturity.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_size_maturity.qassoc", na.strings = "NA", sep = "")
T18d24.totNeonates.df   <- read.csv("/home/gmbaranzoni/Desktop/M06_group_project/Dapnia_project/Dmagna_GWAS/temp_T18d24/temp_T18d24_out.T18d24_totNeonates.qassoc", na.strings = "NA", sep = "")

#### T18d24 Size of the First Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.1b_size.df$P_bonf <- p.adjust(T18d24.1b_size.df$P, method = "bonferroni", n = length (T18d24.1b_size.df$P))
# Add column corrected with FDR method
T18d24.1b_size.df$P_fdr <- p.adjust(T18d24.1b_size.df$P, method = "fdr", n = length (T18d24.1b_size.df$P))

### Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.1b_size.df <- unite (T18d24.1b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.1b_size_man_all <- ggplot(data = T18d24.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size of the First Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.1b_size_man_all)
ggsave("T18d24.1b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.1b_size_man_all <- ggplot(data = T18d24.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size of the First Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.1b_size_man_all)
ggsave("T18d24.1b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.1b_size_top <- T18d24.1b_size.df[order(T18d24.1b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.1b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.1b_size_top <- dplyr::filter(T18d24.1b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.1b_size_man_top <- ggplot(data = T18d24.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size of the First Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.1b_size_man_top)
ggsave("T18d24.1b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.1b_size_man_top <- ggplot(data = T18d24.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size of the First Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.1b_size_man_top)
ggsave("T18d24.1b_size_man_top_fdr.png", device = NULL, width=10, height=7)


#### T18d24 Size of the Second Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.2b_size.df$P_bonf <- p.adjust(T18d24.2b_size.df$P, method = "bonferroni", n = length (T18d24.2b_size.df$P))
# Add column corrected with FDR method
T18d24.2b_size.df$P_fdr <- p.adjust(T18d24.2b_size.df$P, method = "fdr", n = length (T18d24.2b_size.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.2b_size.df <- unite (T18d24.2b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.2b_size_man_all <- ggplot(data = T18d24.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size of the Second Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.2b_size_man_all)
ggsave("T18d24.2b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.2b_size_man_all <- ggplot(data = T18d24.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size of the Second Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.2b_size_man_all)
ggsave("T18d24.2b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.2b_size_top <- T18d24.2b_size.df[order(T18d24.2b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.2b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.2b_size_top <- dplyr::filter(T18d24.2b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.2b_size_man_top <- ggplot(data = T18d24.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size of the Second Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.2b_size_man_top)
ggsave("T18d24.2b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.2b_size_man_top <- ggplot(data = T18d24.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size of the Second Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.2b_size_man_top)
ggsave("T18d24.2b_size_man_top_fdr.png", device = NULL, width=10, height=7)







#### T18d24 Age at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.age_maturity.df$P_bonf <- p.adjust(T18d24.age_maturity.df$P, method = "bonferroni", n = length (T18d24.age_maturity.df$P))
# Add column corrected with FDR method
T18d24.age_maturity.df$P_fdr <- p.adjust(T18d24.age_maturity.df$P, method = "fdr", n = length (T18d24.age_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.age_maturity.df <- unite (T18d24.age_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.age_maturity_man_all <- ggplot(data = T18d24.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Age at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.age_maturity_man_all)
ggsave("T18d24.age_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.age_maturity_man_all <- ggplot(data = T18d24.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Age at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.age_maturity_man_all)
ggsave("T18d24.age_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.age_maturity_top <- T18d24.age_maturity.df[order(T18d24.age_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.age_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.age_maturity_top <- dplyr::filter(T18d24.age_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.age_maturity_man_top <- ggplot(data = T18d24.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Age at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.age_maturity_man_top)
ggsave("T18d24.age_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.age_maturity_man_top <- ggplot(data = T18d24.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Age at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.age_maturity_man_top)
ggsave("T18d24.age_maturity_man_top_fdr.png", device = NULL, width=10, height=7)




#### T18d24 Temperature of Maximum Tolerance ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.ctmax.df$P_bonf <- p.adjust(T18d24.ctmax.df$P, method = "bonferroni", n = length (T18d24.ctmax.df$P))
# Add column corrected with FDR method
T18d24.ctmax.df$P_fdr <- p.adjust(T18d24.ctmax.df$P, method = "fdr", n = length (T18d24.ctmax.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.ctmax.df <- unite (T18d24.ctmax.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.ctmax_man_all <- ggplot(data = T18d24.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Temperature of Maximum Tolerance")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.ctmax_man_all)
ggsave("T18d24.ctmax_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.ctmax_man_all <- ggplot(data = T18d24.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Temperature of Maximum Tolerance FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.ctmax_man_all)
ggsave("T18d24.ctmax_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.ctmax_top <- T18d24.ctmax.df[order(T18d24.ctmax.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.ctmax_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.ctmax_top <- dplyr::filter(T18d24.ctmax.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.ctmax_man_top <- ggplot(data = T18d24.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Temperature of Maximum Tolerance - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.ctmax_man_top)
ggsave("T18d24.ctmax_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.ctmax_man_top <- ggplot(data = T18d24.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Temperature of Maximum Tolerance - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.ctmax_man_top)
ggsave("T18d24.ctmax_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24 First Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.day1b.df$P_bonf <- p.adjust(T18d24.day1b.df$P, method = "bonferroni", n = length (T18d24.day1b.df$P))
# Add column corrected with FDR method
T18d24.day1b.df$P_fdr <- p.adjust(T18d24.day1b.df$P, method = "fdr", n = length (T18d24.day1b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.day1b.df <- unite (T18d24.day1b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.day1b_man_all <- ggplot(data = T18d24.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 First Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day1b_man_all)
ggsave("T18d24.day1b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.day1b_man_all <- ggplot(data = T18d24.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 First Brood Release Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day1b_man_all)
ggsave("T18d24.day1b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.day1b_top <- T18d24.day1b.df[order(T18d24.day1b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.day1b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.day1b_top <- dplyr::filter(T18d24.day1b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.day1b_man_top <- ggplot(data = T18d24.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 First Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day1b_man_top)
ggsave("T18d24.day1b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.day1b_man_top <- ggplot(data = T18d24.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 First Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day1b_man_top)
ggsave("T18d24.day1b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24 Second Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.day2b.df$P_bonf <- p.adjust(T18d24.day2b.df$P, method = "bonferroni", n = length (T18d24.day2b.df$P))
# Add column corrected with FDR method
T18d24.day2b.df$P_fdr <- p.adjust(T18d24.day2b.df$P, method = "fdr", n = length (T18d24.day2b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.day2b.df <- unite (T18d24.day2b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.day2b_man_all <- ggplot(data = T18d24.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day2b_man_all)
ggsave("T18d24.day2b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.day2b_man_all <- ggplot(data = T18d24.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day2b_man_all)
ggsave("T18d24.day2b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.day2b_top <- T18d24.day2b.df[order(T18d24.day2b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.day2b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.day2b_top <- dplyr::filter(T18d24.day2b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.day2b_man_top <- ggplot(data = T18d24.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Second Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day2b_man_top)
ggsave("T18d24.day2b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.day2b_man_top <- ggplot(data = T18d24.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Second Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.day2b_man_top)
ggsave("T18d24.day2b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24 Mortality Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.mort_day.df$P_bonf <- p.adjust(T18d24.mort_day.df$P, method = "bonferroni", n = length (T18d24.mort_day.df$P))
# Add column corrected with FDR method
T18d24.mort_day.df$P_fdr <- p.adjust(T18d24.mort_day.df$P, method = "fdr", n = length (T18d24.mort_day.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.mort_day.df <- unite (T18d24.mort_day.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.mort_day_man_all <- ggplot(data = T18d24.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Mortality Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.mort_day_man_all)
ggsave("T18d24.mort_day_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.mort_day_man_all <- ggplot(data = T18d24.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Mortality Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.mort_day_man_all)
ggsave("T18d24.mort_day_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.mort_day_top <- T18d24.mort_day.df[order(T18d24.mort_day.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.mort_day_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.mort_day_top <- dplyr::filter(T18d24.mort_day.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.mort_day_man_top <- ggplot(data = T18d24.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Mortality Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.mort_day_man_top)
ggsave("T18d24.mort_day_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.mort_day_man_top <- ggplot(data = T18d24.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Mortality Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.mort_day_man_top)
ggsave("T18d24.mort_day_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24 Size at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.size_maturity.df$P_bonf <- p.adjust(T18d24.size_maturity.df$P, method = "bonferroni", n = length (T18d24.size_maturity.df$P))
# Add column corrected with FDR method
T18d24.size_maturity.df$P_fdr <- p.adjust(T18d24.size_maturity.df$P, method = "fdr", n = length (T18d24.size_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.size_maturity.df <- unite (T18d24.size_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.size_maturity_man_all <- ggplot(data = T18d24.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.size_maturity_man_all)
ggsave("T18d24.size_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.size_maturity_man_all <- ggplot(data = T18d24.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.size_maturity_man_all)
ggsave("T18d24.size_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.size_maturity_top <- T18d24.size_maturity.df[order(T18d24.size_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.size_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.size_maturity_top <- dplyr::filter(T18d24.size_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.size_maturity_man_top <- ggplot(data = T18d24.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Size at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.size_maturity_man_top)
ggsave("T18d24.size_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.size_maturity_man_top <- ggplot(data = T18d24.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Size at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.size_maturity_man_top)
ggsave("T18d24.size_maturity_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24 Total Number of Neonates ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24.totNeonates.df$P_bonf <- p.adjust(T18d24.totNeonates.df$P, method = "bonferroni", n = length (T18d24.totNeonates.df$P))
# Add column corrected with FDR method
T18d24.totNeonates.df$P_fdr <- p.adjust(T18d24.totNeonates.df$P, method = "fdr", n = length (T18d24.totNeonates.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24.totNeonates.df <- unite (T18d24.totNeonates.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24.totNeonates_man_all <- ggplot(data = T18d24.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Total Number of Neonates")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.totNeonates_man_all)
ggsave("T18d24.totNeonates_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24.totNeonates_man_all <- ggplot(data = T18d24.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Total Number of Neonates FDR") +
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.totNeonates_man_all)
ggsave("T18d24.totNeonates_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24.totNeonates_top <- T18d24.totNeonates.df[order(T18d24.totNeonates.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24.totNeonates_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24.totNeonates_top <- dplyr::filter(T18d24.totNeonates.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24.totNeonates_man_top <- ggplot(data = T18d24.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24 Total Number of Neonates - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.totNeonates_man_top)
ggsave("T18d24.totNeonates_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24.totNeonates_man_top <- ggplot(data = T18d24.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24 Total Number of Neonates - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24.totNeonates_man_top)
ggsave("T18d24.totNeonates_man_top_fdr.png", device = NULL, width=10, height=7)





