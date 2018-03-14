### Libraries
library(tidyverse)

### Import all the data ####
POP.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_POP.qassoc", na.strings = "NA", sep = "")

### Import all the data from temp_18CD ####
T18CD.1b_size.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_1b_size.qassoc", na.strings = "NA", sep = "")
T18CD.2b_size.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_2b_size.qassoc", na.strings = "NA", sep = "")
T18CD.age_maturity.df <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_age_maturity.qassoc", na.strings = "NA", sep = "")
T18CD.ctmax.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_ctmax.qassoc", na.strings = "NA", sep = "")
T18CD.day1b.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_day1b.qassoc", na.strings = "NA", sep = "")
T18CD.day2b.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_day2b.qassoc", na.strings = "NA", sep = "")
T18CD.mortality.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_mortality.assoc", na.strings = "NA", sep = "")
T18CD.mort_day.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_mort_day.qassoc", na.strings = "NA", sep = "")
T18CD.size_maturity.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_size_maturity.qassoc", na.strings = "NA", sep = "")
T18CD.totNeonates.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_18CD/temp_18CD_out.T18_totNeonates.qassoc", na.strings = "NA", sep = "")

#### POP ####

### Multiple test correction
# Add column corrected with Bonferroni method
POP.df$P_bonf <- p.adjust(POP.df$P, method = "bonferroni", n = length (POP.df$P))
# Add column corrected with FDR method
POP.df$P_fdr <- p.adjust(POP.df$P, method = "fdr", n = length (POP.df$P))

### Create a column with a unique SNP ID by merging the columns "CHR", "BP".
POP.df <- unite (POP.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
POP_man_all <- ggplot(data = POP.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("POP")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(POP_man_all)
ggsave("POP_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
POP_man_all <- ggplot(data = POP.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("POP FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(POP_man_all)
ggsave("POP_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
POP_top <- POP.df[order(POP.df$P),]
# Extract the list of unique contig
top_contig <- unique(POP_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
POP_top <- dplyr::filter(POP.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
POP_man_top <- ggplot(data = POP_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("POP - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(POP_man_top)
ggsave("POP_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
POP_man_top <- ggplot(data = POP_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("POP - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(POP_man_top)
ggsave("POP_man_top_fdr.png", device = NULL, width=10, height=7)



#### T18CD Size of the First Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.1b_size.df$P_bonf <- p.adjust(T18CD.1b_size.df$P, method = "bonferroni", n = length (T18CD.1b_size.df$P))
# Add column corrected with FDR method
T18CD.1b_size.df$P_fdr <- p.adjust(T18CD.1b_size.df$P, method = "fdr", n = length (T18CD.1b_size.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.1b_size.df <- unite (T18CD.1b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.1b_size_man_all <- ggplot(data = T18CD.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size of the First Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.1b_size_man_all)
ggsave("T18CD.1b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.1b_size_man_all <- ggplot(data = T18CD.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size of the First Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.1b_size_man_all)
ggsave("T18CD.1b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.1b_size_top <- T18CD.1b_size.df[order(T18CD.1b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.1b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.1b_size_top <- dplyr::filter(T18CD.1b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.1b_size_man_top <- ggplot(data = T18CD.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size of the First Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.1b_size_man_top)
ggsave("T18CD.1b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.1b_size_man_top <- ggplot(data = T18CD.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size of the First Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.1b_size_man_top)
ggsave("T18CD.1b_size_man_top_fdr.png", device = NULL, width=10, height=7)



#### T18CD Size of the Second Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.2b_size.df$P_bonf <- p.adjust(T18CD.2b_size.df$P, method = "bonferroni", n = length (T18CD.2b_size.df$P))
# Add column corrected with FDR method
T18CD.2b_size.df$P_fdr <- p.adjust(T18CD.2b_size.df$P, method = "fdr", n = length (T18CD.2b_size.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.2b_size.df <- unite (T18CD.2b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.2b_size_man_all <- ggplot(data = T18CD.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size of the Second Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.2b_size_man_all)
ggsave("T18CD.2b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.2b_size_man_all <- ggplot(data = T18CD.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size of the Second Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.2b_size_man_all)
ggsave("T18CD.2b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.2b_size_top <- T18CD.2b_size.df[order(T18CD.2b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.2b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.2b_size_top <- dplyr::filter(T18CD.2b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.2b_size_man_top <- ggplot(data = T18CD.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size of the Second Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.2b_size_man_top)
ggsave("T18CD.2b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.2b_size_man_top <- ggplot(data = T18CD.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size of the Second Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.2b_size_man_top)
ggsave("T18CD.2b_size_man_top_fdr.png", device = NULL, width=10, height=7)







#### T18CD Age at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.age_maturity.df$P_bonf <- p.adjust(T18CD.age_maturity.df$P, method = "bonferroni", n = length (T18CD.age_maturity.df$P))
# Add column corrected with FDR method
T18CD.age_maturity.df$P_fdr <- p.adjust(T18CD.age_maturity.df$P, method = "fdr", n = length (T18CD.age_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.age_maturity.df <- unite (T18CD.age_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.age_maturity_man_all <- ggplot(data = T18CD.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Age at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.age_maturity_man_all)
ggsave("T18CD.age_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.age_maturity_man_all <- ggplot(data = T18CD.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Age at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.age_maturity_man_all)
ggsave("T18CD.age_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.age_maturity_top <- T18CD.age_maturity.df[order(T18CD.age_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.age_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.age_maturity_top <- dplyr::filter(T18CD.age_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.age_maturity_man_top <- ggplot(data = T18CD.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Age at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.age_maturity_man_top)
ggsave("T18CD.age_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.age_maturity_man_top <- ggplot(data = T18CD.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Age at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.age_maturity_man_top)
ggsave("T18CD.age_maturity_man_top_fdr.png", device = NULL, width=10, height=7)




#### T18CD Temperature of Maximum Tolerance ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.ctmax.df$P_bonf <- p.adjust(T18CD.ctmax.df$P, method = "bonferroni", n = length (T18CD.ctmax.df$P))
# Add column corrected with FDR method
T18CD.ctmax.df$P_fdr <- p.adjust(T18CD.ctmax.df$P, method = "fdr", n = length (T18CD.ctmax.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.ctmax.df <- unite (T18CD.ctmax.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.ctmax_man_all <- ggplot(data = T18CD.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Temperature of Maximum Tolerance")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.ctmax_man_all)
ggsave("T18CD.ctmax_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.ctmax_man_all <- ggplot(data = T18CD.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Temperature of Maximum Tolerance FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.ctmax_man_all)
ggsave("T18CD.ctmax_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.ctmax_top <- T18CD.ctmax.df[order(T18CD.ctmax.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.ctmax_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.ctmax_top <- dplyr::filter(T18CD.ctmax.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.ctmax_man_top <- ggplot(data = T18CD.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Temperature of Maximum Tolerance - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.ctmax_man_top)
ggsave("T18CD.ctmax_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.ctmax_man_top <- ggplot(data = T18CD.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Temperature of Maximum Tolerance - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.ctmax_man_top)
ggsave("T18CD.ctmax_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD First Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.day1b.df$P_bonf <- p.adjust(T18CD.day1b.df$P, method = "bonferroni", n = length (T18CD.day1b.df$P))
# Add column corrected with FDR method
T18CD.day1b.df$P_fdr <- p.adjust(T18CD.day1b.df$P, method = "fdr", n = length (T18CD.day1b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.day1b.df <- unite (T18CD.day1b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.day1b_man_all <- ggplot(data = T18CD.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD First Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day1b_man_all)
ggsave("T18CD.day1b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.day1b_man_all <- ggplot(data = T18CD.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD First Brood Release Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day1b_man_all)
ggsave("T18CD.day1b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.day1b_top <- T18CD.day1b.df[order(T18CD.day1b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.day1b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.day1b_top <- dplyr::filter(T18CD.day1b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.day1b_man_top <- ggplot(data = T18CD.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD First Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day1b_man_top)
ggsave("T18CD.day1b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.day1b_man_top <- ggplot(data = T18CD.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD First Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day1b_man_top)
ggsave("T18CD.day1b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD Second Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.day2b.df$P_bonf <- p.adjust(T18CD.day2b.df$P, method = "bonferroni", n = length (T18CD.day2b.df$P))
# Add column corrected with FDR method
T18CD.day2b.df$P_fdr <- p.adjust(T18CD.day2b.df$P, method = "fdr", n = length (T18CD.day2b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.day2b.df <- unite (T18CD.day2b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.day2b_man_all <- ggplot(data = T18CD.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day2b_man_all)
ggsave("T18CD.day2b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.day2b_man_all <- ggplot(data = T18CD.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day2b_man_all)
ggsave("T18CD.day2b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.day2b_top <- T18CD.day2b.df[order(T18CD.day2b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.day2b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.day2b_top <- dplyr::filter(T18CD.day2b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.day2b_man_top <- ggplot(data = T18CD.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Second Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day2b_man_top)
ggsave("T18CD.day2b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.day2b_man_top <- ggplot(data = T18CD.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Second Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.day2b_man_top)
ggsave("T18CD.day2b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD Mortality ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.mortality.df$P_bonf <- p.adjust(T18CD.mortality.df$P, method = "bonferroni", n = length (T18CD.mortality.df$P))
# Add column corrected with FDR method
T18CD.mortality.df$P_fdr <- p.adjust(T18CD.mortality.df$P, method = "fdr", n = length (T18CD.mortality.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.mortality.df <- unite (T18CD.mortality.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.mortality_man_all <- ggplot(data = T18CD.mortality.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Mortality")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mortality_man_all)
ggsave("T18CD.mortality_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.mortality_man_all <- ggplot(data = T18CD.mortality.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Mortality FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mortality_man_all)
ggsave("T18CD.mortality_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.mortality_top <- T18CD.mortality.df[order(T18CD.mortality.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.mortality_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.mortality_top <- dplyr::filter(T18CD.mortality.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.mortality_man_top <- ggplot(data = T18CD.mortality_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Mortality - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mortality_man_top)
ggsave("T18CD.mortality_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.mortality_man_top <- ggplot(data = T18CD.mortality_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Mortality - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mortality_man_top)
ggsave("T18CD.mortality_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD Mortality Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.mort_day.df$P_bonf <- p.adjust(T18CD.mort_day.df$P, method = "bonferroni", n = length (T18CD.mort_day.df$P))
# Add column corrected with FDR method
T18CD.mort_day.df$P_fdr <- p.adjust(T18CD.mort_day.df$P, method = "fdr", n = length (T18CD.mort_day.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.mort_day.df <- unite (T18CD.mort_day.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.mort_day_man_all <- ggplot(data = T18CD.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Mortality Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mort_day_man_all)
ggsave("T18CD.mort_day_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.mort_day_man_all <- ggplot(data = T18CD.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Mortality Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mort_day_man_all)
ggsave("T18CD.mort_day_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.mort_day_top <- T18CD.mort_day.df[order(T18CD.mort_day.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.mort_day_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.mort_day_top <- dplyr::filter(T18CD.mort_day.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.mort_day_man_top <- ggplot(data = T18CD.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Mortality Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mort_day_man_top)
ggsave("T18CD.mort_day_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.mort_day_man_top <- ggplot(data = T18CD.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Mortality Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.mort_day_man_top)
ggsave("T18CD.mort_day_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD Size at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.size_maturity.df$P_bonf <- p.adjust(T18CD.size_maturity.df$P, method = "bonferroni", n = length (T18CD.size_maturity.df$P))
# Add column corrected with FDR method
T18CD.size_maturity.df$P_fdr <- p.adjust(T18CD.size_maturity.df$P, method = "fdr", n = length (T18CD.size_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.size_maturity.df <- unite (T18CD.size_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.size_maturity_man_all <- ggplot(data = T18CD.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.size_maturity_man_all)
ggsave("T18CD.size_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.size_maturity_man_all <- ggplot(data = T18CD.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.size_maturity_man_all)
ggsave("T18CD.size_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.size_maturity_top <- T18CD.size_maturity.df[order(T18CD.size_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.size_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.size_maturity_top <- dplyr::filter(T18CD.size_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.size_maturity_man_top <- ggplot(data = T18CD.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Size at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.size_maturity_man_top)
ggsave("T18CD.size_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.size_maturity_man_top <- ggplot(data = T18CD.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Size at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.size_maturity_man_top)
ggsave("T18CD.size_maturity_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18CD Total Number of Neonates ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18CD.totNeonates.df$P_bonf <- p.adjust(T18CD.totNeonates.df$P, method = "bonferroni", n = length (T18CD.totNeonates.df$P))
# Add column corrected with FDR method
T18CD.totNeonates.df$P_fdr <- p.adjust(T18CD.totNeonates.df$P, method = "fdr", n = length (T18CD.totNeonates.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18CD.totNeonates.df <- unite (T18CD.totNeonates.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18CD.totNeonates_man_all <- ggplot(data = T18CD.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Total Number of Neonates")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.totNeonates_man_all)
ggsave("T18CD.totNeonates_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18CD.totNeonates_man_all <- ggplot(data = T18CD.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Total Number of Neonates FDR") +
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.totNeonates_man_all)
ggsave("T18CD.totNeonates_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18CD.totNeonates_top <- T18CD.totNeonates.df[order(T18CD.totNeonates.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18CD.totNeonates_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18CD.totNeonates_top <- dplyr::filter(T18CD.totNeonates.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18CD.totNeonates_man_top <- ggplot(data = T18CD.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18CD Total Number of Neonates - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.totNeonates_man_top)
ggsave("T18CD.totNeonates_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18CD.totNeonates_man_top <- ggplot(data = T18CD.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18CD Total Number of Neonates - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18CD.totNeonates_man_top)
ggsave("T18CD.totNeonates_man_top_fdr.png", device = NULL, width=10, height=7)





