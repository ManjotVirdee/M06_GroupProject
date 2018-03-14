### Libraries
library(tidyverse)

### Import all the data ####

### Import all the data from T18d24CD ####

T18d24CD.1b_size.df       <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_1b_size.qassoc", sep="")
T18d24CD.2b_size.df       <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_2b_size.qassoc", sep="")
T18d24CD.age_maturity.df  <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_age_maturity.qassoc", sep="")
T18d24CD.ctmax.df         <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_ctmax.qassoc", sep="")
T18d24CD.day1b.df         <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_day1b.qassoc", sep="")
T18d24CD.day2b.df         <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_day2b.qassoc", sep="")
T18d24CD.mort_day.df      <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_mort_day.qassoc", sep="")
T18d24CD.size_maturity.df <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_size_maturity.qassoc", sep="")
T18d24CD.totNeonates.df   <- read.csv("~/Desktop/Dmagna_GWAS/temp_T18d24CD/temp_T18d24CD_out.T18d24_totNeonates.qassoc", sep="")

#### T18d24CD Size of the First Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.1b_size.df$P_bonf <- p.adjust(T18d24CD.1b_size.df$P, method = "bonferroni", n = length (T18d24CD.1b_size.df$P))
# Add column corrected with FDR method
T18d24CD.1b_size.df$P_fdr <- p.adjust(T18d24CD.1b_size.df$P, method = "fdr", n = length (T18d24CD.1b_size.df$P))

### Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.1b_size.df <- unite (T18d24CD.1b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.1b_size_man_all <- ggplot(data = T18d24CD.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size of the First Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.1b_size_man_all)
ggsave("T18d24CD.1b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.1b_size_man_all <- ggplot(data = T18d24CD.1b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size of the First Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.1b_size_man_all)
ggsave("T18d24CD.1b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.1b_size_top <- T18d24CD.1b_size.df[order(T18d24CD.1b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.1b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.1b_size_top <- dplyr::filter(T18d24CD.1b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.1b_size_man_top <- ggplot(data = T18d24CD.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size of the First Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.1b_size_man_top)
ggsave("T18d24CD.1b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.1b_size_man_top <- ggplot(data = T18d24CD.1b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size of the First Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.1b_size_man_top)
ggsave("T18d24CD.1b_size_man_top_fdr.png", device = NULL, width=10, height=7)


#### T18d24CD Size of the Second Brood ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.2b_size.df$P_bonf <- p.adjust(T18d24CD.2b_size.df$P, method = "bonferroni", n = length (T18d24CD.2b_size.df$P))
# Add column corrected with FDR method
T18d24CD.2b_size.df$P_fdr <- p.adjust(T18d24CD.2b_size.df$P, method = "fdr", n = length (T18d24CD.2b_size.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.2b_size.df <- unite (T18d24CD.2b_size.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.2b_size_man_all <- ggplot(data = T18d24CD.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size of the Second Brood")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.2b_size_man_all)
ggsave("T18d24CD.2b_size_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.2b_size_man_all <- ggplot(data = T18d24CD.2b_size.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size of the Second Brood FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.2b_size_man_all)
ggsave("T18d24CD.2b_size_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top 50 contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.2b_size_top <- T18d24CD.2b_size.df[order(T18d24CD.2b_size.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.2b_size_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.2b_size_top <- dplyr::filter(T18d24CD.2b_size.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.2b_size_man_top <- ggplot(data = T18d24CD.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size of the Second Brood - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.2b_size_man_top)
ggsave("T18d24CD.2b_size_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.2b_size_man_top <- ggplot(data = T18d24CD.2b_size_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size of the Second Brood - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.2b_size_man_top)
ggsave("T18d24CD.2b_size_man_top_fdr.png", device = NULL, width=10, height=7)







#### T18d24CD Age at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.age_maturity.df$P_bonf <- p.adjust(T18d24CD.age_maturity.df$P, method = "bonferroni", n = length (T18d24CD.age_maturity.df$P))
# Add column corrected with FDR method
T18d24CD.age_maturity.df$P_fdr <- p.adjust(T18d24CD.age_maturity.df$P, method = "fdr", n = length (T18d24CD.age_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.age_maturity.df <- unite (T18d24CD.age_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.age_maturity_man_all <- ggplot(data = T18d24CD.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Age at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.age_maturity_man_all)
ggsave("T18d24CD.age_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.age_maturity_man_all <- ggplot(data = T18d24CD.age_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Age at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.age_maturity_man_all)
ggsave("T18d24CD.age_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.age_maturity_top <- T18d24CD.age_maturity.df[order(T18d24CD.age_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.age_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.age_maturity_top <- dplyr::filter(T18d24CD.age_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.age_maturity_man_top <- ggplot(data = T18d24CD.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Age at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.age_maturity_man_top)
ggsave("T18d24CD.age_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.age_maturity_man_top <- ggplot(data = T18d24CD.age_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Age at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.age_maturity_man_top)
ggsave("T18d24CD.age_maturity_man_top_fdr.png", device = NULL, width=10, height=7)




#### T18d24CD Temperature of Maximum Tolerance ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.ctmax.df$P_bonf <- p.adjust(T18d24CD.ctmax.df$P, method = "bonferroni", n = length (T18d24CD.ctmax.df$P))
# Add column corrected with FDR method
T18d24CD.ctmax.df$P_fdr <- p.adjust(T18d24CD.ctmax.df$P, method = "fdr", n = length (T18d24CD.ctmax.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.ctmax.df <- unite (T18d24CD.ctmax.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.ctmax_man_all <- ggplot(data = T18d24CD.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Temperature of Maximum Tolerance")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.ctmax_man_all)
ggsave("T18d24CD.ctmax_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.ctmax_man_all <- ggplot(data = T18d24CD.ctmax.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Temperature of Maximum Tolerance FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.ctmax_man_all)
ggsave("T18d24CD.ctmax_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.ctmax_top <- T18d24CD.ctmax.df[order(T18d24CD.ctmax.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.ctmax_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.ctmax_top <- dplyr::filter(T18d24CD.ctmax.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.ctmax_man_top <- ggplot(data = T18d24CD.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Temperature of Maximum Tolerance - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.ctmax_man_top)
ggsave("T18d24CD.ctmax_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.ctmax_man_top <- ggplot(data = T18d24CD.ctmax_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Temperature of Maximum Tolerance - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.ctmax_man_top)
ggsave("T18d24CD.ctmax_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24CD First Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.day1b.df$P_bonf <- p.adjust(T18d24CD.day1b.df$P, method = "bonferroni", n = length (T18d24CD.day1b.df$P))
# Add column corrected with FDR method
T18d24CD.day1b.df$P_fdr <- p.adjust(T18d24CD.day1b.df$P, method = "fdr", n = length (T18d24CD.day1b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.day1b.df <- unite (T18d24CD.day1b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.day1b_man_all <- ggplot(data = T18d24CD.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD First Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day1b_man_all)
ggsave("T18d24CD.day1b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.day1b_man_all <- ggplot(data = T18d24CD.day1b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD First Brood Release Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day1b_man_all)
ggsave("T18d24CD.day1b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.day1b_top <- T18d24CD.day1b.df[order(T18d24CD.day1b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.day1b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.day1b_top <- dplyr::filter(T18d24CD.day1b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.day1b_man_top <- ggplot(data = T18d24CD.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD First Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day1b_man_top)
ggsave("T18d24CD.day1b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.day1b_man_top <- ggplot(data = T18d24CD.day1b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD First Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day1b_man_top)
ggsave("T18d24CD.day1b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24CD Second Brood Release Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.day2b.df$P_bonf <- p.adjust(T18d24CD.day2b.df$P, method = "bonferroni", n = length (T18d24CD.day2b.df$P))
# Add column corrected with FDR method
T18d24CD.day2b.df$P_fdr <- p.adjust(T18d24CD.day2b.df$P, method = "fdr", n = length (T18d24CD.day2b.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.day2b.df <- unite (T18d24CD.day2b.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.day2b_man_all <- ggplot(data = T18d24CD.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day2b_man_all)
ggsave("T18d24CD.day2b_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.day2b_man_all <- ggplot(data = T18d24CD.day2b.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Second Brood Release Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day2b_man_all)
ggsave("T18d24CD.day2b_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.day2b_top <- T18d24CD.day2b.df[order(T18d24CD.day2b.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.day2b_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.day2b_top <- dplyr::filter(T18d24CD.day2b.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.day2b_man_top <- ggplot(data = T18d24CD.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Second Brood Release Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day2b_man_top)
ggsave("T18d24CD.day2b_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.day2b_man_top <- ggplot(data = T18d24CD.day2b_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Second Brood Release Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.day2b_man_top)
ggsave("T18d24CD.day2b_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24CD Mortality Day ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.mort_day.df$P_bonf <- p.adjust(T18d24CD.mort_day.df$P, method = "bonferroni", n = length (T18d24CD.mort_day.df$P))
# Add column corrected with FDR method
T18d24CD.mort_day.df$P_fdr <- p.adjust(T18d24CD.mort_day.df$P, method = "fdr", n = length (T18d24CD.mort_day.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.mort_day.df <- unite (T18d24CD.mort_day.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.mort_day_man_all <- ggplot(data = T18d24CD.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Mortality Day")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.mort_day_man_all)
ggsave("T18d24CD.mort_day_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.mort_day_man_all <- ggplot(data = T18d24CD.mort_day.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Mortality Day FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.mort_day_man_all)
ggsave("T18d24CD.mort_day_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.mort_day_top <- T18d24CD.mort_day.df[order(T18d24CD.mort_day.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.mort_day_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.mort_day_top <- dplyr::filter(T18d24CD.mort_day.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.mort_day_man_top <- ggplot(data = T18d24CD.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Mortality Day - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.mort_day_man_top)
ggsave("T18d24CD.mort_day_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.mort_day_man_top <- ggplot(data = T18d24CD.mort_day_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Mortality Day - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.mort_day_man_top)
ggsave("T18d24CD.mort_day_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24CD Size at Maturity ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.size_maturity.df$P_bonf <- p.adjust(T18d24CD.size_maturity.df$P, method = "bonferroni", n = length (T18d24CD.size_maturity.df$P))
# Add column corrected with FDR method
T18d24CD.size_maturity.df$P_fdr <- p.adjust(T18d24CD.size_maturity.df$P, method = "fdr", n = length (T18d24CD.size_maturity.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.size_maturity.df <- unite (T18d24CD.size_maturity.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.size_maturity_man_all <- ggplot(data = T18d24CD.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size at Maturity")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.size_maturity_man_all)
ggsave("T18d24CD.size_maturity_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.size_maturity_man_all <- ggplot(data = T18d24CD.size_maturity.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size at Maturity FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.size_maturity_man_all)
ggsave("T18d24CD.size_maturity_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.size_maturity_top <- T18d24CD.size_maturity.df[order(T18d24CD.size_maturity.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.size_maturity_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.size_maturity_top <- dplyr::filter(T18d24CD.size_maturity.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.size_maturity_man_top <- ggplot(data = T18d24CD.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Size at Maturity - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.size_maturity_man_top)
ggsave("T18d24CD.size_maturity_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.size_maturity_man_top <- ggplot(data = T18d24CD.size_maturity_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Size at Maturity - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.size_maturity_man_top)
ggsave("T18d24CD.size_maturity_man_top_fdr.png", device = NULL, width=10, height=7)





#### T18d24CD Total Number of Neonates ####

### Multiple test correction
# Add column corrected with Bonferroni method
T18d24CD.totNeonates.df$P_bonf <- p.adjust(T18d24CD.totNeonates.df$P, method = "bonferroni", n = length (T18d24CD.totNeonates.df$P))
# Add column corrected with FDR method
T18d24CD.totNeonates.df$P_fdr <- p.adjust(T18d24CD.totNeonates.df$P, method = "fdr", n = length (T18d24CD.totNeonates.df$P))
# Create a column with a unique SNP ID by merging the columns "CHR", "BP".
T18d24CD.totNeonates.df <- unite (T18d24CD.totNeonates.df, "CHR_BP", c(CHR, BP), sep = "_", remove = FALSE)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal P values
T18d24CD.totNeonates_man_all <- ggplot(data = T18d24CD.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Total Number of Neonates")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.totNeonates_man_all)
ggsave("T18d24CD.totNeonates_man_all.png", device = NULL, width=10, height=7)

### Plots all contigs with all values
# Manhattan plot, all SNPs, normal with FDR correction (0.05 significance level)
T18d24CD.totNeonates_man_all <- ggplot(data = T18d24CD.totNeonates.df) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Total Number of Neonates FDR") +
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.totNeonates_man_all)
ggsave("T18d24CD.totNeonates_man_all_fdr.png", device = NULL, width=10, height=7)



### Selection of top contigs with less than 
# Order the dataframe by lowest P value
T18d24CD.totNeonates_top <- T18d24CD.totNeonates.df[order(T18d24CD.totNeonates.df$P),]
# Extract the list of unique contig
top_contig <- unique(T18d24CD.totNeonates_top$CHR)
# Create a new dataframe keeping only the top contigs (in the previous order)
T18d24CD.totNeonates_top <- dplyr::filter(T18d24CD.totNeonates.df, CHR %in% top_contig[1:50])

### Plots top contigs
# Manhattan plot, top contigs, normal P values
T18d24CD.totNeonates_man_top <- ggplot(data = T18d24CD.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   ggtitle("T18d24CD Total Number of Neonates - Top contigs")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.totNeonates_man_top)
ggsave("T18d24CD.totNeonates_man_top.png", device = NULL, width=10, height=7)

### Plots top contigs
# Manhattan plot, top contigs, FDR
T18d24CD.totNeonates_man_top <- ggplot(data = T18d24CD.totNeonates_top) + theme_bw() +
   geom_point(mapping = aes(x = CHR_BP,
                            y = -log10(P_fdr),
                            fill = as.factor(CHR)),
              shape = 21,
              size = 2) + 
   geom_hline(yintercept = -log10(0.05), colour = "red") +
   ggtitle("T18d24CD Total Number of Neonates - Top contigs FDR")+
   xlab("SNPs") +
   ylab("-log10 (p-value)") +
   theme (axis.text.x = element_blank(), legend.position = "none")

print(T18d24CD.totNeonates_man_top)
ggsave("T18d24CD.totNeonates_man_top_fdr.png", device = NULL, width=10, height=7)





