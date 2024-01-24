## Plot linkage disequilibrium and distance
# 01/24/2024

library(tidyverse)

# Objective: Plot linkage disequilibrium and distance between markers.

# Read in data from TASSEL (calculated full matrix, heterozygous set to missing)
tassel_ld <- read.delim("Data/snps/prelim_merged_5899_LD_full.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
str(tassel_ld)

# Remove NAs for r^2 and distance (ie. not on same chromosome)
t_ld_chr <- tassel_ld %>%
  filter(Dist_bp != "N/A" & R.2 != "NaN") %>%
  mutate(Dist_bp = as.numeric(Dist_bp))

# Test plot
#ld_plot <- ggplot(data = t_ld_chr, mapping = aes(x = Dist_bp, y = R.2)) + geom_point(alpha = 0.05)
#pdf("Figures/ld_dist_5899_chrom_only.pdf", height = 11, width = 8.5)
#ld_plot
#dev.off()

# So many points and so slow!! Commenting out, don't run again.

# Try removing extra information to make plotting faster, and average within windows.
t_ld_chr <- t_ld_chr %>%
  select(Dist_bp, R.2)

# Make bins for intervals of 100 Kb
t_ld_chr$distc <- cut(t_ld_chr$Dist_bp, breaks = seq(from = min(t_ld_chr$Dist_bp)-1, to = max(t_ld_chr$Dist_bp)+1, by = 100000))

# Calculate average r2, count number of observations in bins, extract start position
av_dist_r2 <- t_ld_chr %>%
  filter(!is.na(distc)) %>%
  group_by(distc) %>%
  summarise(mean = mean(R.2), n = n()) %>%
  mutate(start = as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")))

# Plot
av_ld_plot <- ggplot(data = av_dist_r2, mapping = aes(x = start, y = mean, colour = n)) + 
  geom_point() + 
  geom_smooth(colour = "#1F497D") + 
  labs(x = "Distance (Mb)", y = expression("Average Linkage Disequilibrium "~(r^{2})), colour = paste0("Number of", "\n", "Pairwise Marker", "\n", "Combinations")) +
  scale_colour_gradient(low = "#4BACC6", high = "#1F497D") +
  scale_x_continuous(breaks = c(0,1*10^8, 2*10^8, 3*10^8, 4*10^8, 5*10^8), labels = c("0", "100", "200", "300", "400", "500")) + 
  theme_light()
av_ld_plot

ggsave("Figures/avg_ld_dist_5899_per100Kb.svg", height = 11, width = 8.5)
ggsave("Figures/avg_ld_dist_5899_per100Kb.pdf", height = 11, width = 8.5)

# Since plot only shows LD within chromosomes, it would be good to report the average LD across chromosomes (if there is high linkage disequilibrium across chromosomes, family structure is likely strong).

# Filter out NAs, calculate mean r2 for pairwise comparisons of chromosomes
compare_across_chr <- tassel_ld %>%
  filter(R.2 != "NaN") %>%
  group_by(Locus1, Locus2) %>%
  summarise(mean = mean(R.2))

across_chr_plot <- ggplot(data = compare_across_chr, mapping = aes(x = Locus2, y = reorder(Locus1, desc(Locus1)), fill = mean)) + 
  geom_tile(colour = "black") + 
  geom_text(aes(label = paste0(sprintf("%.5f", mean))), colour = "white") +
  coord_fixed() + 
  labs(x = "Chromsome", y = "Chromosome", fill = expression("Average LD"~(r^{2}))) +
  scale_fill_gradient(low = "#4BACC6", high = "#1F497D") + 
  theme(panel.background = element_blank(), axis.ticks = element_blank())
across_chr_plot

ggsave("Figures/avg_ld_across_chrom_5899.svg", height = 8.5, width = 11)
ggsave("Figures/avg_ld_across_chrom_5899.pdf", height = 8.5, width = 11)