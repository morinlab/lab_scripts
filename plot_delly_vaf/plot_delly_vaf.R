
# Import libraries --------------------------------------------------------

library("tidyr")
library("magrittr")
library("dplyr")
library("ggplot2")

# Parse command-line arguments --------------------------------------------

options(echo=TRUE)
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]

# Create data frame -------------------------------------------------------

df <- read.table(input_file, sep="\t", header=TRUE)
tumour_names <- colnames(df)[8:9]
cols <- colnames(df)
cols[8:9] <- c("tumour_1", "tumour_2")
colnames(df) <- cols
head(df)


# Split adta frame into three based on zero values ------------------------

df_t1_zero <- filter(df, tumour_1 == 0)
df_t1_zero <- mutate(df_t1_zero, tumour_1 = -0.05)
df_t2_zero <- filter(df, tumour_2 == 0)
df_t2_zero <- mutate(df_t2_zero, tumour_2 = -0.05)
df_no_zero <- filter(df, tumour_1 != 0, tumour_2 != 0)

# Plot VAFs ---------------------------------------------------------------

plot <- ggplot() + 
    geom_point(data = df_no_zero, aes(tumour_1, tumour_2, colour = factor(type)), alpha = 0.5, size = 1) +
    geom_point(data = df_t1_zero, aes(tumour_1, tumour_2, colour = factor(type)), alpha = 0.5, size = 1,
               position = position_jitter(w = 0.02, h = 0)) +
    geom_point(data = df_t2_zero, aes(tumour_1, tumour_2, colour = factor(type)), alpha = 0.5, size = 1,
               position = position_jitter(w = 0, h = 0.02)) +
    coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
ggsave(file="delly_vaf_t1_vs_t2.pdf", width=6.5, height=5)
