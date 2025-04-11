## code modified from: https://github.com/vincenthahaut/FLASH-Seq/blob/main/R/timing_graph/timing_graph.R

library(gggenes)  # Arrow charts are not well covered by ggplot but worked well with this one.
library(tidyverse)
library(cowplot)

# Remake the timing plot with hands on and hands off times colored
timing <- read_tsv("TimeComparison_20240816.txt")

# 0. Should you produce the hands on time or hands off time plot ? 
HANDSOFF = "ALL"

# 0.5 Grab the possible levels
steps <- factor(timing$Step)

# 1. Tidy the data
mytimes <- timing %>%
  mutate(Step = steps) %>%
  filter(Step != "Start") %>%
  group_by(METHOD) %>%
  mutate(strand = "forward", 
         start = cumsum(duration)-duration,
         start = ifelse(start-10 > 0, start - 10, 0),
         end = cumsum(duration),
         totalTime = sum(duration)) %>%
  ungroup() %>%
  arrange(METHOD, desc(start)) %>%
  # The for loop below mess up the factor levels, easier to do it like this:
  mutate(METHOD = case_when(
    METHOD == "STORM" ~ "STORM-seq",
    METHOD == "VASA" ~ "VASA-seq",
    METHOD == "SST" ~ "Smart-seq-total",
    METHOD == "SS3E" ~ "Smart-seq3xpress",
    METHOD == "SS2" ~ "Smart-seq2",
    METHOD == "snapT" ~ "snapTotal-seq"))

mytimes$METHOD <- factor(mytimes$METHOD,
                         levels = c("STORM-seq",
                                    "Smart-seq-total",
                                    "snapTotal-seq",
                                    "VASA-seq",
                                    "Smart-seq3xpress",
                                    "Smart-seq2"
                         ))
mytimes <- mytimes %>% 
  mutate(METHOD = fct_relevel(METHOD, levels(mytimes$METHOD))) %>%
  arrange(METHOD)

# 2. Add the x-axis
if (HANDSOFF == "ALL") {
  mytimes.tpm <- mytimes
  pA <- ggplot(mytimes.tpm) +
    geom_vline(xintercept = seq(0,1920,60), linetype = "dashed", color = "darkgrey") +
    scale_x_continuous(breaks = seq(0,1920,60), labels = 0:32, "Hours") 
  
}

# 3. Create the theme of the graph
library(viridis)
step_colors <- c("black", "gray")
names(step_colors) <- unique(timing$step)
pA <- pA +theme_cowplot() +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 28),
        legend.position = "none",
        axis.text.x = element_text(size = 22), 
        axis.text.y = element_text(size = 22, face = "bold"), 
        legend.title = element_blank()) +
  scale_fill_manual(values = step_colors) +
  ylab("")

# 4. Add the blocks one-by-one as arrows using gggenes
# The loop allows the addition of each step on top of the previous, creating the white separations
for(row in 1:nrow(mytimes.tpm)){
  pA <- pA +
    geom_gene_arrow(data = mytimes.tpm[row,], aes(xmin = start, xmax = end, y = METHOD, fill = step), 
                    color = "white", 
                    arrow_body_height = unit(6, "mm"), 
                    arrowhead_height = unit(6, "mm"), 
                    arrowhead_width = unit(4, "mm"))
  
}

pA <- pA + scale_y_discrete(limits = rev) + 
  theme(legend.position = "bottom",
        legend.text = element_text(hjust = 0.5,
                                   size = 16))
pA