######### This script is for making figures from Methods ############

# source paths and libraries
source("D:/source/tk_orbitrap_dia/analysis/Utilities.R")
source("D:/projects/pca_urine_spectral_lib/src/Amanda_source.R")
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)
library(SWATH2stats)
library(aLFQ)
library(matrixStats)
library(RColorBrewer)
library(fmsb)
data.path <- "D:/projects/pca_urine_spectral_lib/data/openswath/2020_DIA_MethodOpt_ALL"
fig.path <- "D:/projects/pca_urine_spectral_lib/results/20211130_ALL_MStern_Figures/"


#################### Load data ####################################
replicate_colors <- brewer.pal(8, "Blues")[5:8]

files <- list.files(data.path, pattern = ".RData", full.names = T)
files

load(files[2], environment())
lapply(files[1:2], load, environment())
Methodopt120 <- copy(osw)
remove(osw)
Methodopt120.annotate <- copy(osw.annootate)
remove(osw.annootate)

# load the 60minute data

files
lapply(files[1:2], load, environment())
head(Methodopt60)
head(Methodopt60.annotate)
################# Part 1: Peptide Counts - 120 minute #########################
head(Methodopt120.annotate)
Methodopt120$Condition <- gsub(" 120min", "", Methodopt120$Condition)
Methodopt120.annotate$Condition <- gsub(" 120min", "", Methodopt120.annotate$Condition)
Methodopt120.annotate$Condition <- gsub(" 0[1-3]", "", Methodopt120.annotate$Condition)
Methodopt120.annotate <- Methodopt120.annotate %>%
  mutate(Range = ifelse(grepl("500mz", Condition), "500", "600"),
         WindowSize = str_extract(Condition, "[0-9]+"),
         Placement = ifelse(grepl("Variable", Condition), "variable", NA)) %>%
  mutate(Placement = ifelse(grepl("stagg", Condition), "staggered", Placement)) %>%
  mutate(Placement = ifelse(grepl("static|fixed|stepped", Condition), "fixed", Placement),
         NCE = ifelse(grepl("stepped", Condition), "stepped", "27")) %>%
  group_by(Condition) %>%
  mutate(ID = cur_group_id()) %>% as.data.table()

Methodopt120.annotate$WindowSize <- gsub("25", "24", Methodopt120.annotate$WindowSize)

#save(Methodopt120.annotate, file = files[1])

head(Methodopt120.annotate)

head(Methodopt120)
Methodopt120 <- merge(Methodopt120, Methodopt120.annotate, by = c("filename", "Condition", "Replicate"))
Methodopt120 %>%
  group_by(Condition, Replicate) %>%
  distinct(Sequence) %>% 
  summarise(PeptideCount = n())
# 120min: fixed windows comparison
Methodopt120.annotate %>% filter(Placement == "fixed")
g.staggeredRange <- Methodopt120 %>%
  filter(Placement == "staggered") %>%
  group_by(WindowSize,Range, Replicate, ID) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>%
  group_by(WindowSize, Range) %>% 
  mutate(MedianCount = median(PeptideCount)) %>%
  mutate(WindowSize = paste(WindowSize, "m/z", sep = " ")) %>%
  mutate(Range = paste(Range, "m/z", sep = " "))%>%
  ggplot() +
  geom_col(aes(x = WindowSize, y = MedianCount), fill = NA, colour = "darkgrey", size = 2, width = 0.8) +
  geom_point(aes(x = WindowSize, y = PeptideCount)) +
  # geom_signif(aes(x = WindowSize, y = PeptideCount), test = "t.test",
  #             comparisons=list(c("14 m/z", "24 m/z")), 
  #             map_signif_level = function(p) sprintf("p = %.2g", p),
  #             y_position = 14200, tip_length = 0, vjust=0.2, 
  #             size = 1.2, textsize = 5 ) +
  scale_y_continuous(limits = c(0, 16500), expand = c(0,0)) +
  ggtitle("Staggered window sizes") + ylab("Peptide counts") + xlab("") +
  plot_theme()
g.staggeredRange
g.nce
g.width
g.placement  
g.staggeredSize
ggarrange(g.width,g.nce + theme(axis.text.y = element_blank(), axis.title = element_blank()),
          g.placement + theme(axis.text.y = element_blank())
          + rremove("ylab") + rremove("y.axis") ,
          g.staggeredSize, nrow=2, ncol = 3, label.y = NA)

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staticWindowWidth_14_24_Count.svg"))

g <- Methodopt120 %>%
  filter(Placement == 'fixed' & WindowSize == "14") %>%
  group_by(NCE, Replicate, ID) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>%
  group_by(NCE) %>% 
  mutate(MedianCount = median(PeptideCount)) %>%
  ggplot() +
  geom_col(aes(x = NCE, y = MedianCount), fill = NA, colour = "darkgrey", size = 2, width = 0.8) +
  geom_point(aes(x = NCE, y = PeptideCount)) +
  scale_y_continuous(limits = c(0, 15500), expand = c(0,0)) +
  ggtitle("Collision energy") + ylab("Peptide counts") + xlab("") +
  plot_theme()

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staticNCE_14mz_Count.svg"))

g <- Methodopt120 %>%
  filter(ID %in% c(1,3,10)) %>%
  group_by(Placement, Replicate, ID) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>%
  group_by(Placement) %>% 
  mutate(MedianCount = median(PeptideCount)) %>%
  ggplot() +
  geom_col(aes(x = Placement, y = MedianCount), fill = NA, colour = "darkgrey", size = 2, width = 0.8) +
  geom_point(aes(x = Placement, y = PeptideCount)) +
  scale_y_continuous(limits = c(0, 15500), expand = c(0,0)) +
  ggtitle("Window placement") + ylab("Peptide counts") + xlab("") +
  plot_theme()

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_smallWindow_Placement_Count.svg"), width = 7)

g.staggeredSize <- Methodopt120 %>%
  filter(Placement == 'staggered' & Range == "600") %>%
  group_by(WindowSize, Replicate, ID) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>%
  group_by(WindowSize) %>% 
  mutate(MedianCount = median(PeptideCount)) %>%
  mutate(WindowSize = paste(WindowSize, "m/z", sep = " ")) %>%
  ggplot() +
  geom_col(aes(x = WindowSize, y = MedianCount), fill = NA, colour = "darkgrey", size = 2, width = 0.8) +
  geom_point(aes(x = WindowSize, y = PeptideCount)) +
  # geom_signif(aes(x = WindowSize, y = PeptideCount), test = "t.test",
  #             comparisons=list(c("14 m/z", "24 m/z")), 
  #             map_signif_level = function(p) sprintf("p = %.2g", p),
  #             y_position = 14200, tip_length = 0, vjust=0.2, 
  #             size = 1.2, textsize = 5 ) +
  scale_y_continuous(limits = c(0, 16500), expand = c(0,0)) +
  ggtitle("Window size - staggered") + ylab("Peptide counts") + xlab("") +
  plot_theme()

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staggered_size_Count.svg"), width = 7)


g.range <- Methodopt120 %>%
  mutate(WindowSize = ifelse(WindowSize == "25", "24", WindowSize)) %>%
  filter(Placement == 'staggered') %>%
  group_by(WindowSize, Range, Replicate, ID) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>%
  group_by(WindowSize, Range) %>% 
  mutate(MedianCount = median(PeptideCount)) %>%
  mutate(WindowSize = paste(WindowSize, "m/z", sep = " "),
         Range = paste(Range, "m/z", sep = " ")) %>%
  ggplot() +
  geom_col(aes(x = Range, y = MedianCount, colour = WindowSize),
           fill = NA, size = 1.5, width = 0.8, position = position_dodge(width = 0.9)) +
  geom_point(aes(x = Range, y = PeptideCount, color = WindowSize),
             position = position_dodge(width = 0.9), size = 2) +
  # geom_signif(aes(x = WindowSize, y = PeptideCount), test = "t.test",
  #             comparisons=list(c("14 m/z", "24 m/z")), 
  #             map_signif_level = function(p) sprintf("p = %.2g", p),
  #             y_position = 14200, tip_length = 0, vjust=0.2, 
  #             size = 1.2, textsize = 5 ) +
  scale_y_continuous(limits = c(0, 16500), expand = c(0,0)) +
  ggtitle("Mass range") + ylab("Peptide counts") + xlab("") +
  scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  plot_theme()

g.range

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staggered_RangeCount.svg"), width = 8)


# all peptide counts

g <- Methodopt120 %>% 
  group_by(ID, Replicate) %>%
  distinct(Sequence) %>% 
  summarise(PeptideCount = n()) %>%
  group_by(ID) %>%
  mutate(MedianCount = median(PeptideCount))%>%
  ggplot() +
  geom_col(aes(x = ID, y = MedianCount), width = 0.8, fill = NA, color = "darkgrey", size = 1.5 ) +
  geom_point(aes(x = ID, y = PeptideCount)) +
  scale_y_continuous(limits = c(0, 15500), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1,10)) +
  ylab("Peptide counts") + xlab("Method ID") + plot_theme()

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_allCount.svg"), width = 12)

ggarrange(g.width,g.nce + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
          g.placement + theme(axis.text.y = element_blank())
          + rremove("ylab") + rremove("y.axis") ,
          g.staggeredSize,
          g.range + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
          nrow=2, ncol = 3, align = "hv")
library(cowplot)
library(gridExtra)
g <- arrangeGrob(g.width,g.nce + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
            g.placement + theme(axis.text.y = element_blank())
            + rremove("ylab") + rremove("y.axis") ,
            g.staggeredSize,
            g.range + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
            ncol = 6, nrow = 2, layout_matrix = rbind(c(1,1,2,2,3,3), c(4,4,4,5,5,5)))
g
as_ggplot(g) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 32,
                  x = c(0, 0.3, 0.6, 0, 0.48), y = c(1, 1,1,0.5, 0.5)) 

# plot peptide counts for 60 minutes
Methodopt60 %>% group_by(Replicate,Condition, WindowSize, MassRange, ID) %>%
  distinct(Sequence) %>% summarise(Count = n()) %>%
  group_by(ID) %>% mutate(Median = median(Count)) %>%
  ggplot() +
  geom_col(aes(x = ID, y = Median), fill = "lightgrey", color = "black", size = 1.2, width = 0.8) +
  geom_point(aes(x = ID, y = Count)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 14500)) +
  plot_theme()
# window size
Methodopt60 %>%
  filter(MassRange == 600) %>%
  group_by(Replicate, Condition, ID, WindowSize) %>%
  distinct(Sequence ) %>% summarise(Count = n()) %>%
  mutate(WindowSize = paste(WindowSize, " m/z")) %>%
  group_by(ID) %>% mutate(Median = median(Count), WindowSize = factor(WindowSize)) %>%
  ggplot() + geom_col(aes(x = WindowSize, y = Median/1000), fill = NA, color = 'black', size = 1.2, width = 0.8) +
  geom_point(aes(x = WindowSize, y = Count/1000), size = 3) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 14.500), breaks = seq(0, 14.5, 2)) +
  plot_theme() + theme(axis.title.x = element_blank()) +
  ylab(expression("Peptides (n x "*10^3*")")) + ggtitle("Window size")

# mass range
Methodopt60 %>%
  #filter(WindowSize == "16") %>%
  mutate(MassRange = paste(MassRange, "m/z"),
         WindowSize = paste(WindowSize, "m/z")) %>%
  group_by(Replicate, Condition,MassRange, ID, WindowSize) %>%
  distinct(Sequence) %>% summarise(Count = n()/1000) %>%
  group_by( MassRange, WindowSize) %>%
  mutate(Median = median(Count)) %>%
  ggplot() + geom_col(aes(x = MassRange, y = Median, color = WindowSize),
                      position = position_dodge(width = 0.9), fill = NA, width = 0.8, size = 1.2) +
  geom_point(aes(x = MassRange, y = Count, color = WindowSize), position = position_dodge(width = 0.9), size = 3) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 14.5), breaks = seq(0, 14.5, 2)) +
  plot_theme() + theme(axis.title.x = element_blank()) +
  scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  ggtitle("Mass range") +
  ylab(expression("Peptides (n x "*10^3*")"))


################## Part 1: Peptide Consistency - 120 minute ###################
g.SizeConsistency <- Methodopt120 %>%
  filter(Placement == "fixed" & NCE == "27") %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(WindowSize, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(WindowSize, Sequence, Total) %>%
  summarise(Observation = n()) %>%
  group_by(WindowSize, Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3))) %>%
  mutate(WindowSize = paste(WindowSize, " m/z", "")) %>%
  # filter(Observation == 3) %>%
  ggplot(aes(x = WindowSize, y = Percent, colour = Observation)) +
  # geom_col(fill = NA, color = 'darkgrey', size = 1.5, width = 0.8) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
  xlab("") + ylab("Peptide detection") +
  ggtitle("Window size") +
  plot_theme()
g.SizeConsistency
save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_fixedWindow_size_countInAllTrip.svg"))


g.NCEconsistency <- Methodopt120 %>%
  filter(Placement == "fixed" & WindowSize == "14") %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(NCE, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(NCE, Sequence, Total) %>%
  summarise(Observation = n()) %>%
  group_by(NCE, Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3))) %>%
  #filter(Observation == 3) %>%
  #ggplot(aes(x = NCE, y = PeptideCount)) +
  #geom_col(fill = NA, color = 'darkgrey', size = 1.2, width = 0.8) +
  ggplot(aes(x = NCE, y = Percent, colour = Observation)) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 15000)) +
  xlab("") + ylab("Peptide detection") +
  ggtitle("Collision energy") +
  plot_theme()
g.NCEconsistency
save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_14mz_NCE_countAllReps.svg"))



g.placementConsisten <- Methodopt120 %>%
  filter(ID %in% c(1,3,10)) %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(Placement, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(Placement, Sequence, Total) %>%
  summarise(Observation = n()) %>%
  group_by(Placement, Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3))) %>%
  # filter(Observation == 3) %>%
  # ggplot(aes(x = Placement, y = PeptideCount)) +
  # geom_col(fill = NA, color = 'darkgrey', size = 1.2, width = 0.8) +
  ggplot(aes(x = Placement, y = Percent, colour = Observation)) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 12000)) +
  xlab("") + ylab("Peptide detection") +
  ggtitle("Window placement") +
  plot_theme()
g.placementConsisten
save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_smallWin_Placements_CountinAllRep.svg"), width = 7)
save_svg(as_ggplot(get_legend(g)), file = paste(fig.path, "Replicates_percent_legend.svg"), width = 4, height = 4)

g.staggeredSizeCons <- Methodopt120 %>%
  filter(Placement == "staggered" & Range == '600') %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(WindowSize, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(WindowSize, Sequence, Total) %>%
  summarise(Observation = n()) %>%
  group_by(WindowSize, Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3)),
         WindowSize = paste(WindowSize, " m/z", sep = "")) %>%
  # filter(Observation == 3) %>%
  # ggplot(aes(x = WindowSize, y = PeptideCount)) +
  # geom_col(fill = NA, color = 'darkgrey', size = 1.2, width = 0.8) +
  ggplot(aes(x = WindowSize, y = Percent, colour = Observation)) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 12000)) +
  xlab("") + ylab("Peptide detection") +
  ggtitle("Window size - staggered") +
  plot_theme()
g.staggeredSizeCons
save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staggeredSize_CountAll.svg"), width = 7)
g.rangeCons <- Methodopt120 %>%
  mutate(WindowSize = ifelse(WindowSize == "25", "24", WindowSize)) %>%
  filter(Placement == "staggered") %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(WindowSize, Range, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(WindowSize, Range, Sequence, Total) %>%
  summarise(Observation = n()) %>%
  group_by(WindowSize, Range, Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3)),
         WindowSize = paste(WindowSize, " m/z", sep = ""),
         Range = paste(Range, "m/z", sep = " ")) %>%
  # filter(Observation == 3) %>%
  # ggplot(aes(x = Range, y = PeptideCount, color = WindowSize)) +
  # geom_col(position = position_dodge(width = 0.9), fill = NA, size = 1.2, width = 0.8) +
  ggplot(aes(x = Range, y = Percent, colour = Observation)) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  facet_grid(~WindowSize) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 12000)) + 
  #scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  xlab("") + ylab("Peptide detection") +
  ggtitle("Mass range - staggered") +
  plot_theme()
g.rangeCons

save_svg(g, file = paste(fig.path, "20220420_Fig1_120min_staggered_Range_CountAllRep.svg"), width = 8)
g1 <- arrangeGrob(g.SizeConsistency,g.NCEconsistency + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
                 g.placementConsisten + theme(axis.text.y = element_blank())
                 + rremove("ylab") + rremove("y.axis") ,
                 g.staggeredSizeCons,
                 g.rangeCons + theme(axis.text.y = element_blank(), axis.title.y = element_blank()),
                 ncol = 6, nrow = 2, layout_matrix = rbind(c(1,1,2,2,3,3), c(4,4,4,5,5,5)))
g1
as_ggplot(g1) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 32,
                  x = c(0, 0.3, 0.65, 0, 0.48), y = c(1, 1,1,0.5, 0.5)) 


Methodopt60 %>%
  mutate(MassRange = paste(MassRange, "m/z")) %>%
  #filter(MassRange == 600) %>%
  group_by(ID) %>%
  mutate(Total = length(unique(Sequence))) %>% 
  group_by(WindowSize, MassRange, Replicate) %>%
  distinct(Sequence, Total) %>%
  group_by(WindowSize, Sequence, MassRange, Total) %>%
  summarise(Observation = n()) %>%
  group_by(WindowSize,MassRange,  Observation, Total) %>%
  summarise(PeptideCount = n(),
            Percent = n() / Total)  %>% distinct_all() %>%
  mutate(Observation = factor(Observation, levels = c(1,2, 3))) %>%
  mutate(WindowSize = paste(WindowSize, " m/z", "")) %>%
  # filter(Observation == 3) %>%
  ggplot(aes(x = MassRange, y = Percent, colour = Observation)) +
  # geom_col(fill = NA, color = 'darkgrey', size = 1.5, width = 0.8) +
  geom_col(position = position_stack(), fill = NA, size = 1.2, width = 0.8, show.legend = F) +
  scale_colour_brewer(type = "seq", name = "Replicates \ndetected") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
  facet_wrap(~WindowSize) +
  ylab("Peptide detection (%)") +
  ggtitle("Mass range") + 
  plot_theme() +
  theme(strip.text = element_blank(), axis.title.x.bottom = element_blank())

############## Part 2: Quantitative Comparison - 120 minute ###################
# CV
# Cycle time
# Points per peak

Methodopt120 %>% distinct(Condition)

Methodopt120.short <- reduce_OpenSWATH_output(Methodopt120)
Methodopt120.annotate$BioReplicate <- Methodopt120.annotate$Replicate
setnames(Methodopt120.annotate, 'filename', "Filename")
Methodopt120.annotate$Run <- Methodopt120.annotate$ID
Methodopt120.short <- sample_annotation(as.data.frame(Methodopt120.short), as.data.frame(Methodopt120.annotate))

head(Methodopt120.short)
Methodopt120.short %>% distinct(run_id)

disagg <- disaggregate(Methodopt120.short)

alfq.in <- convert4aLFQ(disagg)
head(alfq.in)

peps <- PeptideInference(alfq.in, transition_topx = 6, transition_strictness = "loose", 
                         transition_summary = "sum", consensus_transitions = FALSE)


Methodopt60.short <- reduce_OpenSWATH_output(Methodopt60)
Methodopt60.annotate$BioReplicate <- Methodopt60.annotate$Replicate
setnames(Methodopt60.annotate, 'filename', "Filename")
Methodopt60.annotate$Run <- Methodopt60.annotate$ID
Methodopt60.short <- sample_annotation(as.data.frame(Methodopt60.short), as.data.frame(Methodopt60.annotate))
head(Methodopt60.short)
disagg <- disaggregate(Methodopt60.short)
head(disagg)
alfq.in <- convert4aLFQ(disagg)
peps <- PeptideInference(alfq.in, transition_topx = 6, transition_strictness = "loose",
                         transition_summary = "sum", consensus_transitions = F)
head(peps)
remove(disagg)
remove(alfq.in)

Methodopt120.pep <-  dcast(peps, peptide_sequence~run_id, drop = FALSE, fill = NaN, sum, 
                           value.var = "peptide_intensity")

Methodopt60.pep <- dcast(peps, peptide_sequence ~ run_id, drop= FALSE, fill = NaN, sum,
                         value.var = 'peptide_intensity')
head(Methodopt60.pep)

head(Methodopt120.pep)
head(Methodopt120.annotate)
Methodopt120.exp <- lapply(unique(Methodopt120.annotate$Condition), function(x){
  id <- paste0(x, "_")
  cols <- colnames(Methodopt120.pep)[grepl(id, colnames(Methodopt120.pep))]
  exp <- subset(Methodopt120.pep, select = c("peptide_sequence", cols))
  exp <- exp[rowSums(is.na(exp)) != 3,]
  exp
})

Methodopt120.exp
Methodopt60.exp <- lapply(unique(Methodopt60.annotate$Condition), function(x){
  id <- paste0(x, "_")
  cols <- colnames(Methodopt60.pep)[grepl(id, colnames(Methodopt60.pep))]
  exp <- subset(Methodopt60.pep, select = c("peptide_sequence", cols))
  exp <- exp[rowSums(is.na(exp)) != 3,]
  exp
})

head(Methodopt60.exp[[1]])

peps %>% 
  mutate(Log2Intensity = log2(peptide_intensity)) %>%
  mutate(Condition = gsub("_(.*)", "", run_id)) %>%
  inner_join(Methodopt120.annotate, by = "Condition") %>%
  ggplot(aes(x = run_id, y = Log2Intensity)) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  geom_boxplot(alpha = 0.6, aes(color = factor(ID)), outlier.shape = NA) +
  scale_color_viridis_d(option = "C", name = "Method ID.") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_blank())

getCV(Methodopt120.exp[[1]][, 2:4])

Methodopt120.CV <- lapply(1:length(Methodopt120.exp), function(x){
  df <- Methodopt120.exp[[x]]
  cvs <- getCV(df[, -1])
  cvs$median <- rowMedians(as.matrix(cvs[, 1:3]))
  cvs$peptide_sequence <- df$peptide_sequence
  cvs$Condition <- gsub("_(.*)", "", colnames(cvs)[1])
  cvs <- subset(cvs, select = c("peptide_sequence","Condition", "CV", "mean", "median"))
  cvs
})

Methodopt60.CV <- lapply(1:length(Methodopt60.exp), function(x){
  df <- Methodopt60.exp[[x]]
  cvs <- getCV(df[, -1])
  cvs$median <- rowMedians(as.matrix(cvs[, 1:3]))
  cvs$peptide_sequence <- df$peptide_sequence
  cvs$Condition <- gsub("_(.*)", "", colnames(cvs)[1])
  cvs <- subset(cvs, select = c("peptide_sequence","Condition", "CV", "mean", "median"))
  cvs
})

head(Methodopt60.CV[[2]])
Methodopt60.CV[[2]] %>% summarise(MedianCV = median(CV, na.rm = T))
  

head(Methodopt120.CV[[1]])
Methodopt120.CV <- do.call('rbind', Methodopt120.CV)
Methodopt120.CV <- as.data.table(Methodopt120.CV)
Methodopt120.annotate$Run <- Methodopt120.annotate$ID
Methodopt120.annotate$WindowSize <- gsub("25", "24", Methodopt120.annotate$WindowSize)

Methodopt60.CV <- do.call('rbind', Methodopt60.CV) %>% as.data.table()
head(Methodopt60.CV)

# plot
Methodopt120.CV %>% 
  inner_join(Methodopt120.annotate[, c("Condition", "Run", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  ggplot(aes(x = Condition, y = CV )) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_classic()
head(Methodopt60.CV)
Methodopt60.exp[[2]] %>% head()


Methodopt60.CV %>%
  inner_join(Methodopt60.annotate, by = "Condition") %>%
  ggplot(aes(x = Condition, y = CV)) +
  geom_boxplot() + theme_classic()

Methodopt60.pep %>% head()
Methodopt60.pep %>% 
  ggplot(aes(x = log2(`16mz staggered 500mz_01_2`), y = log2(`16mz staggered 500mz_02_2`))) +
  geom_point() +
  stat_cor()


Methodopt60.CV %>%
  filter(Condition == "20mz staggered") %>%
  ggplot(aes(x = log2(median), y = CV)) +
  geom_point()



# Plot fixed windows

# change comparisons!!!
Methodopt120.CV %>% 
  inner_join(Methodopt120.annotate[, c("Condition", "Run", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  filter(Placement == "staggered" & WindowSize == "16") %>% 
  mutate(WindowSize = paste(WindowSize, " m/z", sep = "")) %>%
  mutate(Range = paste(Range, " m/z", sep = "")) %>%
  ggplot(aes(x = Range, y = CV)) +
  geom_jitter(alpha = 0.1, width = 0.1) +
  geom_violin(alpha = 0.8, width = 0.8) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.3) +
  #scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  scale_y_continuous(breaks = seq(0, 200, 10), expand = c(0,0), limits = c(0, 170),
                     labels = c(0, rep("", 4), 50, rep("", 4), 100, rep("", 4), 150, rep("", 4), 200)) +
  ylab("CV (%)") + xlab("") + 
  ggtitle("Mass range - 16mz staggered") + # ! change label!!
  plot_theme()
g
g2

save_svg(g, file = paste0(fig.path, "20220421_Fig1_120min_Staggered24mzRanges_CVboxplot.svg")) # change name!!

g4 <- Methodopt120.CV %>% 
  inner_join(Methodopt120.annotate[, c("Condition", "Run", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  filter(Placement == "staggered" & WindowSize == "16") %>%
  #mutate(WindowSize = paste(WindowSize, " m/z", sep = "")) %>%
  dcast(peptide_sequence ~ Range, value.var = "median", fun.aggregate = sum) %>%
  ggplot(aes(x = log2(`600`), y = log2(`500`))) +
  geom_point(alpha = 0.6) +
  stat_cor(cor.coef.name = "rho") +
  xlab("600 m/z") + ylab("500 m/z") + # ! change axis!!
  ggtitle(expression(log[2]*"(Median Intensity)"["16mz staggered"])) +
  plot_theme()
g3
g4
save_svg(g, file = paste0(fig.path, "20220421_Fig1_120min_Staggered24mz_Ranges_CorInt.svg"), width = 7.5)

ggarrange(g, g3, g2, g4, nrow = 2, ncol = 2, widths = c(6, 7.5), heights = c(6,6))
Methodopt120.CV$Condition <- gsub(" NA", "", Methodopt120.CV$Condition)

# plot 60minute comparisons
# plot window size
Methodopt60.CV %>%
  inner_join(Methodopt60.annotate, by = "Condition") %>%
  filter(MassRange == 600) %>%
  mutate(WindowSize = paste(WindowSize, "m/z")) %>%
  ggplot() +
  geom_violin(aes(x = WindowSize, y = CV)) +
  geom_boxplot(aes(x = WindowSize, y = CV), width = 0.2, alpha = 0.6, outlier.shape = NA) +
  plot_theme() + theme(axis.title.x = element_blank()) +
  ylab("CV (%)") + ggtitle("Window size")
# plot mass range
Methodopt60.CV %>%
  inner_join(Methodopt60.annotate, by = "Condition") %>%
  filter(WindowSize == 16) %>%
  mutate(MassRange = paste(MassRange, "m/z")) %>%
  ggplot() +
  geom_violin(aes(x = MassRange, y = CV)) +
  geom_boxplot(aes(x = MassRange, y = CV), width = 0.2, alpha = 0.6, outlier.shape = NA) +
  plot_theme() + theme(axis.title.x = element_blank()) +
  ylab("CV (%)") + ggtitle(expression("Mass range"[`16mz window`]))

# plot peptide abundance comparisons
Methodopt60.CV %>%
  filter(Condition %in% c("16mz staggered", "16mz staggered 500mz")) %>%
  pivot_wider(id_cols = c("peptide_sequence"),
              values_from = "median", names_from = "Condition" ) %>%
  ggplot(aes(x = log2(`16mz staggered`), y = log2(`16mz staggered 500mz`))) +
  geom_point() +
  stat_cor() +
  xlab("600 m/z") + ylab("500 m/z") +
  scale_y_continuous(limits = c(15, 38)) +
  scale_x_continuous(limits = c(15, 38)) +
  plot_theme() + ggtitle(expression(log[2]*"(Median intensity)"[`16mz window`]))


############## Part 2: Cycle time  -120 minute ##################
cycleTime.files <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/MStern_CycleTime/120min/CycleTime/",
                              pattern = ".tsv", full.names = T)
cycleTime.files

cycleTime.filenames <- basename(cycleTime.files)
cycleTime.filenames
cycleTime.filenames <- gsub("_cycletime(.*)", "", cycleTime.filenames)
cycleTime.annotation <- data.frame(filename = cycleTime.filenames, 
                                   stringsAsFactors = F)
cycleTime.anno <- annotate_pyprophet(cycleTime.annotation)
cycleTime.anno$Condition <- gsub(" 120min", "", cycleTime.anno$Condition)
cycleTime.anno

cycleTime.df <- lapply(cycleTime.files, fread)
cycleTime.df <- lapply(seq_along(cycleTime.files), function(x){
  df <- cycleTime.df[[x]]
  df[, Condition := cycleTime.anno$Condition[x]]
  df[, BioReplicate := cycleTime.anno$Replicate[x]]
  df[, filename := cycleTime.anno$filename[x]]
  df
})
cycleTime.df <- do.call('rbind', cycleTime.df)
cycleTime.df <- as.data.table(cycleTime.df)
cycleTime.df
remove(cycleTime.annotation)
cycleTime.df <- cycleTime.df[CycleTime != 0]

# compare and plot
g <- cycleTime.df %>%
  inner_join(Methodopt120.annotate %>% mutate(filename = basename(filename)),
             by = 'filename') %>% 
  filter(Placement == "staggered") %>%
  mutate(WindowSize = paste(WindowSize , "m/z", sep = " ")) %>%
  ggplot(aes(x = Range, y = CycleTime)) +
  geom_boxplot(aes(color = WindowSize), width = 0.6, show.legend = T) +
  #scale_color_manual(values = replicate_colors) + 
  scale_y_continuous(limits = c(0, 4), expand = c(0,0), breaks = seq(0, 4, 0.5)) +
  ggtitle("Window sizes") + ylab("Cycle time (s)") +
  scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  plot_theme() +
  theme(axis.title.x = element_blank())
g
save_svg(g, file = paste0(fig.path, "20220425_Fig1_120min_staggered_RangesCTBoxplot.svg"), width = 7.5)
g <- get_legend(g)
save_svg(as_ggplot(g), file = paste0(fig.path, "20220425_Fig1_120min_staggered_RangesCTBoxplotLegend.svg"), width = 3, height = 3)
head(cycleTime.df)

cycleTime.files <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/MStern_CycleTime/60min/CycleTime/",
                              pattern = ".tsv", full.names = T)
basename(cycleTime.files)
cycleTime.files <- cycleTime.files[3:20]
cycleTime.files

cycleTime.60min <- lapply(cycleTime.files, fread)
head(cycleTime.60min[[2]])
cycleTime.60min <- lapply(1:length(cycleTime.60min), function(x){
  f <- cycleTime.files[x]
  file <- basename(f)
  file <- str_split(file, "_")[[1]]
  replicate <- gsub("\\.(.*)", "",file[grepl("mzML", file)])
  cycleTime.60min[[x]]$Condition <- paste(file[c(2,3,4)], collapse = " ")
  cycleTime.60min[[x]]$Replicate <- replicate
  cycleTime.60min[[x]]
})

cycleTime.60min <- do.call('rbind', cycleTime.60min) %>% as.data.table()
cycleTime.60min$Condition <- gsub(" 60min", "", cycleTime.60min$Condition)
cycleTime.60min %>% group_by(Condition, Replicate) %>% summarise(MedianCT = median(CycleTime))

cycleTime.60min %>%
  mutate(Condition = gsub("staggererd", "staggered", Condition)) %>%
  inner_join(Methodopt60.annotate, by = c("Condition", "Replicate")) %>%
  filter(MassRange == 600) %>%
  group_by(Condition, WindowSize) %>%
  summarise(MedianCT = median(CycleTime)) %>%
  mutate(WindowSize = paste(WindowSize, "m/z")) %>%
  ggplot(aes(x = WindowSize, y = MedianCT)) +
  geom_col(width = 0.8, fill = NA, color = "black", size = 1.2) +
  scale_y_continuous(limits = c(0, 3.50), breaks = seq(0, 3.5, 0.5), expand = c(0,0)) +
  plot_theme() +
  theme(axis.title.x = element_blank()) +
  ylab("Median cycle time (s)") +
  ggtitle("Window size")

cycleTime.60min %>%
  mutate(Condition = gsub("staggererd", "staggered", Condition)) %>%
  inner_join(Methodopt60.annotate, by = c("Condition", "Replicate")) %>%
  filter(MassRange == 600) %>%
  group_by(Condition, WindowSize) %>%
  summarise(MedianCT = median(CycleTime)) %>%
  mutate(WindowSize = paste(WindowSize, "m/z")) %>%
  ggplot(aes(x = WindowSize, y = MedianCT)) +
  geom_col(width = 0.8, fill = NA, color = "black", size = 1.2) +
  scale_y_continuous(limits = c(0, 3.50), breaks = seq(0, 3.5, 0.5), expand = c(0,0)) +
  plot_theme() +
  theme(axis.title.x = element_blank()) +
  ylab("Median cycle time (s)") +
  ggtitle("Window size")


############## Part 3: Points / peak - 120 minute #############3
chromPoints.files <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/MStern_CycleTime/120min/ChromPoints/",
                                pattern = ".tsv", full.names = T)
chromPoints.files <- list.files("D:/projects/pca_urine_spectral_lib/data/openswath/MStern_CycleTime/60min/ChromPoints/",
                                pattern = ".tsv", full.names = T)

chromPoints.60min <- lapply(chromPoints.files, fread)
head(chromPoints.60min[[1]])
chromPoints.files

chromPoints.files
chromPoints.filenames <- basename(chromPoints.files)
chromPoints.filenames
chromPoints.annotation <- data.frame(filename = chromPoints.filenames, 
                                     stringsAsFactors = F)
chromPoints.annotation$filename <- gsub("_ChromPoints(.*)", ".mzML.gz", 
                                        chromPoints.annotation$filename)
chromPoints.annotation <- annotate_pyprophet(chromPoints.annotation)
chromPoints.annotation$Condition <- gsub(" 120min", "", 
                                         chromPoints.annotation$Condition)
chromPoints.annotation
chromPoints.df <- lapply(seq_along(chromPoints.files), function(x){
  df <- fread(chromPoints.files[x])
  df[, Condition := chromPoints.annotation$Condition[x]]
  df[, BioReplicate := chromPoints.annotation$Replicate[x]]
  df[, filename := chromPoints.annotation$filename[x]]
  df
})

chromPoints.df <- do.call('rbind', chromPoints.df) %>% as.data.table()
head(chromPoints.df)

chromPoints.df.peptide
chromPoints.df.peptide <- chromPoints.df %>% 
  group_by(Condition, BioReplicate, Sequence, Charge, filename) %>%
  slice_max(order_by = Points, n = 6, with_ties = F) %>% # take only the best 6 transitions
  summarise(MeanPoints = mean(Points),
            MaxPoints = max(Points),
            Transitions = n()) %>%
  group_by(Condition, BioReplicate, Sequence, filename) %>%
  arrange(desc(MeanPoints)) %>%
  distinct(Condition, BioReplicate, Sequence, .keep_all = T) %>%
  mutate(BioReplicate = as.integer(BioReplicate)) %>% ungroup()

chromPoints.df.peptide <- as.data.table(chromPoints.df.peptide)
g <- chromPoints.df.peptide %>%
  inner_join(Methodopt120.annotate %>% 
               mutate(filename = basename(filename)), by = "filename") %>%
  filter(Placement == "staggered") %>%
  group_by(WindowSize, Range , Replicate) %>%
  mutate(MedianPt = median(MeanPoints)) %>%
  group_by(WindowSize, Range) %>%
  mutate(Median = median(MeanPoints)) %>% distinct(WindowSize,Range,  MedianPt, Median, Replicate) %>%
  mutate(WindowSize = paste(WindowSize, "m/z", sep = " "),
         Range = paste(Range, "m/z", sep = " ")) %>%
  ggplot() +
  geom_col(aes(x = Range, y = Median, color = WindowSize), width = 0.8, fill = NA, size = 2, 
           position = position_dodge(width = 0.9), show.legend = T) +
  geom_point(aes(x = Range, y = MedianPt, color = WindowSize),
             position = position_dodge(width = 0.9), size = 2.5, show.legend = T) +
  scale_y_continuous(limits = c(0, 11), expand = c(0,0)) +
  scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  ylab("Points per peak") + 
  ggtitle("Mass range") +
  plot_theme() +
  theme(axis.title.x = element_blank())

g.ChromPoints

save_svg(g, file = paste0(fig.path, "20220425_Fig1_120min_Staggered_SizesRangeChromP.svg"), width = 7.5)
save_svg(as_ggplot(get_legend(g)), 
         file = paste0(fig.path, "20220425_Fig1_120min_Staggered_SizesRangeChromPLegend.svg"),
         width = 4, height = 3)

chromPoints.60min <- lapply(1:length(chromPoints.60min), function(x){
  filename <- basename(chromPoints.files[x])
  filename <- str_split(filename, "_")[[1]]
  condition <- paste(filename[2:4], collapse = " ")
  condition <- gsub(" 60min", "", condition)
  df <- chromPoints.60min[[x]]
  df[, Condition := condition]
  df$Replicate <- filename[grepl("^0", filename)]
  df
})

head(chromPoints.60min[[1]])

chromPoints.60min <- do.call('rbind', chromPoints.60min) %>% as.data.table()
chromPoints.60min.peptide <- chromPoints.60min %>% 
  group_by(Condition, Sequence, Charge, Replicate) %>%
  slice_max(order_by = Points, n = 6, with_ties = F) %>% # take only the best 6 transitions
  summarise(MeanPoints = mean(Points),
            MaxPoints = max(Points),
            Transitions = n()) %>%
  group_by(Condition, Replicate) %>%
  arrange(desc(MeanPoints)) %>%
  distinct(Condition, Sequence, Replicate, .keep_all = T)

chromPoints.60min.peptide$Condition <- gsub("staggererd", "staggered", chromPoints.60min.peptide$Condition)

setnames(chromPoints.60min.peptide, "Replicate", "BioReplicate")
chromPoints.60min.peptide %>%
  inner_join(Methodopt60.annotate, by = c("Condition", "BioReplicate")) %>%
  mutate(ID = factor(ID)) %>%
  ggplot(aes(x = ID, y = MeanPoints,color = BioReplicate)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(position = position_dodge(width= 0.9), width = 0.3) +
  theme_classic()

chromPoints.60min.peptide %>%
  inner_join(Methodopt60.annotate, by = c("Condition", "BioReplicate")) %>%
  group_by(Sequence, ID, WindowSize, MassRange) %>%
  summarise(MedianPt = median(MeanPoints)) %>%
  filter(MassRange == 600) %>%
  mutate(WindowSize = paste(WindowSize, "m/z")) %>%
  ggplot(aes(x = WindowSize, y = MedianPt)) +
  geom_boxplot() +
  plot_theme() +
  theme(axis.title.x = element_blank()) +
  ggtitle("Window size") + ylab("Points per peak")

############### Part 3: Summary plots ####################

# outcomes: median peptide count, median points per peak, mean cycle time, median CV
# 
library(fmsb)
cv.sum <- Methodopt120.CV %>% 
  inner_join(Methodopt120.annotate[, c("Condition", "Run", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  group_by(ID, NCE, WindowSize, Placement, Range, Condition) %>%
  summarise(medianCV = median(CV, na.rm = T)) %>% ungroup()

chromPoints.sum <- chromPoints.df.peptide %>% 
  inner_join(Methodopt120.annotate %>% 
               mutate(filename = basename(filename)), by = "filename") %>%
  group_by(ID, NCE, WindowSize, Placement, Range) %>% 
  summarise(MedianPt = median(MeanPoints)) %>% ungroup()

ct.sum <- cycleTime.df %>%
  inner_join(Methodopt120.annotate %>% mutate(filename = basename(filename)),
             by = 'filename') %>% 
  group_by(ID, NCE, WindowSize, Placement, Range) %>%
  summarise(MedianCT = median(CycleTime)) %>% ungroup()

counts.sum <- Methodopt120 %>%
  mutate(WindowSize = gsub("25", "24", WindowSize) ) %>%
  group_by(WindowSize, Replicate, ID, Placement, Range, NCE, Condition) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()) %>% 
  group_by(ID, WindowSize, NCE, Placement, Range) %>%
  summarise(medianCount = median(PeptideCount)) %>% ungroup()

sumDF <- Reduce(function(x, y) merge(x,y, by = c("ID", "NCE", "WindowSize", "Placement", "Range"), all = T),
                list(cv.sum, chromPoints.sum, ct.sum, counts.sum))

sumDF$medianCount <- as.numeric(sumDF$medianCount)
sumDF$medianCV <- as.numeric(sumDF$medianCV)
sumDF$MedianCT <- as.numeric(sumDF$MedianCT)
sumDF$MedianPt <- as.numeric(sumDF$MedianPt)

str(sumDF)
?radarchart
radarchart(df, axistype = 0)

sumDF[11,] <- c("MAX", NA, NA, NA, NA, "MAX", 20, 11, 3.6, 16000)
sumDF[12, ] <- c("MIN", NA, NA, NA, NA, "MIN", 8, 6, 1.5, 6500)


radarchart(df)
radarchart(sumDF[c(1,11, 12), c('medianCV', 'MedianPt', "MedianCT", "medianCount")])

sumDF <- as.data.table(sumDF)
sumDF <- sumDF[order(-ID)]
sumDF
sumDF[1:2, ] <- sumDF[2:1, ]
sumDF$medianCV <- as.numeric(sumDF$medianCV)
sumDF$medianCount[2] <- 12000
library(BoutrosLab.plotting.general)

show.available.palettes()
color <- default.colours(number.of.colours = 12, palette.type = "pastel")
display.colours(color)
# compare window placements
# compare mass range
# REVERSE ORDER
color <- color[-c(1,7)]
sumDF$color <- c(NA, NA, rev(color))
radarchart(sumDF[c(1,2,4,7,9), 7:10], 
           pcol = sumDF$color[c(4,7,9)], axistype = 0,
           caxislabels = NA,
           pfcol = scales::alpha(sumDF$color[c(4,7,9)], 0.3), plwd = 5, plty = 1,
           # Customize the grid
           cglcol = "darkgrey", cglty = 1, cglwd = 4,
           axislabcol ="darkgrey",
           vlabels = c("CV", "Points \nper peak", "Cycle time", "Peptide\ncount"),
           vlcex = 5) 

#display.colours(sumDF$color[c(9, 11, 12)])

# Customize the axis)

legend(x = -1.55, y = -0.45,
  legend = paste(sumDF$WindowSize[c(4,7,9)], "m/z"), horiz = F,
  bty = "n", pch = 20 , col = sumDF$color[c(4,7,9)],
  text.col = "black", cex = 2, pt.cex = 2.2
)
l

print(l)
grid.newpage()
grid.draw(l)

##################### part 4: summarise barplot ##################3

bestOpt <- data.frame(Method = c("DIA-45min", "DDA"), PeptideCount = c(12882, 6331), stringsAsFactors = F)
bestOpt %>% 
  mutate(yLabel = PeptideCount / 1000) %>%
  ggplot(aes(x = Method, y = yLabel)) +
  geom_col(aes(fill = Method, color = Method), width = 0.6, size = 2, show.legend = F) +
  scale_color_manual(breaks = c("DDA", "DIA-45min"), values = c("#FC8D62","#66C2A5")) +
  scale_fill_manual(breaks = c("DDA", "DIA-45min"), values = c("#FC8D624C","#66C2A54C")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 14), breaks = seq(0, 144, 2.500)) +
  plot_theme() + ylab(expression("Peptides (n) x"*10^3))

sumDF %>% filter(grepl("[0-9]", ID)) %>% 
  mutate(ID = as.integer(ID)) %>%
  ggplot(aes(x = ID, y = medianCount)) + 
  geom_col(color = "black", fill = "grey", width = 0.8, size = 1.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 16500)) + 
  scale_x_reverse(breaks = seq(1, 10)) +
  theme_light(base_size = 24) +
  coord_flip() + ylab("Peptides (n)")

Methodopt120.CV %>% 
  inner_join(Methodopt120.annotate[, c("Condition", "Run", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  ggplot(aes(x = ID, y = CV)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.8, aes(group = ID)) +
  #scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  scale_y_continuous(breaks = seq(0, 200, 10), limits = c(0, 100)) +
  scale_x_reverse(breaks = seq(1,10)) +
  ylab("CV (%)") + xlab("") + 
  theme_light(base_size = 24) + coord_flip()

sumDF %>% 
  filter(grepl("[0-9]", ID)) %>% 
  mutate(ID = as.integer(ID)) %>%
  ggplot(aes(x = ID, y = MedianCT)) + 
  geom_col(color = "black", fill = "grey", width = 0.8, size = 1.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4)) + 
  scale_x_reverse(breaks = seq(1, 10)) +
  theme_light(base_size = 24) +
  coord_flip() + ylab("Peptides (n)")

sumDF %>% filter(grepl("[0-9]", ID)) %>% 
  mutate(ID = as.integer(ID)) %>%
  ggplot(aes(x = ID, y = MedianPt)) + 
  geom_point(color = default.colours(1, palette.type = "spiral.dusk"), fill = "grey", size = 4) +
  geom_segment(aes(yend = MedianPt, xend = ID, y = 0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12)) + 
  scale_x_reverse(breaks = seq(1, 10)) +
  theme_light(base_size = 24) +
  coord_flip() + ylab("Points \n per peak")




g


df2 <- sumDF %>%
  filter(grepl("[0-9]", ID)) %>% 
  mutate(ID = as.integer(ID),
         WindowSize = as.numeric(WindowSize),
         medianCV = as.numeric(medianCV),
         Range = as.numeric(Range)) %>%
  select(ID, WindowSize, Range, medianCV, MedianPt, MedianCT, medianCount) %>%
  pivot_longer(
    cols = c(WindowSize, Range, medianCV, MedianPt, MedianCT, medianCount),
    names_to = "Parameters",
    values_to = "value"
  )
head(df2)

ggdotchart(df2, x = "ID", y = "value", 
           group = "Parameters", color = "Parameters", palette = "jco",
           add = "segment", position = position_dodge(0.3),
           sorting = "descending", facet.by = "Parameters",
           rotate = TRUE, legend = "none", scale = "free_x")

g.count <- Methodopt120 %>%
  group_by(ID, Replicate, Condition) %>%
  distinct(Sequence) %>%
  summarise(PeptideCount = n()/1000) %>%
  group_by(ID) %>%
  mutate(MedianCount = median(PeptideCount)) %>%
  ggplot() +
  geom_col(aes(x = ID, y = MedianCount), fill = "lightgrey", color = 'black', width = 0.8, size = 1.5) +
  geom_point(aes(x = ID, y = PeptideCount), size = 2.5) +
  scale_y_continuous(limits = c(0, 16.5), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1, 10)) +
  theme_classic(base_size = 24) +
  #theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  ylab(expression(atop("Peptides", paste("(n x", 10^{3}, ")"))))
g.count
g.CV <- Methodopt120.CV %>%
  inner_join(Methodopt120.annotate[, c("Condition", "ID", "NCE", "WindowSize","Placement", "Range")], 
             by = "Condition") %>%
  ggplot(aes(x = ID, y = CV, group = ID)) +
  #geom_violin() +
  geom_boxplot(width = 0.8, outlier.shape = NULL) +
  scale_x_reverse(breaks = seq(1, 10)) + 
  theme_classic(base_size = 24) + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  ylab("CV \n (%)") + coord_flip()

g.CT <- cycleTime.df %>%
  inner_join(Methodopt120.annotate %>% mutate(filename = basename(filename)),
             by = 'filename') %>%   
  group_by(ID) %>%
  summarise(MedianCT = median(CycleTime)) %>%
  ggplot(aes(x = ID, y = MedianCT)) +
  geom_col(fill = 'lightgrey', color = 'black', size = 1.5, width = 0.8) +
  #scale_color_manual(values = replicate_colors) + 
  scale_y_continuous(limits = c(0, 4), expand = c(0,0), breaks = seq(0, 4, 0.5), labels = c(0, "", 1, "", 2, "", 3, "", 4)) +
  ylab("Cycle \n time (s)") +
  theme_classic(base_size = 24) +
  scale_x_reverse(breaks = seq(1, 10)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + coord_flip()


g.ChromPoints <- chromPoints.df.peptide %>%
  inner_join(Methodopt120.annotate %>% 
               mutate(filename = basename(filename)), by = "filename") %>%
  ggplot() +
  geom_violin(aes(x = ID, y = MeanPoints, group = ID)) +
  geom_boxplot(aes(x = ID , y = MeanPoints, group = ID), outlier.shape = NA, width = 0.3) +
  scale_x_continuous(breaks = seq(1, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 10), 
                     labels = c(0, rep("", 4), 50, rep("", 4), 100, "", "" )) +
  #scale_color_viridis_d(option = "magma", direction = 1, end = 0.8, name = "Window width") +
  ylab("Points \n per peak") + 
  theme_classic(base_size = 24) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + coord_flip()

g.count
g.CV
g.CT
g.ChromPoints

ggarrange(g.Windows, g.count, g.CV, g.CT, g.ChromPoints, nrow = 1, align = "h")

sumDF
sumDF$Windows <- c(40, 31, 40, 38, 32, 31, 26, 25, 25, 22)
g.Windows <- sumDF %>%
  ggplot(aes(x = ID, y = Windows)) +
  geom_col(fill = "#65B4A2", color = "black", width = 0.8, size = 1.2) +
  scale_x_reverse(breaks = seq(1, 10)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 42)) +
  theme_classic(base_size = 24) + ylab("Number of \n windows") +
  coord_flip() + xlab("Method ID ")
g.Windows

# make covariate bars 
# covariate for window placement (variable, staggered, fixed), maybe window size, mass range, collision energy
sumDF <- as.data.table(sumDF)
sumDF <- sumDF[order(-medianCount)]
sumDF
sumDF$Rank <- as.character(seq(1:nrow(sumDF)))
covariate.data <- data.frame(
  WindowPlacement = sumDF$Placement,
  MassRange = sumDF$Range,
  WindowSize = sumDF$WindowSize,
  CollisionEnergy = sumDF$NCE)
covariate.numeric <- data.frame(lapply(covariate.data, as.character), stringsAsFactors = FALSE)

covariate.numeric
covariate.numeric <- covariate.numeric %>%
  mutate(WindowPlacement = case_when(WindowPlacement == 'fixed' ~ 1,
                                     WindowPlacement =='variable' ~ 2,
                                     WindowPlacement =='staggered' ~3),
         MassRange = case_when(MassRange == '600' ~ 4, MassRange == '500' ~ 5),
         WindowSize = case_when(WindowSize =='14' ~ 6,
                                WindowSize == '16' ~ 7,
                                WindowSize == '20' ~ 8,
                                WindowSize == '24' ~ 9,
                                WindowSize == '31' ~ 10),
         CollisionEnergy = case_when(CollisionEnergy == '27' ~ 11,
                                     CollisionEnergy == 'stepped' ~ 12)) %>%
  as.data.frame

covariate.numeric <- covariate.numeric[,order(ncol(covariate.numeric):1)]

placement.colors <- c("seagreen3" ,  "violetred3",  "#154389")
range.colors <- c("#FDFCB7" ,"#87B3C4")
size.colors <- default.colours(5, 'spiral.morning')
energy.colors <- c("white", "darkgrey")
covariate.bar <- create.heatmap(
  x = data.matrix(covariate.numeric),
  clustering.method = "none",
  print.colour.key = FALSE,
  # set colour scheme
  total.colours = 13,
  colour.scheme = c(placement.colors, range.colors, size.colors, energy.colors),# add row lines
  grid.row = TRUE,
  grid.col = TRUE,
  row.colour = "black",
  col.colour = "black",
  yaxis.lab = T
)

covariate.bar

placement.leg <- list(colours = placement.colors, labels = c("fixed", "variable", "staggered"), title = "Window\nplacement")
size.leg <- list(colours = size.colors, labels = c("14", "16", "20", "24", "Multi"), title = "Window \nsize (m/z)")
range.leg <- list(colours = range.colors, labels = c("600", "500"), title = "Mass \nrange (m/z)")
energy.leg <- list(colours = energy.colors, labels = c("27", "stepped"), title = "NCE")

covariate.legends <- list(legend = placement.leg, 
                          legend = range.leg,
                          legend = size.leg,
                          legend = energy.leg)

legend1 <- legend.grob(
  covariate.legends,
  size = 1.25,
  label.cex = 0.75,
  title.cex = 0.75,
  layout = c(1,length(covariate.legends)),
  # add black box around legend
  border = list(col = "black", lwd = 3, lty = 1),
  border.padding = 1.5
)

sumDF$Rank <- as.integer(sumDF$Rank)
sumDF$Rank <- as.numeric(sumDF$Rank)
sumDF$Rank <- factor(sumDF$Rank)
g.count
bpg.count <- create.scatterplot(formula = medianCount ~ Rank, 
                                data = sumDF, 
                                type = c("p",'h'), xat = seq(1, 10))
bpg.count <- create.barplot(medianCount ~ Rank, data = sumDF, xat = seq(1, 10))
bpg.count
bpg.CT <- create.barplot(MedianCT~Rank,data = sumDF,  xat = sumDF$Rank)
bpg.CT 

sumChrom <- merge(chromPoints.df.peptide, Methodopt120.annotate %>%
                    mutate(filename = basename(filename)), by = "filename")

head(sumChrom)

sumChrom <- merge(sumChrom, sumDF, by.x = c("Range", "WindowSize", "Placement", "NCE", "ID"),
                  by.y = c("Range", "WindowSize", "Placement", "NCE", "ID"))

bpg.PointsBoxPlot <- create.boxplot(MeanPoints~factor(Rank), data = sumChrom, add.stripplot = TRUE,
               points.pch = 1 , jitter.factor = 0.5)
bpg.PointsBoxPlot
bpg.PointsVioin <- create.violinplot(MeanPoints~factor(Rank), data = sumChrom)

sumCV <- merge(Methodopt120.CV, Methodopt120.annotate %>%
                 distinct(Condition, Range,WindowSize, Placement, NCE, ID), by = "Condition")

sumCV <- merge(sumCV, sumDF, by = c("Range", "WindowSize", "Placement", "NCE", "ID"))
sumCV$Condition.x %>% unique()

bpg.CV <- create.boxplot(CV ~ factor(Rank), data = sumCV, add.stripplot = T, points.pch = 1, jitter.factor = 0.5)



# ultimately:
# 1. peptide count
# 2. Cycle time
# 3. Points per peak
# 4. CVs
# Covariate


create.multiplot(plot.objects = list(covariate.bar,bpg.CV, bpg.PointsVioin, bpg.CT, bpg.count),
                 plot.layout = c(1,5),
                 xaxis.cex = 1.3,yaxis.cex = 1.3,
                 main.x = 'Rank',
                 main.key.padding = 4,
                 ylab.cex = 1.3,
                 ylab.padding = 6,
                 xaxis.alternating = 0, yaxis.alternating = 0, 
                 yat = list(1:5, seq(0, 200 , 50), seq(0, 150, 40), seq(0, 3.8, 0.5), seq(0, 17000, 2000) ),
                 yaxis.labels = list(NULL, seq(0, 150, 50), seq(0, 120, 40), seq(0, 3.8, 0.5), seq(0, 17000, 2000)),
                 ylab.label = rev(c("\t", "CV \n(%)","\t", "Points\n per peak","\t", "Cycle \ntime (s)", "\t",
                                    'Peptide \ncount')),
                 x.relation = 'free', y.relation = 'free',
                 panel.heights = c(1, 1, 1,1, 0.3), 
                 y.spacing = c(0.5, 0.5, 0.5, 0.5),
                 print.new.legend = T,
                 legend = list(right = list(
                   fun = legend1,
                   x = 0.5,
                   y = -0.1,
                   corner = c(0.5, 0.5))),
                 filename = "D:/projects/pca_urine_spectral_lib/results/20211130_ALL_MStern_Figures/20220705_dia_120min_Method_Covariate.pdf",
                 width = 8.5,
                 height = 11)

