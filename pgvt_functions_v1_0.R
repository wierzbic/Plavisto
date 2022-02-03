# Plavisto - Plastid Genome Visualization Tool
# Version 1.0 (02.03.2022)
# Copyright (C) 2022  Andrzej Wierzbicki
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

load_data_tmp <- function(file, strand, name) #function to load data from input files into a data frame
{
  #browser()
  input_data <- fread(paste("./datasets/", file, sep = ""), sep='\t', stringsAsFactors = FALSE, drop = c(1,2,3))
  if(strand == '-') {    strand_value <- -1
  } else{                strand_value <- 1 }
  if(ncol(input_data) == 2) {                                                       #test if there are 5 columns in the input file (no sd)
    coverage <- data.frame(input_data$V4, input_data$V5, name, strand_value, 0, row.names = NULL)  #sd is set at 0
    colnames(coverage) <- c('position', 'data', 'name', 'strand', 'sd')
  } else {
    coverage <- data.frame(input_data$V4, input_data$V5, name, strand_value, input_data$V6, row.names = NULL) #sd loaded from file
    colnames(coverage) <- c('position', 'data', 'name', 'strand', 'sd')
  }
  return(coverage)
}

#Allow caching of datasets to memory, shared between users
load_data <- memoise(load_data_tmp, cache = cachem::cache_mem(max_size = 64 * 1024^2)) #Cache size limit

reselect_datasets_tree <- function(data_list_item, tree_list, coverage_plot_globallist, wc_globallist){   #change displayed datasets based on tree selection take a data frame to modify
  #browser()
  data_to_load <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE, na.strings=c("NA", ""))
  time <- Sys.time()
  existing_data <- unique(coverage_plot_globallist[[data_list_item]]$name)      #find names of all existing datasets
  kept_data <- tree_list[tree_list %in% existing_data]          #names of datasets to keep
  added_data <- tree_list[!tree_list %in% existing_data]        #names of datasets to add
  removed_data <- existing_data[!existing_data %in% tree_list]  #names of datasets to remove
  if(length(removed_data) > 0){
    coverage_plot_globallist[[data_list_item]] <- filter(coverage_plot_globallist[[data_list_item]], name %in% kept_data)  #remove deselected data
    cat("Removed datasets: ")
    for(d in removed_data) cat(paste(d, ", "))
    cat(paste("time: ", Sys.time() - time, "\n"))  
  }
  if(length(added_data) > 0){
    #browser()
    for(i in added_data){
      if(!is.na(data_to_load$filename_plus[data_to_load$samplename == i])) coverage_plot_globallist[[data_list_item]]  <- bind_rows(coverage_plot_globallist[[data_list_item]], load_data(data_to_load$filename_plus[data_to_load$samplename == i], "+", i))
      if(!is.na(data_to_load$filename_minus[data_to_load$samplename == i])) coverage_plot_globallist[[data_list_item]] <- bind_rows(coverage_plot_globallist[[data_list_item]], load_data(data_to_load$filename_minus[data_to_load$samplename == i], "-", i))
    }
    wc_globallist[[data_list_item]] <- summarize(group_by(coverage_plot_globallist[[data_list_item]], name), count = sum(data))  #replace the wc counts for the list object using dplyr functions
    cat("Added datasets: ")
    for(d in added_data) cat(paste(d, ", "))
    cat(paste("time: ", Sys.time() - time, "\n"))   
  }
  return(list(coverage_plot_globallist, wc_globallist))
}

plot_navi <- function(position_s, position_e, annotation, operons)
{
  position_df <- data.frame(cbind(position_s, position_e))
  ggplot() +
    geom_segment(aes(x = genome_regions[1], y = 0, xend = genome_regions[2], yend = 0), lineend = "round", linejoin = "mitre", size = 1.4, color = "#807dba", arrow = arrow(length = unit(0.2, "cm"))) +
    geom_segment(aes(x = genome_regions[4], y = 0, xend = genome_regions[3], yend = 0), lineend = "round", linejoin = "mitre", size = 1.4, color = "#807dba", arrow = arrow(length = unit(0.2, "cm"))) +
    
    geom_rect(data=annotation[annotation$V4 == "r", ], aes(xmin=V1, xmax=V2, ymin=0.5*V8, ymax=V8), fill = "#FC8D62") +
    geom_rect(data=annotation[annotation$V4 == "t", ], aes(xmin=V1, xmax=V2, ymin=0.5*V8, ymax=V8), fill = "#E78AC3") +
    geom_rect(data=annotation[annotation$V4 == "p", ], aes(xmin=V1, xmax=V2, ymin=0.5*V8, ymax=V8), fill = "#66C2A5") +
    geom_rect(data=annotation[annotation$V4 == "o", ], aes(xmin=V1, xmax=V2, ymin=0.5*V8, ymax=V8), fill = "#FFD92F") +
    geom_rect(data=operons, aes(xmin=V1, xmax=V2, ymin=0.25*V7, ymax=0.25*V7), color = "#CC79A7", size=1) +
    geom_rect(aes(xmin=position_s, xmax=position_e, ymin=-1, ymax=1), fill="grey20", alpha=0.5, color="black") +
    scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
    scale_x_continuous(limits = c(0, 154477), expand=c(0,0)) +
    theme_void() +
    theme(text = element_text(size=15), legend.position="none", panel.background = element_rect(fill = 'white', colour = 'white')) 
}

plot_MNase <- function(data_list_item, seq_start, seq_end, bin_size, zoom, log_transform, normalization, click_equalize, equalize_position, click_yvalue, differential_display, coverage_plot_globallist, wc_globallist, annotation, operons)
{
  #Response to empty dataset
  if(is.null(coverage_plot_globallist[[data_list_item]])) {                                   
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste0("Use the Datasets panel to select data for plot #", data_list_item), cex = 1.6, col = "black")
    return(0)
  }
  #browser()
  #select and scale down a subset of the dataset for further manipulations and display using dplyr
  coverage_plot_tmp <- filter(coverage_plot_globallist[[data_list_item]], position >= seq_start & position < seq_end)        #select a subset
  coverage_plot_tmp$name <- factor(coverage_plot_tmp$name,  levels = unique(coverage_plot_tmp$name)) #correct the list of levels to maintain proper order later
  #coverage_plot_tmp <- group_by(coverage_plot_tmp, name, position)                              #combine data from BOTH strands
  #coverage_plot_tmp <- summarize(coverage_plot_tmp, data = sum(data), sd=sum(sd))               #
  if(seq_end - seq_start > bin_size){                                   #test if number of data points is greater than max number of bins
    coverage_plot_tmp <- mutate(coverage_plot_tmp, window = ntile(position, n = bin_size))     #add column with bin numbers
    coverage_plot_tmp <- group_by(coverage_plot_tmp, name, window)                             #subset based on name and bin number
    coverage_plot_lim <- summarize(coverage_plot_tmp, enrichment = mean(data), sd = mean(sd), position = mean(position)) #calculate sum value and assign to min position
  } else {
    coverage_plot_lim <- coverage_plot_tmp
    coverage_plot_lim$enrichment <- coverage_plot_lim$data
  }
  
  #data normalization
  if(normalization == 0){        #normalization based on total read counts
    dataset1_read_count <- as.numeric(unlist(wc_globallist[[data_list_item]][1, 2]))                          #Read count of dataset #1
    coverage_plot_lim <- left_join(coverage_plot_lim, wc_globallist[[data_list_item]], by = "name")           #For each data point get the normalozation value based on name match
    coverage_plot_lim$normalized <- coverage_plot_lim$enrichment / coverage_plot_lim$count * dataset1_read_count    #calculate normalized signal levels
    coverage_plot_lim$sd <- coverage_plot_lim$sd / coverage_plot_lim$count * dataset1_read_count    #calculate normalized sd values
  } else if(normalization == 1){        #normalization based on local read counts
    wc2 <- sapply(unique(coverage_plot_lim$name), function(x) sum(coverage_plot_lim$enrichment[coverage_plot_lim$name == x]))
    coverage_plot_lim$normalized <- coverage_plot_lim$enrichment / wc2[coverage_plot_lim$name] * wc2[wc2 > 0][1]    #normalize all datasets to the first dataset using total read counts
    coverage_plot_lim$sd <- coverage_plot_lim$sd / wc2[coverage_plot_lim$name] * wc2[wc2 > 0][1]    #normalize all sd to the first dataset using total read counts
  } else if(normalization == 3){        #no normalization
    coverage_plot_lim$normalized <- coverage_plot_lim$enrichment
  }
  #display differential data relative to the first dataset
  if(differential_display != 0){
    datasets_list <- unique(coverage_plot_lim$name)     #get list of unique sample names
    coverage_plot_ref <- filter(coverage_plot_lim, name == datasets_list[1])  #gather reference data, first dataset, selected strand
    coverage_plot_lim <- left_join(coverage_plot_lim, coverage_plot_ref, by = "position")                         #match reference value with each position in coverage_plot_lim
    coverage_plot_lim$normalized <- (coverage_plot_lim$normalized.x + 1) / (coverage_plot_lim$normalized.y + 1)   #divide datapoints by reference data 
    coverage_plot_lim$sd <- coverage_plot_lim$normalized * sqrt( (coverage_plot_lim$sd.x / coverage_plot_lim$normalized.x)^2 + (coverage_plot_lim$sd.y / coverage_plot_lim$normalized.y)^2 )  #calculate error propagation
    coverage_plot_lim$sd[coverage_plot_lim$name.x == datasets_list[1]] <- 0 #remove sd data for the reference dataset
    coverage_plot_lim$name <- coverage_plot_lim$name.x      
  }
  
  #logaritmic or linear scale
  if(log_transform == 0) {
    max_x <- max(coverage_plot_lim$normalized, na.rm=TRUE) / zoom
    min_x <- min(coverage_plot_lim$normalized, na.rm=TRUE) / zoom
    }
  if(log_transform == 1) {
    coverage_plot_lim$logv <- log2(coverage_plot_lim$normalized)
    coverage_plot_lim$logv_min <- log2(coverage_plot_lim$normalized - coverage_plot_lim$sd)
    coverage_plot_lim$logv_max <- log2(coverage_plot_lim$normalized + coverage_plot_lim$sd)
    max_x <- max(log2(coverage_plot_lim$normalized), na.rm=TRUE) / zoom
    #min_x <- min(log2(coverage_plot_lim$normalized), na.rm=TRUE) / zoom
    min_x <- min(log2(coverage_plot_lim$normalized[coverage_plot_lim$normalized > 0]), na.rm=TRUE) / zoom
  }
  height_x <- max_x - min_x
  
  ggplot() +
  geom_segment(aes(x = genome_regions[1], y = max_x + height_x * 0.2, xend = genome_regions[2], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bcbddc", arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = genome_regions[4], y = max_x + height_x * 0.2, xend = genome_regions[3], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bcbddc", arrow = arrow(length = unit(0.3, "cm"))) +
  geom_hline(yintercept = (max_x + height_x * 0.1), color = "grey", show.legend=FALSE) +
  geom_rect(data=annotation[annotation$V4 == "r", ], aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.05 * V8)), fill = "#FC8D62", size=1, show.legend=FALSE) +
  geom_rect(data=annotation[annotation$V4 == "t", ], aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.05 * V8)), fill = "#E78AC3", size=1, show.legend=FALSE) +
  geom_rect(data=annotation[annotation$V4 == "p", ], aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.05 * V8)), fill = "#66C2A5", size=1, show.legend=FALSE) +
  #geom_rect(data=annotation[annotation$V4 == "o", ], aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.05 * V8)), fill = "#FFD92F", size=1, show.legend=FALSE) +
  #geom_rect(data=operons, aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1)), color = "#CC79A7", size=1.5, show.legend=FALSE) +
  geom_text(data=annotation[order(annotation$V2-annotation$V1), ], aes(x=V1+(V2-V1)/2, y=(max_x + height_x * 0.1 + height_x * 0.08 * V8), label=V5), size=4, show.legend=FALSE, check_overlap = TRUE) +
      
  {if(log_transform == 1) geom_hline(yintercept = 0, color = "#bdbdbd", show.legend=FALSE)} +    
  {if(log_transform == 0) geom_line(data=coverage_plot_lim, aes(x=position, y=(normalized), group=name, color=name), size=1)} + 
  {if(log_transform == 0) geom_point(data=coverage_plot_lim, aes(x=position, y=(normalized), group=name, color=name), size=0.5)} +
  {if(log_transform == 1) geom_point(data=coverage_plot_lim, aes(x=position, y=(logv), group=name, color=name), size=0.5)} +
  {if(log_transform == 0) geom_ribbon(data=coverage_plot_lim, aes(x=position, ymax=normalized+sd, ymin=normalized-sd, fill=name), alpha=0.2)} +
  {if(log_transform == 1) geom_ribbon(data=coverage_plot_lim, aes(x=position, ymax=logv_max, ymin=logv_min, fill=name), alpha=0.2)} +
  {if(differential_display != 0) geom_text(aes(label = paste("Reference: ", datasets_list[1]), y = min_x, x = seq_end), size = 5, vjust = "inward", hjust = "inward", show.legend=FALSE)} +
    
  coord_cartesian(ylim=c(min_x, max_x + height_x * 0.2), xlim=c(seq_start, seq_end)) +    
  scale_x_continuous(expand = c(0, 0)) +  
  theme(text = element_text(size=15), legend.position=c(0, 0.8), legend.justification = c(0, 1), 
          legend.title = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA), legend.background=element_rect(colour = NA, fill = NA),
          panel.background = element_rect(fill = 'white', colour = 'white'))#, plot.margin = margin(-1, 0, 0, 0, "cm") )
}