# Plavisto - Plastid Genome Visualization Tool
# Version 1.0 (01.19.2022)
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

number_of_plots <- 2
plot_height <- 200
genome_length <<- 154477 #https://doi.org/10.1093/dnares/6.5.283
genome_regions <<- c(84170, 84170+26264, 84170+26264+17780, 154477) #Genomic regions IR1start, IR1end, IR2start, IR2end
options(dplyr.summarise.inform = FALSE)


library("dplyr")
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("jsTreeR")) {
  install.packages("jsTreeR", dependencies = TRUE)
  library(jsTreeR)
}
if (!require("tidyr")) {
  install.packages("tidyr", dependencies = TRUE)
  library(tidyr)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("memoise")) {
  install.packages("memoise", dependencies = TRUE)
  library(memoise)
}