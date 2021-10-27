# Analysis of cell distribution within the layers of the human cortex
# input data: QuPath .txt file and data from previous R analysis
#
# AKA
#
# Alessia's second attempt to write code from scratch...

install.packages("tidyverse")
library(tidyverse)
library(reshape)
library(data.table)

##### 0. Calculate layers thickness from reference image #####

## If you can, import layers_proportion.csv instead

## Otherwise..
## Define thickness of each cortical layer, calculate proportion and save as csv file
reference_layers <- function(){
  # extract layer thickness from quantification from ImageJ
  layers_thickness <- read_csv("Cortical layers thickness.csv", 
                               col_names = TRUE)
  
  layers_thickness_lenght <- select(layers_thickness, Length)
  layers_thickness_lenght
  
  # add extra columns for layer number and percentage
  layers_thickness_lenght$Layer <- "whole cortex"
  layers_thickness_lenght$Percentage <- 100.00
  
  # calculate total cortex thickness
  tot_cortex <- layers_thickness_lenght[[1,1]]
  
  # calculate real total cortex thickness (sum of all layers)
  real_tot_cortex <- sum(layers_thickness_lenght$Length)-
    layers_thickness_lenght[[1,1]]
  
  ## DON'T RUN
  # calculate single layers thickenss
  #for (layer in c(2:nrow(layers_thickness_lenght))){
  #  layerx <- layers_thickness_lenght[[layer, 1]]
  #  layer_name <- paste0("layer_", layer-1)
  #  assign(layer_name, layerx)
  #}
  
  # calculate proportion of each layer compared to total and insert in tibble
  for (layer in c(2:nrow(layers_thickness_lenght))){
    layerx <- layers_thickness_lenght[[layer, 1]]/real_tot_cortex*100
    layer_name <- paste0("layer_", layer-1)
    assign(layer_name, layerx)
    layers_thickness_lenght[[layer, 2]] <- layer_name
    layers_thickness_lenght[[layer, 3]] <- layerx
  }
  
  # now the proportion are correct as the sum of each percentage = 100
  sum(layers_thickness_lenght$Percentage)
  # = 200 
  
  # we should also substitute the whole cortex value
  layers_thickness_lenght[[1,1]] <- real_tot_cortex
  
  # just checking all is fine
  if ((sum(layers_thickness_lenght$Length)-real_tot_cortex) == real_tot_cortex) {
    print("all good")
  } else {
    print("not good")
  }
  
  # Remove whole cortex row and length in pixels (first column)
  layers_thickness_lenght <- layers_thickness_lenght[-1,-1]
  
  # Export as .csv file 
  write_csv(layers_thickness_lenght, "layers_proportion.csv")
}
reference_layers()

##### 1. Calculate size/proportion of each layer in the cortex ######

## Calculate for our acquired image the proportion of each layer
layers_size_cum <- function(){
  # Import allcells_xvalue file
  allcells_xvalues <- read.csv(paste0("tables/", marker_1_name, " and ", marker_3_name, " cells x values.csv"))
  
  # Import reference layer thickness
  layers_proportion <- read.csv("layers_proportion.csv")

  # add to the table the total length of the cortex in the experiment
  exp_tot_cortex <- max(allcells_xvalues$Distance)
  
  # add new column in layers_proportion
  layers_proportion$Lenght_layers_exp <- 0
  
  # calculate size of each layer in the experiment
  for (layer in c(1:nrow(layers_proportion))){
    layer_perc_exp <- ((exp_tot_cortex/100)*layers_proportion[[layer,2]])
    layers_proportion[[layer, 3]] <- layer_perc_exp
  }
  
  # just checking all is fine
  sum_layers <- sum(layers_proportion$Lenght_layers_exp)
  sum_layers <- round(sum_layers, 1)
  
  if ((sum_layers) == exp_tot_cortex) {
    print("the sum of each layer length is equal to the total lenght of the cortex")
  } else {
    print("NOT GOOD!")
  }
  
  # create a new column for cumulative distance
  layers_proportion$Lenght_layers_exp_cum <- layers_proportion$Lenght_layers_exp[[1]]
  
  # fill up the column 
  for (layer in 2:nrow(layers_proportion)){
    layers_proportion$Lenght_layers_exp_cum[[layer]] <- 
      (layers_proportion$Lenght_layers_exp_cum[[layer-1]] + 
         layers_proportion$Lenght_layers_exp[[layer]])
  }
  
  # check everything is fine
  if (sum_layers == (round(layers_proportion$Lenght_layers_exp_cum[[6]],1))) {
    print("The cumulative lenght of layers is equal to the total lenght of the cortex")
  } else {
    print("NOT GOOD!")
  }
  
  # Remove not needed columns
  layers_proportion <- layers_proportion[,-3]
  layers_proportion <- layers_proportion[,-2]
  layers_proportion
}
layers_proportion <-layers_size_cum()

##### 2. Calculate number of cells per layer #####
cells_in_layers <- function(){
  # create a table for the future data
  cells_per_layer <- layers_proportion
  # add 6 columns to tables
  for (i in 3:8){
    cells_per_layer[,i] <- NA
  }
  # change column names
  colnames(cells_per_layer) <- c("Layers", 
                                 "Length",
                                 paste0("tot ", marker_1_name),
                                 paste0(marker_1_name, "/", substr(marker_3_name,1,nchar(marker_3_name)-1), "-"), 
                                 marker_2_name, 
                                 paste0("tot ", marker_3_name),
                                 paste0(marker_3_name, "/", substr(marker_1_name,1,nchar(marker_1_name)-1), "-"), 
                                 marker_4_name)
  # add cells totals to each column + subset # Layer_1
  layer = 1
  for (cell_type in c(3:ncol(cells_per_layer))){
    cells_layer <- subset(allcells_xvalues, Distance <= cells_per_layer[[layer,2]])
    cells_per_layer[layer, 3] <-  sum(cells_layer$Cell_type == marker_1_name | cells_layer$Cell_type == marker_2_name)
    cells_per_layer[layer, 4] <-  sum(cells_layer$Cell_type == marker_1_name)
    cells_per_layer[layer, 5] <-  sum(cells_layer$Cell_type == marker_2_name)
    cells_per_layer[layer, 6] <-  sum(cells_layer$Cell_type == marker_3_name | cells_layer$Cell_type == marker_4_name)
    cells_per_layer[layer, 7] <-  sum(cells_layer$Cell_type == marker_3_name)
    cells_per_layer[layer, 8] <-  sum(cells_layer$Cell_type == marker_4_name)
    layer_name <- paste0("cell_xvalue_layer", layer)
    assign(layer_name, cells_layer)
  }
  # add cells totals to each column + subset # all other layers
  for (layer in c(2:nrow(cells_per_layer))) {
    for (cell_type in c(3:ncol(cells_per_layer))){
      cells_layer <- subset(allcells_xvalues, Distance <= cells_per_layer[[layer,2]] & Distance >= cells_per_layer[[layer-1,2]])
      cells_per_layer[layer, 3] <-  sum(cells_layer$Cell_type == marker_1_name | cells_layer$Cell_type == marker_2_name)
      cells_per_layer[layer, 4] <-  sum(cells_layer$Cell_type == marker_1_name)
      cells_per_layer[layer, 5] <-  sum(cells_layer$Cell_type == marker_2_name)
      cells_per_layer[layer, 6] <-  sum(cells_layer$Cell_type == marker_3_name | cells_layer$Cell_type == marker_4_name)
      cells_per_layer[layer, 7] <-  sum(cells_layer$Cell_type == marker_3_name)
      cells_per_layer[layer, 8] <-  sum(cells_layer$Cell_type == marker_4_name)
      layer_name <- paste0("cell_xvalue_layer", layer)
      assign(layer_name, cells_layer)
    }
  }
  # check that all is fine (number of each cells type per layer equal to its total)
  
  if (results_1[[1,2]] == sum(cells_per_layer[,3]) &
      results_1[[2,2]] == sum(cells_per_layer[,5]) &
      results_2[[1,2]] == sum(cells_per_layer[,6]) &
      results_2[[2,2]] == sum(cells_per_layer[,8])) {
    print("The number of each cells type per layer equal to its total")
  } else {
    print("NOT GOOD!!")
  }
  # check that all is fine (total number of cells per layer, equals number of cells in subset file)
  if (nrow(cell_xvalue_layer1) == sum(cells_per_layer[1,c(3,6)]) &
      nrow(cell_xvalue_layer2) == sum(cells_per_layer[2,c(3,6)]) &
      nrow(cell_xvalue_layer3) == sum(cells_per_layer[3,c(3,6)]) &
      nrow(cell_xvalue_layer4) == sum(cells_per_layer[4,c(3,6)]) &
      nrow(cell_xvalue_layer5) == sum(cells_per_layer[5,c(3,6)]) &
      nrow(cell_xvalue_layer6) == sum(cells_per_layer[6,c(3,6)])) {
    print("The total number of cells per layer, equals number of cells in subset file")
  } else {
    print("NOT GOOD!!")
  }
  cells_per_layer
  write_csv(cells_per_layer, paste0("tables/", marker_1_name, " and ", marker_3_name, " cells per layers.csv"))
}
all_markers_layers <- cells_in_layers()
# all_markers_layers <- cells_per_layer

##### 3. Plot distribution of cell types per layer #####
make_plots_layers <- function(){
  # create a new tibble with data from layers only (remove whole cortex and thickness of layers)
  all_markers_layers_clean <- all_markers_layers[,-2]
  # transform tibble "wide format" to "long format" for plotting
  long <- melt(setDT(all_markers_layers_clean), variable.name = "Cell_type")
  # plot
  plot <- ggplot(data = long, 
                 mapping = aes(y = reorder(Layers, desc(Layers)), x = value, fill = Cell_type)) +
    geom_col() +
    labs(title = "Distibution of cells along cortical layers", 
         x = "no. of cells",
         y = "") +
    facet_wrap(~Cell_type)
  print(plot)
  # save as PDF
  pdf(paste0("plots/", marker_1_name, " and ", marker_3_name, " cells in layers (separated).pdf"))
  print(plot)
  dev.off()
}

make_plots_layers()
