# Analysis of cell number and distribution along the layers of the human cortex
# input data: QuPath .txt file
#
# AKA
#
# Alessia's first attempt to write code from scratch...

install.packages("tidyverse")
library(tidyverse)
library(gridExtra)

#### 1. Convert .txt files to data frames and .csv files ####
txttodfs <- function(){
  filelist = list.files(path = "txt",  pattern = ".txt") # Identify files to be converted
  result <- data.frame()
  for (i in 1:length(filelist)){ #Loop through the numbers of files 
    input <- filelist[i] # Extract file to work on
    output <- paste0(gsub("\\.txt$", "", input), ".csv") # Change end of file name
    print(paste("Processing the file:", input)) # Show file being processed
    data = read.delim(paste0(file.path(("txt"), input)), header = TRUE)   # make data frame from .txt file
    dfname <- paste0(gsub("\\.txt$", "", input))
    assign(paste('df', sep="_", dfname),data) # create a data frame to hold results
    result <- bind_rows(result, data)
  }
  result
}

all_markers = txttodfs()

#### 2. Separate data frame by marker and save their names ####

# get markers names
markers <- unique(all_markers$Class)
markers <- markers[markers != ""] # remove empty values

# subset for markers
for (i in 1:length(markers)) {
  marker_id <- markers[i]
  print(paste0("Marker #", i , " is: ", marker_id))
  subset_marker <- subset(all_markers, subset = all_markers$Name == markers[i])
  assign(paste0("marker_", i), subset_marker)
  assign(paste0("marker_", i, "_name"), marker_id)
}
rm(subset_marker)
rm(i)
rm(marker_id)
print("Now you need to merge markers")

if(!exists("marker_1_name")) { marker_1_name <- "missing" }
if(!exists("marker_2_name")) { marker_2_name <- "missing" }
if(!exists("marker_3_name")) { marker_3_name <- "missing" }
if(!exists("marker_4_name")) { marker_4_name <- "missing" }

if(!exists("marker_1")) { marker_1 <- data.frame()}
if(!exists("marker_2")) { marker_2 <- data.frame()}
if(!exists("marker_3")) { marker_3 <- data.frame()}
if(!exists("marker_4")) { marker_4 <- data.frame()}

#### 3. Merge markers with same marker ####
merge_markers <- function(){
  msg <- "Do you want to merge markers 1+2 and 3+4?"
  reply <- askYesNo(msg, default = TRUE, 
                    prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
  if (reply){
    print("Ok I will merge markers 1+2 and 3+4")
    ## Merge markers 1+2
    # Names of the first two markers to merge (1+2)
    markers_1_2_names <<- markers[str_detect(markers, paste0("^", markers[1]))]
    
    # merge data frames in a new data frame
    markers_1_2 <<- bind_rows(marker_1, marker_2)
    
    ## Merge markers 3+4
    # Names of the last two markers to merge (3+4)
    markers_3_4_names <<- markers[str_detect(markers, paste0("^", markers[3]))]
    
    # merge data frames in a new data frame
    markers_3_4 <<- bind_rows(marker_3, marker_4)
  }
}
merge_markers()

if(!exists("markers_1_2_names")) { markers_1_2_names <- "missing" }
if(!exists("markers_3_4_names")) { markers_3_4_names <- "missing" }

#### 4. Calculate number of cells, percentages and write a result table #### 
 
# Create a csv folder for the result file
dir.create("tables")

make_tables <- function(){
  # Make a table for results for markers 1 an 2
  numbers_markers_1_2 <- function(){
    { if (marker_1_name != "missing") {
      cols_title <- c(paste0("Total ", markers[1]), paste0("Total ", markers[2]), paste0("% ", markers[2]))
      cols_values <- c(nrow(markers_1_2), 
                       nrow(marker_2), 
                       (nrow(marker_2)/nrow(markers_1_2))*100)
      results_1_2 <- data.frame(cols_title, cols_values)
      colnames(results_1_2) <- c(" ", " ")
      pdf(paste0("tables/", marker_1_name, " cells summary.pdf"))
      grid.table(results_1_2)
      dev.off()
    } else {
      print("Marker 1 and 2 missing")
      }
    }
    results_1_2
  }
  results_1 <<- numbers_markers_1_2()
  
  # Make a table for results for markers 3 and 4
  numbers_markers_3_4 <- function(){
    { if (marker_3_name != "missing") {
      cols_title <- c(paste0("Total ", markers[3]), paste0("Total ", markers[4]), paste0("% ", markers[4]))
      cols_values <- c(nrow(markers_3_4), 
                       nrow(marker_4), 
                       (nrow(marker_4)/nrow(markers_3_4))*100)
      results_3_4 <- data.frame(cols_title, cols_values)
      colnames(results_3_4) <- c(" ", " ")
      pdf(paste0("tables/", marker_3_name, " cells summary.pdf"))
      grid.table(results_3_4)
      dev.off()
    } else {
      print("Marker 3 and 4 missing")
    }
    }
    results_3_4
  }
  results_2 <<- numbers_markers_3_4()
  }

make_tables()

#### 5. Check calculations are correct ####
check_numbers <- function() {
  if (marker_1_name != "missing") {
    if (nrow(marker_1)+nrow(marker_2) == nrow(markers_1_2) &
       nrow(markers_1_2) == results_1[1,2] &
       nrow(marker_2) == results_1[2,2]
  ){
    print("all good")
  } else {
    print("NOT GOOD!")
  }
  } else {
    print("Marker 1 missing")
  }
  if (marker_3_name != "missing") {
    if ( nrow(marker_3)+nrow(marker_4) == nrow(markers_3_4) &
         nrow(markers_3_4) == results_2[1,2] &
         nrow(marker_4) == results_2[2,2]
    ){
      print("all good")
    } else {
      print("NOT GOOD!")
    }
  } else {
    print("Marker 3 missing")
  }
}
  
check_numbers()

#### 6. Calculate and plot distribution ####

make_plots <- function(){
  # create a directory for saving pdf files
  dir.create("plots")
  # Extract column of X values of distribution along the cortical layers
  allcells_xvalues <<- tibble("Cell_type" =  all_markers$Name,
                                 "Distance" =  all_markers$Centroid.X.µm)
  markers_1_2_xvalues <<- tibble("Cell_type" =  markers_1_2$Name,
                                "Distance" =  markers_1_2$Centroid.X.µm)
  markers_3_4_xvalues <<- tibble("Cell_type" =  markers_3_4$Name,
                                 "Distance" =  markers_3_4$Centroid.X.µm)
  
  # Create the histogram plot for marker_1_2 and save as PDF
 if (marker_1_name != "missing") {
   plot <- ggplot(data = markers_1_2_xvalues, 
                  mapping = aes(x = Distance,
                                fill = Cell_type)) +
     labs(title = "Distibution of cells along the cortex", 
          x = "Distance from top of the cortex (um)", 
          y = "no. of cells") +
     geom_histogram(binwidth = 200, 
                    alpha = 0.6) +
     coord_flip() +
     scale_x_reverse()
   
   plot_separated <- plot + facet_wrap(~Cell_type)
   
   pdf(paste0("plots/all_", marker_1_name, "_cells.pdf"))
   print(plot)
   dev.off()
   
   pdf(paste0("plots/all_", marker_1_name, "_cells (separated).pdf"))
   print(plot_separated)
   dev.off()
 } 

  # Create the histogram plot for marker_3_4 and save as PDF
  if (marker_3_name != "missing") {
    plot <- ggplot(data = markers_3_4_xvalues, 
                   mapping = aes(x = Distance,
                                 fill = Cell_type)) +
      labs(title = "Distibution of cells along the cortex", 
           x = "Distance from top of the cortex (um)", 
           y = "no. of cells") +
      geom_histogram(binwidth = 200, 
                     alpha = 0.6) +
      coord_flip() +
      scale_x_reverse()
    
    plot_separated <- plot + facet_wrap(~Cell_type)
    
    pdf(paste0("plots/all_", marker_3_name, "_cells.pdf"))
    print(plot)
    dev.off()
    
    pdf(paste0("plots/all_", marker_3_name, "_cells (separated).pdf"))
    print(plot_separated)
    dev.off()
  } else {
    print("Marker 3 and 4 missing")
    }
}
make_plots()

#### 7. Check cortical thicken and recalculate, if necessary ######

# What's the total length of the cortex?
check_lenght <- function(){
  print(paste0("Total length of the cortex is ", max(allcells_xvalues$Distance), "um"))
  if ((max(allcells_xvalues$Distance)) > 3300) {
    msg <- (paste0("Cortex is ", max(allcells_xvalues$Distance), "um thick. Do you want to cut the max of the cortex to 3300um?"))
    reply <<- askYesNo(msg,  prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
  } else {
    print("Cortex leght is less than 3300")
  }
# Sub setting cells beyond 3300um
  if (reply) {
    print("Ok, I will cut to max cortical distance of 3300um")
    all_markers <<- subset(all_markers, all_markers$Centroid.X.µm < 3300)
    marker_1 <<- subset(marker_1, marker_1$Centroid.X.µm < 3300)
    marker_2 <<- subset(marker_2, marker_2$Centroid.X.µm < 3300)
    marker_3 <<- subset(marker_3, marker_3$Centroid.X.µm < 3300)
    marker_4 <<- subset(marker_4, marker_4$Centroid.X.µm < 3300)
    # Merge markers
    merge_markers()
    # Make a table 
    make_tables()
    } else {
      print("ok, I won't change the data")
      }
  # Check everything is fine
  check_numbers()
  # Make plots
  make_plots()
  print("Data has been recalculated for maximum cortex thickness of 3300um and pdf re-printed")
} 
check_lenght()

#### 7. Save #####
write_csv(allcells_xvalues, paste0("tables/", marker_1_name, " and ", marker_3_name, " cells x values.csv"))







