# wheres_my_neuron

Welcome to Alessia's first contribution to the bioinformatic field!

I am very new to coding, please be gentle. Which also means, all supporting description here are coding dumbies-proof! 

This simple code was written to specifically process the .txt output file of the image analysis software QuPath, for quantifying the total number of cells detected for a defined marker, the proportion of them co-expressing a second defined marker, as well as their distribution within the tissue (using X coordinate). 


1. Download QuPath

QuPath is a fantastic open source software for image analysis that you can download here: https://qupath.github.io/

2. Analyse your image in QuPath

Your image can be anything, but wheres_my_neuron has been currently only tested on immunofluorescence staining using two antibodies (markers). 

Once opened a project and uploaded your image, you can detect cells using the "Cell detection" function: 
https://qupath.readthedocs.io/en/stable/docs/tutorials/cell_detection.html
The marker you use to detect cells in the first place will be your marker_1

Then, you can classify your detected cells for their co-expression of a second marker of interest (marker_2), using the "Cell classification" tool: https://qupath.readthedocs.io/en/stable/docs/tutorials/cell_classification.html

3. Export your analysis

At the end of this analysis you will have all your detected objects classified for being only positive for marker_1, or double positive for marker_1 and marker_2. 
Make sure you give marker_1 and marker_2 the correct names as they will be used by wheres_my_neuron to make plots and tables. 

wheres_my_neuron can support up to teo combination of markers (marker_1 positive only, marker_1 and marker_2 double positive, marker_3 positive only, marker_3 and marker_4 double positive)

You can export the .txt file with information from all your detected objects via the measurement table.

4. Prepare folders for data processing

Now copy your .txt file to a "txt" folder in the same folder as the wheres_my_neuron.R file. You're ready to go!
