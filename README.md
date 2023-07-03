# Tensor-Decomposition-for-Analysis-of-fNIRS-Signals
Abstract: The analysis of functional near-infrared spectroscopy (fNIRS) signals has not kept pace with the increased use of fNIRS in the behavioral and brain sciences. The popular grand averaging method collapses the oxygenated hemoglobin data within a predefined time of interest window and across multiple channels within a region of interest, potentially leading to a loss of important temporal and spatial information. On the other hand, the tensor decomposition method can reveal patterns in the data without making prior assumptions of the hemodynamic response and without losing temporal and spatial information. The aim of the current study was to examine whether the tensor decomposition method could identify significant effects and novel patterns. We used  tensor decomposition (i.e., canonical polyadic and Tucker decompositions) to analyze the significant differences in the hemodynamic response patterns across conditions. We will add a link to our paper that has more details about the synthesized data and the methods.

Prerequisites: To run this code, you need Matlab and Tensorlab 3.0 [Computer software] (2016), retrieved from https://www.tensorlab.net/. Tensorlab 3.0 must be in the tensor decomposition folder. You will also need fNIRS data as input. The data should be in a txt format, with one file for each subject. Each row of the file is a time point and each column is a channel. All fNIRS data needs to be saved in a folder named "Text Files" with one subfolder for each condition. For a different dataset, you should change the directory and other parameters in "Data directory and hyperparameters" section in the main code for each method. After running CPD or TD as follows, the results are saved in the "Results" folder.

Notes for running the code: There is a folder to run a different type of tensor decomposition (e.g., CPD and Tucker decomposition [TD]). The headers (e.g., 1 and 2) are the main steps of running the tensor decompositions and the subheaders (e.g., 1.1, 1.2, etc.) are supporting functions required from the according headers. The folder "Text Files" contains a sample of the data used in our work for demonstration only. The whole data set is available on reasonable request. 

For CPD, run  Tensor_Construction_Tensor_Decomposition.m then Excluding_Irrelevant_Components_Determination_of_TOI_and_ROI.m:

	1- Tensor_Construction_Tensor_Decomposition.m: forms the tensor and decomposes the tensor into components
	
 		1.1 - Text Files: a folder with fNIRS data. Within Text Files is one subfolder for each condition. The files in the folder have to be in .txt format. 
		1.2 - load_data.m: Reads channel indices from the .txt file
   		1.3 - read_condition_tensor.m: Reads all the .txt files for one condition that is in one folder.
   		1.4 - tensorlab_2016-03-28: package downloaded from Tensorlab 3.0 [Computer software]. (2016). Retrieved from https://www.tensorlab.net/
		1.5 - CPD_plots: Generate a figure for each CPD component contains subplots of temporal (time), spectral (frequency), spatial (channels), and subjects. Plotting channels activations on the left and right hemisphere. Taking condition indeces to use in bar plots.
		1.6 - left pic.gif: .gif file of the left hemisphere
		1.7 - right pic.gif: .gif file of the right hemisphere
		1.8 - CPD_plots_mean.m:  Generate a figure for each CPD component contains subplots of temporal (time), spectral (frequency), spatial (channels), and subjects. Plotting channels activations on the left and right hemisphere. Taking condition indeces to use in bar plots. Plotting the mean of the signatures of each group of subjects
	2 - Excluding_Irrelevant_Components_Determination_of_TOI_and_ROI.m: Applies 2-way ANOVA on the decomposed subjects' signatures after applying tensor-decomposition. Excludes irrelevant components. Takes the all the statistically siginifcant components that are highly corrlated and weights them in proportion of the occurance. Generate a figure for each CPD component contains subplots of temporal (time) spectral (frequency), spatial (channels), and subjects. Plotting channels activations on the left and right hemisphere. Plotting the mean of each condition. 


For TD, run  Tensor_Construction_Tensor_Decomposition.m then Excluding_Irrelevant_Components_Determination_of_TOI_and_ROI.m:

	1 - Tensor_Construction_Tensor_Decomposition.m: forms the tensor and decomposes the tensor into components
		1.1 - Text Files: a folder with fNIRS data. Within Text Files is one subfolder for each condition. The files in the folder have to be in .txt format. 
		1.2 - load_data.m: Reads channel indices from the .txt file
   		1.3 - read_condition_tensor.m: Reads all the .txt files for one condition that is in one folder.
   		1.4 - tensorlab_2016-03-28: package downloaded from Tensorlab 3.0 [Computer software]. (2016). Retrieved from https://www.tensorlab.net/
		1.5 - TD_plots: Generate a figure for each component [temporal (time) spatial (channels), and subjects]. Plotting channels activations on the left and right hemisphere. Taking condition indeces to use in bar plots.
		1.6 - left pic.gif: .gif file of the left hemisphere
		1.7 - right pic.gif: .gif file of the right hemisphere
		1.8 - TD_plots_mean.m:  Generate figures for component [temporal (time), spatial (channels), and subjects.] Combines multiple component into one plot for easy viewing. Plotting channels activations on the left and right hemisphere. Taking condition indeces to use in bar plots. Plotting the mean of the signatures of each group of subjects
	2 - Excluding_Irrelevant_Components_Determination_of_TOI_and_ROI.m: perform ANOVA 
		2.1 - TD_anova_plot.m: makes individual plots for each mode 
