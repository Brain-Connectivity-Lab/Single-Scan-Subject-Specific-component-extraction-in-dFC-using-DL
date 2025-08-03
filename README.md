# Single-Scan-Subject-Specific-component-extraction-in-dFC-using-DL
Code for the paper Single Scan, Subject-Specific component extraction in dynamic functional connectivity using Dictionary Learning. 

## How to run?
First, download all the codes and files in this repo.

## To just view the figures provided in the manuscript
run the results_and_figures.m file.
  i) Before running change the working directory path. Put the path of the dir where the repo was downloaded.
  ii) Just by running thr results_and_figures.m file now you should be able to genrate the figures corresponding to the Schefer_HCP_419 roi atlas.

## To recreate everything, you will have to do the following
1) Download the HCP data and run divide_in_3_and_get_reg_ts_files.m
     i) Downlooad HCP data with ICA fix from connectome db (1200 Subject Data release).
     ii) Open the divide_in_3_and_get_reg_ts_files.m file and 
     iii) put the path of the folder corresponding to the downloaded files in data_folder
     iv) Give a path to an empty working directory folder.
     v) This will create the average timeseries for every subject for the atlas choosen, randomly divide the subjects in 3 groups and save the .mat files in the working directory
