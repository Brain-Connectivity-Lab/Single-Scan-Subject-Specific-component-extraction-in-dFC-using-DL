# Single-Scan-Subject-Specific-component-extraction-in-dFC-using-DL
Code for the paper Single Scan, Subject-Specific component extraction in dynamic functional connectivity using Dictionary Learning. 

These scripts were written by Pratik Jain. If you have any questions, email jainpratik412[at]gmail[dot]com

## How to run?
First, download all the codes and files in this repo.

## To just view the figures provided in the manuscript
run the results_and_figures.m file.

  i) Before running, change the working directory path. Put the path of the directory where the repo was downloaded.
  
  ii) Just by running the results_and_figures.m file, now you should be able to generate the figures corresponding to the Schefer_HCP_419 roi atlas.
  

## To recreate everything, you will have to do the following

### 1) Download the HCP data and run divide_in_3_and_get_reg_ts_files.m
   
     i) Download HCP data with ICA fix from connectome db (1200 Subject Data release).

     ii) Open the divide_in_3_and_get_reg_ts_files.m file and
   
     iii) put the path of the folder corresponding to the downloaded files in data_folder
   
     iv) Give a path to an empty working directory folder.
   
     v) This will create the average timeseries for every subject for the atlas chosen, randomly divide the subjects into 3 groups, and save the .mat files in the working directory
   
### 2) Run the parameter_search.m file
   
     i) Change the working dir path and run.
  
     ii) This script uses the average timeseries created from the previous script divide_in_3_and_get_reg_ts_files.m, and performs a grid search to find the best parameters for the specified sub_block, atlas.

### 3) Run the Create_avg_para_search_across_blocks.m
   
     i) Change the working dir path and run.
   
     ii)This Script uses the results of the parameter check from the previous script parameter_search.m, and averages the optimal parameters across the blocks

### 4) Run the main_results_5x2cv.m
   
     i) Change the working dir path and run.
   
     ii)This Script uses the optimal parameters from the previous script Create_avg_para_search_values_across_blocks.m and performs 5x2 CV

### 5) Run the restricted_data_analysis.m

    If you have the restricted data CSV, only then you will be able to run this.
     i) Change the working dir path and run.
   
     ii)This Script uses the optimal parameters from the script Create_avg_para_search_values_across_blocks.m and performs 5x2 CV on the monozygotic, dizygotic, and unrelated subjects

### 6) Run the from hcp_dict_nki_5x2cv.m

    i) Change the working dir path.
    ii) Download the NKI Rockland sample data and extract the timeseries corresponding to the atlas.
    iii) I have uploaded the NKI mat files for the s400_hcp atlas that you can use to see how the shape of the averate timeseries should be for the code to work. 

### 7) Run the results_and_figures.m file

    Once all the above codes are executed, run the results_and_figures.m file. It will generate the figures for the dataset that you have..
     
