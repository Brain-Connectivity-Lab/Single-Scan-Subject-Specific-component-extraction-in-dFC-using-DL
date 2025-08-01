function CM_run_COBE_Data = CBL_Obtain_COBE_Normal_Corr_Data( Time_Course_run,Component_remv ) %Corr_Common_cell_run , Corr_Orig_run

%Getting the components to be removed and the applying COBE
CM_run_COBE_Data = my_cobec(Time_Course_run, Component_remv);

end

