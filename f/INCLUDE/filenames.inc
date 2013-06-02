********************************************************************
*                                                                  *
*                  Include file                                    *
*  the filenames for randomnumbers, configuration ...              *
*  unit numbers for the files                                      *
********************************************************************

       character*(49) random_number_file, u_conf_file, a_conf_file,
     &                action_file, joblog_file, lastconfig_file, 
     &                corr_file, afield_file, checkp_file, 
     &                polyak_file, corra2_file

       character*(9)  inputfile

       integer        random_number_unit, u_unit, a_unit,
     &                action_unit, joblog_unit,input_unit,
     &                corr_unit, afield_unit, checkp_unit, 
     &                polyak_unit, corra2_unit

***********************************************************************
*
*                     THE COMMON BLOCK
*
***********************************************************************

       common /filenames/ random_number_file, u_conf_file, a_conf_file,
     &                    action_file, joblog_file,inputfile, 
     &                    corr_file, afield_file, checkp_file,
     &                    polyak_file, corra2_file
       common /units/ random_number_unit, u_unit, a_unit,
     &                action_unit, joblog_unit, input_unit,
     &                corr_unit, afield_unit, checkp_unit,
     &                polyak_unit, corra2_unit   

***********************************************************************
***********************************************************************







