clear;clc;
pps_path_mode =  'moderate' ;% 'conservative' 'ambitious'
main_run_SS_pps(pps_path_mode);
main_run_trans_pps(pps_path_mode);
% clear;
% pps_path_mode = 'conservative' ;
% main_run_SS_pps(pps_path_mode);
% main_run_trans_pps(pps_path_mode);
clear;
pps_path_mode = 'ambitious' ;
main_run_SS_pps(pps_path_mode);
main_run_trans_pps(pps_path_mode);


% matlab -batch "pps_path_mode = 'moderate' ;main_run_SS_pps(pps_path_mode);main_run_trans_pps(pps_path_mode);"
% matlab -batch "pps_path_mode ='conservative' ;main_run_SS_pps(pps_path_mode);main_run_trans_pps(pps_path_mode);"
% matlab -batch "pps_path_mode = 'ambitious' ;main_run_SS_pps(pps_path_mode);main_run_trans_pps(pps_path_mode);"