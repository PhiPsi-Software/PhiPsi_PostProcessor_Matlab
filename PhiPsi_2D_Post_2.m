%     .................................................
%             ____  _       _   ____  _____   _        
%            |  _ \| |     |_| |  _ \|  ___| |_|       
%            | |_) | |___   _  | |_) | |___   _        
%            |  _ /|  _  | | | |  _ /|___  | | |       
%            | |   | | | | | | | |    ___| | | |       
%            |_|   |_| |_| |_| |_|   |_____| |_|       
%     .................................................
%     PhiPsi:     a general-purpose computational      
%                 mechanics program written in Fortran.
%     Website:    http://phipsi.top                    
%     Author:     Fang Shi  
%     Contact me: shifang@ustc.edu.cn   

%-------------------------------------------------------------------
%--------------------- PhiPsi_Post_Plot2 ---------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
clear all; close all; clc; format compact;  format long;
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack Defor_Factor
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves num_Na_Crack
global Plot_Aperture_Curves Plot_Pressure_Curves Plot_Velocity_Curves Num_Step_to_Plot
global Key_TipEnrich Plot_Quantity_Curves Plot_Concentr_Curves Itera_Num
global Na_Crack_X Na_Crack_Y num_Na_Crack Itera_HF_Time Key_HF_Analysis
global num_Hole Hole_Coor   
global num_Circ_Inclusion num_Poly_Inclusion Key_Time_String
global Key_Animation Key_Ani_Ave
global Time_Delay


%-------------------------- Settings -------------------------------
% Add path of source files.
addpath('src_fcw')
addpath('src_geom3d')
addpath('src_meshes3d')
addpath('src_phipsi_post_animate')
addpath('src_phipsi_post_cal')
addpath('src_phipsi_post_main')
addpath('src_phipsi_post_plot')
addpath('src_phipsi_post_read')
addpath('src_phipsi_post_tool')

% Set default figure colour to white.
set(0,'defaultfigurecolor','w')

% Set default figure visible off.
set(0,'DefaultFigureVisible','off')

% Output information of matlab command window to log file.
diary('Command Window.log');        
diary on;
% Display welcome information.
Welcome
       
tic;
Tclock=clock;
Tclock(1);

disp([' >> Start time is ',num2str(Tclock(2)),'/',num2str(Tclock(3)),'/',num2str(Tclock(1))...
     ,' ',num2str(Tclock(4)),':',num2str(Tclock(5)),':',num2str(round(Tclock(6))),'.'])
disp(' ') 

% Make the "patch" method supported by "getframe", added in version 4.8.10
% See more: http://www.mathworks.com/support/bugreports/384622
opengl('software')      

%----------------------- Pre-Processing ----------------------------
disp(' >> Reading input file....') 


% -------------------------------------
%   Set color and font
% -------------------------------------                           
PhiPsi_Color_and_Font_Settings

% -------------------------------------
%   Start Post-processor.      
% -------------------------------------   
Key_PLOT   = zeros(6,15);                                   % Initialize the Key_PLOT

%###########################################################################################################
%##########################            User defined part        ############################################
%###########################################################################################################
% Filename='example';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\example';

% Filename='XFEM';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\XFEM';
% Filename='XFEM';Work_Dirctory='X:\PhiPsi_fortran_work_LEFM\PhiPsi_work\XFEM';


% Filename='Brief_intro_1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Brief_intro_1';
% Filename='test_1x1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_1x1';
% Filename='test_1x2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_1x2'; %!!!!!!!!!
% Filename='test_2x2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_2x2';
% Filename='test_3x3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_3x3';
% Filename='test_3x3';Work_Dirctory='X:\PhiPsi fortran work1\PhiPsi_work\test_3x3';
% Filename='test_3x2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_3x2';
% Filename='test_3x2';Work_Dirctory='X:\PhiPsi fortran work1\PhiPsi_work\test_3x2';
% Filename='test_5x5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_5x5';
% Filename='test_5x5_compression';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_5x5_compression';
% Filename='test_7x7';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_7x7';
% Filename='test_11x11';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_11x11';
% Filename='test_13x13';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_13x13';
% Filename='test_23x23';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_23x23';
% Filename='test_23x23';Work_Dirctory='X:\PhiPsi fortran work1\PhiPsi_work\test_23x23';
% Filename='exa_cracks_meet_holes';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\exa_cracks_meet_holes';
% Filename='Crack_Meet_Hole';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole';
% Filename='Crack_Meet_Hole_f';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_f';
% Filename='Crack_Meet_Hole_ff';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_ff';
% Filename='Crack_Meet_Hole_m';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_m';
% Filename='Crack_Meet_Hole_New';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_New';
% Filename='Crack_Meet_Hole_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_fine';
% Filename='Arc_Crack_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Arc_Crack_test';
% Filename='Arc_Crack_test2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Arc_Crack_test2';
% Filename='Crack_Meet_Inclusion_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Inclusion_test';
% Filename='Crack_Meet_Inclusion_31x31';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Inclusion_31x31';
% Filename='Crack_Meet_Inclusion_51x51';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Inclusion_51x51';
% Filename='Crack_Meet_Inclusion_101x101';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Inclusion_101x101';
% Filename='test_55x55';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\test_55x55';
% Filename='Cohesive_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Cohesive_test';
% Filename='Cohesive_test_Tang';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Cohesive_test_Tang';
% Filename='Three_point_bending_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Three_point_bending_test';
% Filename='Multi_Cracks';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Multi_Cracks';
% Filename='exa_random_multi_cracks';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\exa_random_multi_cracks';
% Filename='exa_random_multi_cracks_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\exa_random_multi_cracks_fine';
% Filename='Multi_Cracks_Ya';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Multi_Cracks_Ya';
% Filename='Cross_shaped_Crack';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Cross_shaped_Crack';
% Filename='exa_crossing_cracks';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\exa_crossing_cracks';
% Filename='Snow_shaped_Crack';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Snow_shaped_Crack';
% Filename='uniaxial_tension_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\uniaxial_tension_test';
% Filename='hole_regular_mesh';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\hole_regular_mesh';
% Filename='hole_regular_mesh_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\hole_regular_mesh_tension';
% Filename='K_benchmark_Bigger';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\K_benchmark_Bigger';
% Filename='static_random';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\static_random';
% Filename='Penalty_Test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test';
% Filename='Penalty_Test_11x11';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_11x11';
% Filename='Penalty_Test_11x11_2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_11x11_2';
% Filename='Penalty_Test_11x11_3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_11x11_3';
% Filename='Penalty_Test_11x11_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_11x11_tension';
% Filename='Penalty_Test_21x21';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_21x21';
% Filename='Penalty_Test_77x77_bi-tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Penalty_Test_77x77_bi-tension';
% Filename='Simple';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Simple';
% Filename='Shear';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Shear';
% Filename='brazil_coarse_mesh';Work_Dirctory='D:\PhiPsi_Project\PhiPsi_work\brazil_coarse_mesh';
% Filename='brazil_half';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\brazil_half';
% Filename='brazil_half_with_crack';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\brazil_half_with_crack';
% Filename='brazil_fine_mesh';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\brazil_fine_mesh';
% Filename='recta_uniaxial-tension_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\recta_uniaxial-tension_test';
% Filename='HF_pressure_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_pressure_test';
% Filename='HF_pressure_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_pressure_test';
% Filename='HF_pressure_test2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_pressure_test2';
% Filename='Junction';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Junction';
% Filename='HF_Sys_5x5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_5x5';
% Filename='HF_Sys_5x5_1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_5x5_1';
% Filename='HF_Sys_11x11';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_11x11';
% Filename='E_1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\E_1';
% Filename='E_1_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\E_1_Fixed';
% Filename='HF_Sys_11x11';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_11x11';
% Filename='HF_31x31_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_31x31_fixed';
% Filename='HF_Sys_17x17';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17';
% Filename='HF_Sys_17x17_0.25';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_0.25';
% Filename='HF_Sys_17x17_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_Fixed';
% Filename='HF_Sys_17x17_InSitu_0_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_InSitu_0_5';
% Filename='HF_Sys_17x17_InSitu_4_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_InSitu_4_5';
% Filename='HF_Sys_17x17_InSitu_5_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_InSitu_5_5';
% Filename='HF_Sys_17x17_Static';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_17x17_Static';
% Filename='Paper1_example_1_2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper1_example_1_2';
% Filename='Example_2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Example_2';
% Filename='HF_Sys_21x21';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_21x21';
% Filename='HF_Sys_23x23';Work_Dirctory='X:\PhiPsi fortran work1\PhiPsi_work\HF_Sys_23x23';
% Filename='HF_Sys_25x25';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25';
% Filename='HF_Sys_25x25_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_fixed';
% Filename='HF_Sys_25x25_Right_Force';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_Right_Force';
% Filename='HF_Sys_25x25_Right_Neg_Force';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_Right_Neg_Force';
% Filename='HF_Sys_25x25_Top_Neg_Force';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_Top_Neg_Force';
% Filename='HF_Sys_25x25_Geostress';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_Geostress';
% Filename='HF_Sys_25x25_InSitu_4_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_InSitu_4_5';
% Filename='HF_Sys_25x25_InSitu_0_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_InSitu_0_5';
% Filename='HF_Sys_45x45';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_45x45';
% Filename='HF_Sys_45x45_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_45x45_fixed';
% Filename='HF_Sys_45x45_InSitu_4_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_45x45_InSitu_4_5';
% Filename='HF_Sys_55x55';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55';
% Filename='HF_Sys_55x55_unixal_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_unixal_tension';
% Filename='HF_Sys_55x55_tension_and_compres';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_tension_and_compres';
% Filename='HF_Sys_55x55_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_fixed';
% Filename='HF_Sys_55x55_InSitu_4_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_InSitu_4_5';
% Filename='HF_Sys_55x55_InSitu_0_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_InSitu_0_5';
% Filename='HF_Sys_55x55_InSitu_7_10';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_InSitu_7_10';
% Filename='HF_Sys_55x55_InSitu_9_10';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_55x55_InSitu_9_10';
% Filename='HF_Sys_65x65';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_65x65';
% Filename='HF_Sys_65x65_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_65x65_fixed';
% Filename='HF_Sys_45x80';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_45x80';
% Filename='HF_Sys_75x150';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_75x150';
% Filename='HF_Sys_77x77_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_77x77_fixed';
% Filename='HF_Sys_77x77_InSitu_5_4.5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_77x77_InSitu_5_4.5';
% Filename='HF_Sys_99x99_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_99x99_fixed';
% Filename='HF_Sys_100x180';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_100x180';
% Filename='HF_Sys_25x25_Geostress';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_25x25_Geostress';
% Filename='HF_Sys_33x33';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_33x33';
% Filename='HF_Sys_33x33_Static';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_33x33_Static';
% Filename='HF_Sys_79x49';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_79x49';
% Filename='HF_Sys_100x100';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_100x100';
% Filename='HF_Sys_Youn';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_Youn';
% Filename='Excav_Sys_25x25_InSitu_4_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Excav_Sys_25x25_InSitu_4_5'; %开挖(生死单元)
% Filename='Excav_25x25_InSitu_y_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Excav_25x25_InSitu_y_5'; %开挖(生死单元)
% Filename='n2xx';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\n2xx';
% Filename='n2xxf';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\n2xxf';
% Filename='HF_Sys_200x_300_80x80_FZ_5_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\HF_Sys_200x_300_80x80_FZ_5_5';
% Filename='contact_reduction_test_15x15';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\contact_reduction_test_15x15';
% Filename='contact_reduction_test_30x30';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\contact_reduction_test_30x30';
% Filename='contact_reduction_test_50x50';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\contact_reduction_test_50x50';
% Filename='Paper02_01';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_01';
% Filename='Paper02_03';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_03';
% Filename='Paper02_03_test1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_03_test1';
% Filename='XiaoLong_test_7x7';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\XiaoLong_test_7x7';
% Filename='Paper02_04';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04';
% Filename='Paper02_04_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fixed';
% Filename='Paper02_04_InSitu_5_4';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_InSitu_5_4';
% Filename='Paper02_04_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine';
% Filename='Paper02_04_fine_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_fixed';
% Filename='Paper02_04_fine_InSitu_5_4';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4';
% Filename='Paper02_04_fine_InSitu_5_4_test2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test2';
% Filename='Paper02_04_fine_InSitu_5_4_test3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test3';
% Filename='Paper02_04_fine_InSitu_5_4_test4';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test4';
% Filename='Paper02_04_fine_InSitu_5_4_test5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test5';
% Filename='Paper02_04_fine_InSitu_5_4_test6';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test6';
% Filename='Paper02_04_fine_InSitu_5_4_test7';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test7';
% Filename='Paper02_04_fine_InSitu_5_4_test8';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test8';
% Filename='Paper02_04_fine_InSitu_5_4_test9';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4_test9';
% Filename='Paper02_04_fine_InSitu_5_4.5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_4.5';
% Filename='Paper02_04_fine_InSitu_5_5';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_04_fine_InSitu_5_5';
% Filename='Paper02_new_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_new_fine';
% Filename='Paper02_new_fine_InSitu_5_4';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper02_new_fine_InSitu_5_4';

% Filename='Paper03_Verification';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper03_Verification';
% Filename='Paper03_Sensitivity_Basecase_HF';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper03_Sensitivity_Basecase_HF';
% Filename='Paper03_Sensitivity_Basecase_HF2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper03_Sensitivity_Basecase_HF2';
% Filename='Paper03_Sensitivity_Basecase_HF3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper03_Sensitivity_Basecase_HF3';
% Filename='Paper03_Sensitivity_Basecase_HF3';Work_Dirctory='X:\PhiPsi fortran work_c\PhiPsi_work\Paper03_Sensitivity_Basecase_HF3';
% Filename='Fatigue_Compression';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Fatigue_Compression';
% Filename='exa_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\exa_tension';
% Filename='Hole';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Hole';
% Filename='Holes_and_cracks';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Holes_and_cracks';
% Filename='exa_tension';Work_Dirctory='X:\PhiPsi_work\exa_tension';
% Filename='openmp_test_100x100';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\openmp_test_100x100';
% ---------------------以下是场问题相关-------------------------
% Filename='Static_Field_test_1_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Static_Field_test_1_fine';
% Filename='Static_Field_shale_gas';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Static_Field_shale_gas';
% Filename='Transient_Field_shale_gas';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_gas';
% Filename='Transient_Field_shale_gas2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_gas2';
% Filename='Transient_Field_shale_gas3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_gas3';
% Filename='Transient_Field_book_analytical';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_book_analytical';
% Filename='Transient_Field_shale_paper04';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_paper04';
% Filename='Transient_Field_shale_paper04_so';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_paper04';
% Filename='Transient_Field_shale_paper04-1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_paper04-1';
% Filename='Transient_Field_shale_paper04-2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_paper04-2';
% Filename='Transient_Field_shale_paper04-3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_shale_paper04-3';
% Filename='Transient_Field_thermal';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Transient_Field_thermal';

% ---------------------以下是复合材料问题相关-------------------------
% Filename='compos_16mmx4mm_all_fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_16mmx4mm_all_fixed';
% Filename='compos_16mmx4mm_free';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_16mmx4mm_free';
% Filename='compos_16mmx4mm_uni_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_16mmx4mm_uni_tension';
% Filename='compos_16mmx4mm_uni_tens_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_16mmx4mm_uni_tens_fine';  %1600单元
% Filename='compos_16mmx4mm_uni_tens_fine2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_16mmx4mm_uni_tens_fine2';  %4556单元
% Filename='compos_10mmx10mm_uni_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_10mmx10mm_uni_tension';
% Filename='compos_10mmx10mm_free';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_10mmx10mm_free';
% Filename='compos_10mmx10mm_uni_tension_fin';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_10mmx10mm_uni_tension_fin';
% Filename='compos_1mx1m_uni_tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_1mx1m_uni_tension';
% Filename='compos_4mmx8mm_uni_tens';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_4mmx8mm_uni_tens';
% Filename='compos_1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_1';
% Filename='compos_2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_2';
% Filename='compos_3';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\compos_3';
% ---------------------以下是动态分析问题相关-------------------------
% Filename='beam';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\beam';
% Filename='Dynamic_3_points_beam_Initial_V';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dynamic_3_points_beam_Initial_V';
% Filename='Dynamic_3_points_beam_Initial_V_cor';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dynamic_3_points_beam_Initial_V_cor';
% Filename='EQ_dam_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\EQ_dam_test';
% Filename='EQ_dam_test_fine1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\EQ_dam_test_fine1';
% Filename='Dynamic';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dynamic';
% Filename='Dynamic_Coarse';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dynamic_Coarse';
% ---------------------以下是PhiPsi网站算例相关-------------------------
% Filename='exa_inclusions';Work_Dirctory='X:\PhiPsi_work\exa_inclusions';
% Filename='exa_tension';Work_Dirctory='X:\PhiPsi_work\exa_tension';
% Filename='exa_earthquake';Work_Dirctory='X:\PhiPsi_work\exa_earthquake';
% Filename='exa_hydraulic_fracturing';Work_Dirctory='X:\PhiPsi work\exa_hydraulic_fracturing';
% Filename='exa_cohesive_crack';Work_Dirctory='X:\PhiPsi_work\exa_cohesive_crack';
% Filename='exa_multi_cracks';Work_Dirctory='X:\PhiPsi_work\exa_multi_cracks';
% Filename='exa_crossing_cracks';Work_Dirctory='X:\PhiPsi_work\exa_crossing_cracks';
% Filename='exa_cracks_meet_holes';Work_Dirctory='X:\PhiPsi_work\exa_cracks_meet_holes';
% Filename='exa_random_multi_cracks';Work_Dirctory='X:\PhiPsi_work\exa_random_multi_cracks';
% Filename='exa_plastic';Work_Dirctory='X:\PhiPsi_work\exa_plastic';
% Filename='Crack_Meet_Hole_ff';Work_Dirctory='X:\PhiPsi_work\Crack_Meet_Hole_ff';
% Filename='Crack_Meet_Hole_ff';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Crack_Meet_Hole_ff';
% Filename='exa_cr_emerge_from_ran_hole';Work_Dirctory='X:\PhiPsi_work\exa_cr_emerge_from_ran_hole';
% Filename='exa_cr_emerge_from_ran_hole2';Work_Dirctory='X:\PhiPsi_work\exa_cr_emerge_from_ran_hole2';
% Filename='exa_compression_shear';Work_Dirctory='X:\PhiPsi Work\exa_compression_shear';
% ---------------------以下是书上的算例相关-------------------------
% Filename='T_Shaped_Crack';Work_Dirctory='X:\PhiPsi_work\T_Shaped_Crack';
% Filename='Crossing_Crack';Work_Dirctory='X:\PhiPsi_work\Crossing_Crack';
% Filename='Snow_Shaped_Crack';Work_Dirctory='X:\PhiPsi_work\Snow_Shaped_Crack';
% Filename='Crack_Meet_Hole';Work_Dirctory='X:\PhiPsi_work\Crack_Meet_Hole';
% Filename='Hole';Work_Dirctory='X:\PhiPsi_work\Hole';
% Filename='Inclusion';Work_Dirctory='X:\PhiPsi_work\Inclusion';
% Filename='FEM';Work_Dirctory='X:\PhiPsi_work\FEM';
% Filename='XFEM';Work_Dirctory='X:\PhiPsi_work\XFEM';
% Filename='Center_Crack';Work_Dirctory='X:\PhiPsi_work\Center_Crack';
% Filename='Edge_Crack';Work_Dirctory='X:\PhiPsi_work\Edge_Crack';
% Filename='Propagating_Crack';Work_Dirctory='X:\PhiPsi_work\Propagating_Crack';
% Filename='Multi_Crack';Work_Dirctory='X:\PhiPsi_work\Multi_Crack';
% Filename='Dynamic';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dynamic';
% Filename='Dynamic_and_Propagation';Work_Dirctory='X:\PhiPsi_work\Dynamic_and_Propagation';
% Filename='Plastic_Notched_Tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Plastic_Notched_Tension';
% Filename='Plastic_Tension';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Plastic_Tension';
% Filename='Plastic_Plane_Strain';Work_Dirctory='X:\PhiPsi_work\Plastic_Plane_Strain';
% Filename='Plastic_Strip_Footing';Work_Dirctory='X:\PhiPsi_work\Plastic_Strip_Footing';
% Filename='hydraulic_fracturing_T';Work_Dirctory='X:\PhiPsi work\hydraulic_fracturing_T';
% Filename='hydraulic_fracturing_T';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\hydraulic_fracturing_T';
% Filename='Contact';Work_Dirctory='X:\PhiPsi_work\Contact';
% Filename='Contact_and_Propagation';Work_Dirctory='X:\PhiPsi_work\Contact_and_Propagation';
% Filename='Contact_and_Propagation';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Contact_and_Propagation';

% ---------------------分子动力学模拟相关-------------------------
% Filename='MD_test1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\MD_test1';
% Filename='MD_test';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\MD_test1';
% ---------------------近场动力学模拟相关-------------------------
% Filename='PD_test1';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\PD_test1';
% ---------------------塑性相关-------------------------
% Filename='notch_tension_plastic';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\notch_tension_plastic';
% Filename='plate_with_a_hole';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\plate_with_a_hole';
% Filename='Plastic_Plane_Strain';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Plastic_Plane_Strain';
% Filename='notch_tension_plastic_cor';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\notch_tension_plastic_cor';
% ---------------------新-------------------------
% Filename='model_feikang_slope';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\model_feikang_slope';
% Filename='Contact';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Contact';
% Filename='Three_point_beam_5x1.2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Three_point_beam_5x1.2';
% Filename='Three_point_beam_5x1.2_fine';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Three_point_beam_5x1.2_fine';
% Filename='Three_point_beam_5x1.2_finer';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Three_point_beam_5x1.2_finer';
% Filename='Three_point_beam_5x1.2_finer2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Three_point_beam_5x1.2_finer2'; %与Three_point_beam_5x1.2_finer完全相同
% Filename='Paper05_test1_100mmx100mm';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper05_test1_100mmx100mm';
% Filename='Paper05_test1_100mmx100mm_c';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper05_test1_100mmx100mm_c';
% Filename='Paper05_test1_100mmx100mm_cc';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Paper05_test1_100mmx100mm_cc';
% Filename='Dam_JiaXin_D';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\Dam_JiaXin_D';
Filename='paper8_2d_edge_crack';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_Work\paper8_2d_edge_crack';

Num_Step_to_Plot      = -999             ;%后处理结果计算步号或者破裂步号(对于水力压裂分析),若为-999,则全部绘制
                                        % 水力压裂分析不能设为-999
Defor_Factor          = 1000                ;%变形放大系数
Key_TipEnrich         = 1                ;%裂尖增强方案：1,标准；2,仅F1; 3,光滑过渡Heaviside; 4,粘聚裂尖
Key_HF_Analysis       = 1                ;%水力压裂分析
%-------------------------------
Key_Animation= [0 2 0 0 0];           % 1:displacement(2,Gauss);2:stress(2,Gauss);3:deformation;4:场变量(1:仅场变量,2:仅流量矢量,3:场变量和流量矢量);5:MD
Key_Ani_Ave  = [1 1 1 0 0];           % Key for uniform maxima and minima(1:displacement;2:stress;3:Blank;4:场变量(仅支持平均显示,即Key_Ani_Ave(4)==1));5:Blank 
Time_Delay   = 0.2;%0.2;                 % Delay time of each frame of the gif animation,default: 0.025
Key_Time_String = 1;                  %时间的单位: =1,则为s;=2,min;=3,hour;=4,day;=5,month;=6,year
%-------------------------------
% 第1行,有限元网格: Mesh(1),Node(2),El(3),Gauss poin ts(4),
%                   5: 裂缝及裂缝坐标点(=1,绘制裂缝;=2,绘制裂缝及坐标点),
%                   6: 计算点及其编号(=1,计算点;=2,计算点和编号),
%                   7: 裂缝节点(计算点)相关(=1,节点集增强节点载荷;=2,计算点净水压;=3,计算点流量;=4,计算点开度),
%                   增强节点(8),网格线(9),
%                   支撑剂(10),单元应力状态是否σ1-σ3>Tol(11),天然裂缝(12),单元接触状态(13),裂缝编号(14),Fracture zone(15)
% 第2行,网格变形图: Deformation(1),Node(2),El(3),Gauss points(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   FsBs(7=1or2or3),Deformed and undefor(8),Blank(9),支撑剂(10),Blank(11),Blank(12),
%                   单元接触状态(13),增强节点(14),Fracture zone(15)
% 第3行,应力云图:   Stress Contour(1,2:Gauss points),(2:1,Only Mises stress;2,仅x应力;3,仅y应力;4,仅剪应力),主应力(3:1,主应力;2,主应力+方向;3,仅最大主应力),
%                   Crack(5:1,line;2,shape),Scaling Factor(6),FsBs(7=1or2or3),
%                   undeformed or Deformed(8),mesh(9),支撑剂(10)Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% 第4行,位移云图:   Plot Dis-Contour(1,2:Gauss points),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   FsBs(7=1or2or3),,undeformed or Deformed(8),mesh(9),支撑剂(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% 第5行,场问题云图: Plot Contour(1:场值,2:流量矢量,3:场值+流量矢量),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   Blank(7),undeformed or Deformed(8),mesh(9),支撑剂(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% 第6行,分子动力学: Plot Contour(1:位置,2:速度,3:),轨迹(2=1,全部轨迹;>1,轨迹包括的计算步数),Blank(3),Blank(4),Blank(4),Scaling Factor(6),
%                   Blank(7),Original location(8),plot a mesh(9),Blank(10),Blank(11),Blank(12),Blank(13),Blank(14),Box(15)
%                         1   2   3   4   5   6              7   8   9  10  11  12  13  14   15
Key_PLOT(1,:)         = [ 1,  0,  0,  1,  2,  0,             4,  0,  1  ,1  ,0  ,1  ,1  ,0  ,1];  
Key_PLOT(2,:)         = [ 1,  0,  0,  0,  2,  Defor_Factor,  3,  1,  1  ,1  ,1  ,1  ,1  ,0  ,1];  
Key_PLOT(3,:)         = [ 2,  3,  0,  0,  2,  Defor_Factor,  0,  1,  0  ,0  ,0  ,1  ,0  ,0  ,1];  
Key_PLOT(4,:)         = [ 2,  0,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,1  ,0  ,1  ,0  ,0  ,1];
Key_PLOT(5,:)         = [ 0,  0,  0,  0,  1,  Defor_Factor,  0,  0,  0  ,1  ,0  ,1  ,0  ,0  ,1];   
Key_PLOT(6,:)         = [ 1, 50,  0,  0,  0,  1,             0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,1];  
%##########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_Post_2_Go
