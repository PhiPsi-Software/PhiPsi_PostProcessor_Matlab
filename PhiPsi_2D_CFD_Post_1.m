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
%--------------------- PhiPsi_Post_CFD_Plot1 ---------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
clear all; close all; clc; format compact;  format long;
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack Defor_Factor
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves num_Na_Crack
global Plot_Aperture_Curves Plot_Pressure_Curves Plot_Velocity_Curves Num_Step_to_Plot
global Plot_Tan_Aper_Curves
global Key_TipEnrich Plot_Quantity_Curves Plot_Concentr_Curves     
global num_Hole Plot_Wpnp_Curves Plot_Wphp_Curves
global num_Circ_Inclusion num_Poly_Inclusion
global Key_Gas_Prod_rate Key_Gas_Production Key_One_Node_Pres CFD_Contour_fine

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
Key_PLOT   = zeros(2,15);                                   % Initialize the Key_PLOT

%###########################################################################################################
%##########################            User defined part        ############################################
%###########################################################################################################
Filename='CFD-fish-tail';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\CFD-fish-tail';

CFD_Contour_fine      = 200;          %云图的精细度,越大越精细,但越慢,默认200,取样精度:(W+L)/2/CFD_Contour_fine
Num_Step_to_Plot      = -999        ;%后处理结果计算步号(若-999,则绘制最后一步的)


% 第1行,有限元网格: Mesh(1),Mesh points(3),Mesh points number(3),boundary condition points(4),blank
% 第2行,云图:   Contour(1),What to plot(2=1,pressure;=2),Plot mesh(3),blank
%                         1   2   3   4   5   6    7   8   9  10  11  12  13  14   15
Key_PLOT(1,:)         = [ 0,  0,  0,  1,  0,  0,   0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,0];  
Key_PLOT(2,:)         = [ 1,  1,  0,  0,  0,  0,   0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,0];   

%###########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_CFD_Post_1_Go
