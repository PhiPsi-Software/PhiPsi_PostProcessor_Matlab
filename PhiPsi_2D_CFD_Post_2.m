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
global Time_Delay CFD_Contour_fine


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
Filename='CFD-fish-tail';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\CFD-fish-tail';
CFD_Contour_fine      = 200;          %云图的精细度,越大越精细,但越慢,默认200,取样精度:(W+L)/2/CFD_Contour_fine
Num_Step_to_Plot      = -999             ;%后处理结果计算步号或者破裂步号,若为-999,则全部绘制
%-------------------------------
Key_Animation= [1 0 0 0 0];           % 1:presure;2:stress(2,Gauss);3:deformation;4:场变量(1:仅场变量,2:仅流量矢量,3:场变量和流量矢量);5:MD
Key_Ani_Ave  = [0 0 0 0 0];           % Key for uniform maxima and minima(1:displacement;2:stress;3:Blank;4:场变量(仅支持平均显示,即Key_Ani_Ave(4)==1));5:Blank 
Time_Delay   = 0.025;%0.2;                 % Delay time of each frame of the gif animation,default: 0.025
Key_Time_String = 1;                  %时间的单位: =1,则为s;=2,min;=3,hour;=4,day;=5,month;=6,year

%-------------------------------
% 第1行,有限元网格: Mesh(1),Mesh points(3),Mesh points number(3),boundary condition points(4),blank
% 第2行,云图:   Contour(1),What to plot(2=1,pressure;=2),Plot mesh(3),blank
%                         1   2   3   4   5   6    7   8   9  10  11  12  13  14   15
Key_PLOT(1,:)         = [ 1,  0,  0,  1,  0,  0,   0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,0];  
Key_PLOT(2,:)         = [ 1,  1,  0,  0,  0,  0,   0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,0];   

%##########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_CFD_Post_2_Go
