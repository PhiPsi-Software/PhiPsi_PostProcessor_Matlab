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
%--------------------- PhiPsi_Post_Plot1 ---------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
%clear all; close all; clc; format compact;  format long;


function PhiPsi_2D_Post_1_Function(Input_Defor_Factor,Input_Num_Step_to_Plot,Input_Filename,...
                                   Input_Work_Dirctory,Input_Key_PLOT, ...
                                   Input_KIKII_Crack,Input_KIKII_Tip,Input_Contour_String)
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
global Key_Gas_Prod_rate Key_Gas_Production Key_One_Node_Pres

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
Version='1.9.1';Date='September 14, 2018';

disp(['  PhiPsi Post Processor 1.'])  
disp([' -----------------------------------------------------------------------']) 
disp([' > RELEASE INFORMATION:                                                 ']) 
disp(['   PhiPsi Post Processor 1 is used for plotting deformed or undeformed  ']) 
disp(['   mesh, contours of displacements and stresses at specified substep.   ']) 
disp([' -----------------------------------------------------------------------']) 
disp([' > AUTHOR: Shi Fang, Huaiyin Institute of Technology                    ']) 
disp([' > WEBSITE: http://PhiPsi.top                                           ']) 
disp([' > EMAIL: shifang@ustc.edu.cn                                           ']) 
disp([' -----------------------------------------------------------------------']) 
disp(['  '])     
       
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
%PhiPsi_Color_and_Font_Settings                      
PhiPsi_Color_and_Font_Settings_Function(Input_Contour_String)   

% -------------------------------------
%           从GUI传入数据.      
% -------------------------------------   
Defor_Factor      = Input_Defor_Factor;
Num_Step_to_Plot  = Input_Num_Step_to_Plot;
Filename          = Input_Filename;
Work_Dirctory     = Input_Work_Dirctory;
Key_PLOT          = zeros(7,15);
Key_PLOT(1:7,1:15)= Input_Key_PLOT(1:7,1:15);

PhiPsi_Post_1_Go
