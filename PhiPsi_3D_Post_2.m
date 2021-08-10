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
% Filename='exa_3D_crack';Work_Dirctory='x:\PhiPsi_work\exa_3D_crack';
% Filename='exa_3D_hollow_cylinder';Work_Dirctory='x:\PhiPsi_work\exa_3D_hollow_cylinder';
% Filename='3D_Block_11x11x11';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11';
% Filename='3D_Block_Tension';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension';
% Filename='3D_Block_21x21x21_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Geo_Insitu'; %21x21x21
% Filename='exa_3D_block_tension';Work_Dirctory='X:\PhiPsi Work\exa_3D_block_tension';
% Filename='3D_Block_Tension_Fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_Fine';
% Filename='3D_beam_sifs';Work_Dirctory='X:\PhiPsi Work\3D_beam_sifs';
% Filename='3D_Beam_ThreePointBending';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Beam_ThreePointBending';
Filename='paper8_3d_edge_crack';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\paper8_3d_edge_crack';

Num_Step_to_Plot      = -999             ;%后处理结果计算步号或者破裂步号(对于水力压裂分析),若为-999,则全部绘制
                                        % 水力压裂分析不能设为-999
Defor_Factor          = 1                ;%变形放大系数
Key_TipEnrich         = 1                ;%裂尖增强方案：1,标准；2,仅F1; 3,光滑过渡Heaviside; 4,粘聚裂尖
Key_HF_Analysis       = 1                ;%水力压裂分析
%-------------------------------
Key_Animation= [0 0 1 0 0];           % 1:displacement(2,Gauss);2:stress(2,Gauss);3:deformation;4:场变量(1:仅场变量,2:仅流量矢量,3:场变量和流量矢量);5:MD
Key_Ani_Ave  = [1 1 1 0 0];           % Key for uniform maxima and minima(1:displacement;2:stress;3:Blank;4:场变量(仅支持平均显示,即Key_Ani_Ave(4)==1));5:Blank 
Time_Delay   = 0.8;%0.2;                 % Delay time of each frame of the gif animation,default: 0.025
Key_Time_String = 1;                  %时间的单位: =1,则为s;=2,min;=3,hour;=4,day;=5,month;=6,year
%-------------------------------
% 第1行,有限元网格: Mesh(1),Node(2),El(3),Gauss points(4),
%                   5: 裂缝及裂缝坐标点(=1,绘制裂缝;=2,绘制裂缝及坐标点),
%                   6: 计算点及其编号(=1,计算点;=2,计算点和编号),
%                   7: 裂缝节点(计算点)相关(=1,节点集增强节点载荷;=2,计算点净水压;=3,计算点流量;=4,计算点开度),
%                   增强节点(8),网格线(9),
%                   支撑剂(10),单元应力状态是否σ1-σ3>Tol(11),天然裂缝(12),单元接触状态(13),裂缝编号(14),Fracture zone(15)
% 第2行,网格变形图: Deformation(1),Node(2),El(3),Gauss points(4),Crack(5:1,area;2,area with node),Scaling Factor(6),
%                   FsBs(7=1or2or3),Deformed and undefor(8),网格线(9),支撑剂(10),Blank(11),Blank(12),
%                   单元接触状态(13),增强节点(14),Fracture zone(15)
% 第3行,应力云图:   Stress Contour(1,2:Gauss points),(2:1,Only Mises stress;2,仅x应力;3,仅剪应力),主应力(3:1,主应力;2,主应力+方向;3,仅最大主应力),
%                   Crack(5:1,line;2,shape),Scaling Factor(6),
%                   undeformed or Deformed(8),mesh(9),支撑剂(10)Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% 第4行,位移云图:   Plot Dis-Contour(1,2:Gauss points),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   Blank(7),undeformed or Deformed(8),mesh(9),支撑剂(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% 第5行,场问题云图: Plot Dis-Contour(1),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   Blank(7),undeformed or Deformed(8),mesh(9),支撑剂(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
%                         1   2   3   4   5   6              7   8   9  10  11  12  13  14   15
Key_PLOT(1,:)         = [ 0,  0,  2,  0,  1,  0,             0,  0,  0  ,0  ,0  ,0  ,0  ,0  ,0];  
Key_PLOT(2,:)         = [ 1,  0,  0,  0,  1,  Defor_Factor,  0,  1,  0  ,0  ,0  ,0  ,0  ,0  ,0];  
Key_PLOT(3,:)         = [ 0,  0,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,0  ,0  ,1  ,0  ,0  ,1];  
Key_PLOT(4,:)         = [ 0,  0,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,1  ,0  ,1  ,0  ,0  ,1];
Key_PLOT(5,:)         = [ 0,  0,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,1  ,0  ,1  ,0  ,0  ,1]; 
%##########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_Post_2_Go_3D
