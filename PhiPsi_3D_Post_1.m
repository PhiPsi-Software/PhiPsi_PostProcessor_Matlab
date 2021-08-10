

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
%--------------------- PhiPsi_Post_Plot ----------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
clear all; close all; clc; format compact;  format long;
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot
global Key_TipEnrich

% Number of Gauss points of enriched element (default 64) for integral solution 2.
Num_Gauss_Points = 64;       

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
% rotate3d on;
% Display welcome information.

set(gcf,'renderer','opengl')

Welcome  
       
tic;
Tclock=clock;
Tclock(1);

disp([' >> Start time is ',num2str(Tclock(2)),'/',num2str(Tclock(3)),'/',num2str(Tclock(1))...
     ,' ',num2str(Tclock(4)),':',num2str(Tclock(5)),':',num2str(round(Tclock(6))),'.'])
disp(' ') 

% Make the "patch" method supported by "getframe", added in version 4.8.10
% See more: http://www.mathworks.com/support/bugreports/384622
% opengl('software') 
opengl hardware     

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
% Filename='exa_3D_hollow_cylinder';Work_Dirctory='x:\PhiPsi work\exa_3D_hollow_cylinder';
% Filename='3D_Block_10x5x3';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_10x5x3';
% Filename='3D_Block_11x11x11';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11';
% Filename='3D_Block_11x11x11_no_force';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11_no_force';
% Filename='3D_Block_11x11x11_Pressure';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11_Pressure';
% Filename='3D_Block_23x23x5_Pressure';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_23x23x5_Pressure';
% Filename='3D_Block_23x23x9_Pressure';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_23x23x9_Pressure';
% Filename='3D_Block_23x23x23_Pressure';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_23x23x23_Pressure';
% Filename='3D_Block_23x23x23_Pressure_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_23x23x23_Pressure_Fixed';
% Filename='3D_Block_25x25x25_Pressure_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_25x25x25_Pressure_Fixed';
% Filename='3D_Block_29x29x29_Pressure_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_29x29x29_Pressure_Fixed';
% Filename='3D_Block_35x35x5_Pressure';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_35x35x5_Pressure';
% Filename='3D_Block_21x21x11_uneven_Pressur';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x11_uneven_Pressur';
% Filename='3D_Block_Paper06_Example01';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example01';
% Filename='3D_Block_Paper06_Example01_fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example01_fine';
% Filename='3D_Block_Paper06_Example03';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example03';
% Filename='3D_Block_Paper06_Example03-new';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example03-new-bak-gap=5.0 m-nodal-stress';
% Filename='exa_3D_block_tension';Work_Dirctory='X:\PhiPsi Work\exa_3D_block_tension';     
% Filename='exa_3D_block_tension_cir';Work_Dirctory='X:\PhiPsi Work\exa_3D_block_tension_cir';
% Filename='exa_3D_block_tension_cir';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\exa_3D_block_tension_cir';
% Filename='3D_Block_Tension';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension';  %17x17x17
% Filename='3D_Block_Tension_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_Coarse';  %11x11x1
% Filename='3D_Block_Tension_Coarse_fix_Bot';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_Coarse_fix_Bot';  %11x11x1
% Filename='3D_Block_Tension_Fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_Fine';   %21x21x21
% Filename='3D_Block_Tension_More_Fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_More_Fine'; %30x30x30
% Filename='3D_Block_Tension_More_Fine2';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_More_Fine2'; %35x35x35 = 4.2875万个单元(element size=1)
% Filename='3D_Block_Tension_More_Fine3';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_More_Fine3'; %35x35x35 = 12.5万个单元(element size=0.7)
% Filename='3D_Block_21x21x21_Geostatic_test';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Geostatic_test'; %21x21x21
% Filename='3D_Block_Compress_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Compress_Coarse';  %11x11x1
% Filename='3D_Block_Compress_x_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Compress_x_Coarse';  %11x11x1
% Filename='3D_Block_Compress_y_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Compress_y_Coarse';  %11x11x1
% Filename='3D_Block_11x11x11_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11_Geo_Insitu'; %11x11x11
% Filename='3D_Block_21x21x21_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Geo_Insitu'; %21x21x21
% Filename='3D_Block_21x21x21_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Geo_Insitu-2'; %21x21x21
% Filename='3D_Cylinder_Compression';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Cylinder_Compression'; 
% Filename='3D_Cylinder_Compression_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Cylinder_Compression_Coarse'; 
% Filename='3D_Block_31x31x31_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_31x31x31_Geo_Insitu'; %31x31x31
% Filename='3D_Block_33x33x33_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_33x33x33_Geo_Insitu'; %33x33x33
% Filename='3D_Block_35x35x35_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_35x35x35_Geo_Insitu'; %35x35x35
% Filename='3D_beam_sifs';Work_Dirctory='X:\PhiPsi Work\3D_beam_sifs';
% Filename='3D_Beam_ThreePointBending';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Beam_ThreePointBending';
% Filename='3D_Block_Paper06_Example03-new';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example03-new-bak-gap=5.0m-2';
% Filename='3D_Block_Paper06_Example03-new';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Paper06_Example03-new-bak-gap=5.0m-3';
% Filename='Song_3D_test01';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Song_3D_test01';
% Filename='Song_3D_test01';Work_Dirctory='x:\PhiPsi_Project-good\PhiPsi_work\Song_3D_test01';
% Filename='Song_3D_test01_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Song_3D_test01_Coarse';
% Filename='Song_3D_test01_More_Coarse';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Song_3D_test01_More_Coarse';
% Filename='Song_3D_test01_fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Song_3D_test01_fine';
% Filename='Song_3D_Small';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Song_3D_Small';
% Filename='Bend_3D_Skew_Crack';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Bend_3D_Skew_Crack';
% Filename='Bend_3D_Skew_Crack_Fine1';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Bend_3D_Skew_Crack_Fine1';
% Filename='Bend_3D_Skew_Crack_80';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\Bend_3D_Skew_Crack_80';
% Filename='3D_Tetrahedron_Elements';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Tetrahedron_Elements';
Filename='paper8_3d_edge_crack';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\paper8_3d_edge_crack';

Num_Step_to_Plot      =  -999                 ;%-999     %后处理结果计算步号
Defor_Factor          =  1              ;   %变形放大系数


% 第1行,有限元网格: Plot Mesh(=1),
%                   2: Node(=1,only nodes;=2,plot nodes + node number)
%                   3:Element(=1,plot solid element;=2,plot only element number;=3,plot element number + element center),Gauss points(4),
%                   5: 裂缝及裂缝坐标点(=1,绘制裂缝面;=2,绘制裂缝面及坐标点;=3,绘制裂缝面及坐标点编号),
%                   6: 计算点(流体节点)及其编号
%                      =1,流体单元和节点;          =2,流体单元计算点和编号;   =3,流体单元(彩色显式);
%                      =4,流体单元(彩色显式)及编号;=5,流体单元(彩色显式)及编号+计算点及编号),
%                   7: 流体单元计算点以及Gauss点相关信息
%                      =1, 节点集增强节点载荷(矢量);= 2, 计算点净水压(矢量); = 3, 计算点流量(矢量);  =4, 计算点开度(矢量);
%                      =12,流体单元净水压(云图);    =13,流体单元流量(云图);  =14,流体单元开度(云图),
%                      =21,流体单元Gauss点(单点积分)外法线向量;=22,流体单元Gauss点(单点积分)局部坐标系;=25,流体单元Gauss点(单点积分)接触力
%                   8:增强节点(=1,增强节点; ),9:网格线(9=1),
%                   10:绘制指定单元及其节点号(起始单元号),11:绘制指定单元及其节点号(结束单元号),
%                   天然裂缝(12),单元接触状态(13),裂缝编号(14),15(Blank)
% 第2行,网格变形图: Deformation(1),Node(2),El(3),Gauss points(4)
%                   5:Crack Shape. =1,fracture area; =2,fracture volume
%                   6:Scaling Factor.
%                   7:Boundary Conditions. =1 or =2 or =3. 
%                   8:变形前的网格. =1,变形前的全部网格；=2,变形前的网格边界
%                   9:网格线(=0,外边界线框;=1,全部网格;=2,增强单元网格+外边界线框;=3,表面网格+外边界线框;=4,增强单元网格+表面网格),
%                   10:=1,裂缝边界节点的局部坐标系;=2,裂尖增强单元的基准线(baseline);=3,裂尖增强单元的基准线(baseline)+基准局部坐标系;=4,裂尖增强节点对应的单元(参考单元)
%                   11:=1,裂缝面前缘节点主应力矢量;
%                   12:=1,绘制裂尖应力计算点,CFCP=2时有效;=2,绘制计算球内的高斯点,CFCP=2时有效;=3,绘制增强节点的位移矢量,
%                   13:单元接触状态
%                   14:增强节点
%                   15:Fracture zone
% 第3行,应力云图:   Nodal Stress Contour(1:x slice;2:y slice;3:z slice;=4:节点值),(2:=1,Only Mises stress;=2,仅x应力;=3,仅y应力;=4,仅z应力),
%                   主应力(3:1,主应力;2,主应力+方向;3,仅最大主应力),x or y or z slice location(4),
%                   Crack(5:=1,line;=2,shape;3,不显示增强单元),Scaling Factor(6),
%                   undeformed or Deformed(8),mesh(9),Blank(10)Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第4行,位移云图:   Nodal Displacement Contour(1:x slice;2:y slice;3:z slice;=4:节点值),(2=1,only x;=2,only y;=3,only z),Blank(3),Blank(4),
%                   Crack(5:=1,plot crack;3,不显示增强单元),Scaling Factor(6),
%                   Blank(7),Blank(8),mesh(9),Blank(10),Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第5行,裂缝云图:   Plot Crack Contour(=1,aperture contour;=2,),
%                   2:=1,results of fluid elements and nodes; =2, results of discrete crack surface
%                   3-Blank,4-Blank,
%                   5:编号绘制:=1,绘制裂缝面离散点编号(or流体节点编号); =2,绘制裂缝面离散单元编号(or流体单元编号); =3, nodes+elements
%                   6:Scaling Factor,7-Blank,8-Blank,
%                   9:网格线(=0,外边界线框;=1,全部网格;=2,增强单元网格+外边界线框;=3,表面网格+外边界线框;=4,增强单元网格+表面网格),Blank(10),
%                   Blank(11),Blank(12),Blank(13),Blank(14),Blank(15)
%                         1   2   3    4   5   6              7    8   9  10      11    12  13  14      15
Key_PLOT(1,:)         = [ 1,  1,  0,   0,  3,  0,             0,   1,  1  ,198   ,0   ,0  ,0  ,0  ,  8.0];  
Key_PLOT(2,:)         = [ 1,  0,  0,   0,  1,  Defor_Factor,  3,   2,  1  ,1     ,0     ,0  ,0  ,1  ,    0];  
Key_PLOT(3,:)         = [ 0,  3,  0,   0,  1,  Defor_Factor,  0,   1,  0  ,0     ,0     ,0  ,0  ,0  ,  0.5];  
Key_PLOT(4,:)         = [ 0,  2,  0,   0,  3,  Defor_Factor,  0,   0,  0  ,1     ,0     ,0  ,0  ,0  ,    5];
Key_PLOT(5,:)         = [ 0,  1,  0,   0,  3,  Defor_Factor,  0,   0,  0  ,0     ,0     ,0  ,0  ,0  ,    0];   

%###########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_Post_1_Go_3D
