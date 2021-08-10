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

% This  file defines controlling parameters of PhiPsi.

global Size_Font Color_Mesh Color_Vector 
global Color_outline_Udefor Color_Backgro_Defor_1 Elem_Fontsize Node_Fontcolor
global Node_Fontsize Elem_Fontcolor Num_Contourlevel
global Color_Contourlevel
global Color_Crack Width_Crack Key_Contour_Metd
global Key_Figure_off
global Output_Freq
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Color_Backgro_Defor_1 Color_Backgro_Defor_2 Color_Backgro_Defor_3 Color_Backgro_Defor_4
global Color_Backgro_Defor_5 Color_Backgro_Defor_6 Color_Backgro_Defor_7
global Color_Backgro_Defor_8 Color_Backgro_Defor_9 Color_Backgro_Defor_10
global Color_Inclusion Key_Flipped_Gray
global Num_Accuracy_Contour Key_Posivite_Compressive
global Num_Colorbar_Level Colorbar_Font
global Title_Font Key_Figure_Control_Widget
 
Key_Figure_off        = 0;                                 % Key for close the graphics window.
Key_Contour_Metd      = 2;                                 % 1: Slow and fine(not available for animation); 2: fast but coarse(default).
Output_Freq           = 1;                                 % Frequency of output: every n iterations(includes the first) or 0: only last
Color_Crack           = 'black';                           % Color of cracks
Width_Crack           = 3.0 ;                              % Width of lines of cracks
Color_Mesh            = 'black';                           % Color of mesh of the contours and vectors display,default:'black'
Color_Vector          = 'blue';                            % Color of the vectors
% Color_Backgro_Mesh_1  = [0.92,0.92,0.92];                % Color of mesh background of material 1,[0.49,1.00,0.83]
Color_Backgro_Mesh_1  = [0.49,1.00,0.83];                  % Color of mesh background of material 1,[0.49,1.00,0.83]
% Color_Backgro_Mesh_1  = [30/255,144/255,1];
% Color_Backgro_Mesh_1  = [30/255,144/255,1];
% Color_Backgro_Mesh_1  = [245/255,245/255,245/255];

Color_Backgro_Mesh_2  = [0.94,0.50,0.83];                  % Color of mesh background of material 2
Color_Backgro_Mesh_3  = [0.29,0.30,0.53];                  % Color of mesh background of material 3
Color_Backgro_Mesh_4  = [0.39,0.70,0.43];                  % Color of mesh background of material 4
Color_Backgro_Mesh_5  = [0.66,0.58,0.50];                  % Color of mesh background of material 5
Color_Backgro_Mesh_6  = [0.59,1.00,0.23];                  % Color of mesh background of material 6
Color_Backgro_Mesh_7  = [0.69,1.00,0.13];                  % Color of mesh background of material 7
Color_Backgro_Mesh_8  = [0.69,0.80,0.13];                  % Color of mesh background of material 8
Color_Backgro_Mesh_9  = [0.69,0.70,0.13];                  % Color of mesh background of material 9
Color_Backgro_Mesh_10 = [0.69,0.60,0.13];                  % Color of mesh background of material 10

Color_Backgro_Defor_1 = [0.49,1.00,0.83];                  % Color of deformation background of material 1
Color_Backgro_Defor_2 = [0.94,0.50,0.83];                  % Color of deformation background of material 2
Color_Backgro_Defor_3 = [0.29,0.30,0.53];                  % Color of deformation background of material 3
Color_Backgro_Defor_4 = [0.39,0.70,0.43];                  % Color of deformation background of material 4
Color_Backgro_Defor_5 = [0.66,0.58,0.50];                  % Color of deformation background of material 5
Color_Backgro_Defor_6 = [0.59,1.00,0.23];                  % Color of deformation background of material 6
Color_Backgro_Defor_7 = [0.69,1.00,0.13];                  % Color of deformation background of material 7
Color_Backgro_Defor_8 = [0.69,0.80,0.13];                  % Color of deformation background of material 8
Color_Backgro_Defor_9 = [0.69,0.70,0.13];                  % Color of deformation background of material 9
Color_Backgro_Defor_10= [0.69,0.60,0.13];                  % Color of deformation background of material 10
Color_outline_Udefor  = 'red';                             % Color of undeformation background
Color_Crushed_ele     = 'red';        

Color_Inclusion       = [160/255,102/255,211/255'];        %Color of inclusion
%Color_Inclusion       = [227/255,23/255,13/255'];          %Color of inclusion

Node_Fontsize         = 8;                                 % Font Size of node number,Default:8
Elem_Fontsize         = 8;                                 % Font Size of element number display,Default:8
Node_Fontcolor        = 'black';                           % Font colour of node number display                                                                                            
Elem_Fontcolor        = 'black';                           % Font colour of element number display
Color_Contourlevel    = 'jet';                             % ColorContourLevel,'jet','hot','cool','gray','spring','summer'
% Color_Contourlevel    = 'summer';                          % ColorContourLevel,'jet','hot','cool','gray','spring','summer'
% Color_Contourlevel    = 'gray';                            % ColorContourLevel,'jet','hot','cool','gray','spring','summer'
Num_Colorbar_Level    = 8;                                 % Number of colorbar level of values, default: 8, 2021-07-27
Size_Font             = 10;                                % Front size, Default=12; 
Title_Font            = 'Helvetica';                       % Title Font,  'Helvetica','Times New Roman', 2021-07-29
Colorbar_Font         = 'Helvetica';                         % Colorbar front, 'Helvetica','Times New Roman'
Key_Flipped_Gray      = 0;                                 % Flipped gray contour level,more bigger,more deeper.
Num_Contourlevel      = 50;                                % Number of contour levels,Default:12  
% Num_Contourlevel      = 10                                 % Number of contour levels,Default:12  
Num_Accuracy_Contour  = 10;                                % Number of sample point every one element length for contouring. 
Key_Posivite_Compressive = 1;                              % if=1 then posivite value indicates compressive stress, for 2D stress contours only.  
Key_Figure_Control_Widget= 1;                              % Active figure control widget fcw.m (default 0) (2021-08-01).