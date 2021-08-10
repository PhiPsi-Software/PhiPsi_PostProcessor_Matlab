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

function Plot_CFD_Contours(isub)
% This function plots contours for CFD analysis.

global Mesh_Points Mesh_Connects Key_POST_HF
global Num_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor
global Key_PLOT aveg_area_ele Outline
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global num_Crack num_of_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Elem_Material Num_Step_to_Plot
global Na_Crack_X Na_Crack_Y num_Na_Crack
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl
global Color_Inclusion
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl Elem_Type_Incl POS_Incl
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Yes_Field_Problem
global Field_Boundary_Value
global Field_Boundary_Qn
global num_Cross Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Arc_Crack_Coor Yes_Arc_Crack
global Key_Data_Format
global num_HC HC_Coor Enriched_Node_Type_HC Node_HC_elem POS_HC Elem_Type_HC
global Mesh_Points number_mesh_points number_mesh_connect Mesh_Connects CFD_Contour_fine
global Key_Contour_Metd Num_Contourlevel Key_Flipped_Gray Color_Contourlevel Color_Mesh Inline
global Title_Font Key_Figure_Control_Widget

% Check if *.pres file exists or not.
Check_exist = [Full_Pathname,'.pres_',num2str(isub)];
if exist(Check_exist,'file') ==0
    disp(['    Error :: Can not find pres files.'])
	warndlg('Can not find pres files!','ERROR')
else
    Mesh_Points_Pres   = load(Check_exist);
end


disp('      Resample pressure....') 
%Key_Contour_Metd = 1
if 	Key_Contour_Metd==1
	[X,Y,c_Pressure] = griddata(Mesh_Points(:,2),Mesh_Points(:,3),Mesh_Points_Pres,...
			   unique(Mesh_Points(:,2)),unique(Mesh_Points(:,3))');
elseif 	Key_Contour_Metd==2
    % Get resample coors.
	delta = ((Max_X_Coor - Min_X_Coor)+(Max_Y_Coor-Min_Y_Coor))/2/300;
	gx    = Min_X_Coor:delta:Max_X_Coor; 
	gy    = Min_Y_Coor:delta:Max_Y_Coor;
	[c_Pressure,X,Y] = Tools_gridfit(Mesh_Points(:,2),Mesh_Points(:,3),Mesh_Points_Pres,gx,gy);
end	


			
% New figure.
Tools_New_Figure
hold on;
title('CFD pressure','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
delta = ((Max_X_Coor - Min_X_Coor)+(Max_Y_Coor-Min_Y_Coor))/2/CFD_Contour_fine;
axis([Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor]);
disp('      Contouring pressure....') 
% Set colormap.
if Key_Flipped_Gray==0
	colormap(Color_Contourlevel)
elseif Key_Flipped_Gray==1
	colormap(flipud(gray))
end
colorbar('FontName',Title_Font,'FontSize',Size_Font);
contourf(X,Y,c_Pressure,Num_Contourlevel,'LineStyle','none')

% Plot the outline of the mesh
line([Mesh_Points(Outline(:,1),2) Mesh_Points(Outline(:,2),2)]', ...
	 [Mesh_Points(Outline(:,1),3) Mesh_Points(Outline(:,2),3)]','LineWidth',0.5,'Color','black')		
%外边界填充
Tools_Fillout(Mesh_Points(Outline(:,1),2),Mesh_Points(Outline(:,1),3),[Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor],'w'); 
%内边界填充
if isempty(Inline)==0
    fill(Mesh_Points(Inline(:,1),2),Mesh_Points(Inline(:,1),3),'w');
end	
		
% Plot the mesh		
if Key_PLOT(2,3)==1		
	for iconnect =1:number_mesh_connect
		NN = [Mesh_Connects(iconnect,2) Mesh_Connects(iconnect,3) Mesh_Connects(iconnect,4) Mesh_Connects(iconnect,2)];   
		xi = Mesh_Points(NN',2);                                           % Deformed x-coordinates of nodes
		yi = Mesh_Points(NN',3);                                           % Deformed y-coordinates of nodes
		plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
	end
end


% Save pictures.
Save_Picture(c_figure,Full_Pathname,'pres')
