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

function Plot_CFD_Mesh(isub)
% This function plots the initial geometry.

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
global Mesh_Points number_mesh_points number_mesh_connect Mesh_Connects
global Title_Font Key_Figure_Control_Widget

% Preparing.
xi = zeros(3,number_mesh_connect); yi = xi;
disp(['      ----- Plotting mesh......'])
for iconnect = 1:number_mesh_connect
    NN = [Mesh_Connects(iconnect,2) Mesh_Connects(iconnect,3) Mesh_Connects(iconnect,4) ];                             % Nodes for current element
	xi(:,iconnect) = Mesh_Points(NN',2);                                     % Initial x-coordinates of nodes
	yi(:,iconnect) = Mesh_Points(NN',3);                                     % Initial y-coordinates of nodes	
end

% New figure.
Tools_New_Figure
hold on;
title('CFD mesh','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
delta = ((Max_X_Coor - Min_X_Coor)+(Max_Y_Coor-Min_Y_Coor))/2/50;
axis([Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor]);

% Plot mesh.  
% patch(xi,yi,Color_Backgro_Mesh_1) 
patch(xi,yi,[255/255,250/255,240/255]); 

% Plot mesh points.
if Key_PLOT(1,2)==1  
	for j =1:number_mesh_points
		x_point= Mesh_Points(j,2);                                          
		y_point = Mesh_Points(j,3);                                          
		plot(x_point,y_point,'bo','Color','blue','MarkerSize',3)
	end
end
% Plot mesh points number.
if Key_PLOT(1,3)==1  
	for j =1:number_mesh_points
		x_point = Mesh_Points(j,2);                                          
		y_point = Mesh_Points(j,3);                                          
        text(Mesh_Points(j,2)+0.05*delta,Mesh_Points(j,3),1,num2str(j),...
		              'FontName',Title_Font,'FontSize',Node_Fontsize,'color',Node_Fontcolor)		
	end
end

% Plot boundary condition points.
if Key_PLOT(1,4)==1  
	for j =1:number_mesh_points
	  
		x_point= Mesh_Points(j,2);                                          
		y_point = Mesh_Points(j,3);      	
        if Mesh_Points(j,4) == 1		
		    plot(x_point,y_point,'.','Color','c','MarkerSize',12)  %青色
		elseif Mesh_Points(j,4) == 2
		    plot(x_point,y_point,'.','Color','y','MarkerSize',12)  %黄色
		elseif Mesh_Points(j,4) == 3
		    plot(x_point,y_point,'.','Color',[124/255,252/255,0/255],'MarkerSize',12)  %草绿色
		elseif Mesh_Points(j,4) == 4
		    plot(x_point,y_point,'.','Color',[1,0,0],'MarkerSize',12)  %大红色
		elseif Mesh_Points(j,4) == 5
		    plot(x_point,y_point,'.','Color','m','MarkerSize',12)	%洋红色		
		elseif Mesh_Points(j,4) == 6
		    plot(x_point,y_point,'.','Color','b','MarkerSize',12)	%蓝色		
        end			
	end
end

            
% Plot the node numbers.
% if Key_PLOT(1,1) ==1 && Key_PLOT(1,2) == 1
    % disp(['      ----- Plotting node number...'])
    % for iNode = 1:Num_Node
        % text(Mesh_Points(iNode,1)+0.05*delta,Mesh_Points(iNode,2),1,num2str(iNode),...
		              % 'FontName',Title_Font,'FontSize',Node_Fontsize,'color',Node_Fontcolor)
    % end
% end

% Plot the element numbers.
% if Key_PLOT(1,1) ==1 && Key_PLOT(1,3) == 1
    % disp(['      ----- Plotting element number...'])
    % for iconnect = 1:Num_Elem
        % NN = [Mesh_Connects(iconnect,1) Mesh_Connects(iconnect,2) ...
	          % Mesh_Connects(iconnect,3) Mesh_Connects(iconnect,4)];
        % XN = Mesh_Points(NN,1);
        % YN = Mesh_Points(NN,2);
        % text(mean(XN),mean(YN),1,num2str(iconnect),'FontName',Title_Font,'FontSize',Elem_Fontsize,'color',Elem_Fontcolor)
    % end
% end

% Save pictures.
Save_Picture(c_figure,Full_Pathname,'mesh')
