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

function Plot_Deformation(DISP,FORCE_Matrix,Boundary_X_Matrix,Boundary_Y_Matrix,Post_Enriched_Nodes ...
                          ,isub,Crack_X,Crack_Y,POS,Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
					      Coors_Vertex,Coors_Junction,Coors_Tip,Crack_Tip_Type,Shaped_Crack_Points,Crushed_element)
% This function plots the initial geometry.

global Node_Coor Elem_Node Outline
global Num_Node Num_Elem
global Key_PLOT num_Crack Color_Mesh
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global Color_outline_Udefor num_of_Material
global aveg_area_ele Elem_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Defor_1 Color_Backgro_Defor_2 Color_Backgro_Defor_3 Color_Backgro_Defor_4
global Color_Backgro_Defor_5 Color_Backgro_Defor_6 Color_Backgro_Defor_7
global Color_Backgro_Defor_8 Color_Backgro_Defor_9 Color_Backgro_Defor_10
global Color_Crushed_ele
global Ele_Weibull_Par_TS Defor_Factor
global Na_Crack_X Na_Crack_Y Num_Step_to_Plot
global num_Na_Crack Key_HF_Analysis
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl
global num_Circ_Inclusion Circ_Inclu_Coor
global Color_Inclusion
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global num_Cross Cross_Coor Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Enriched_Node_Type_Incl Elem_Type_Incl POS_Incl
global Arc_Crack_Coor Yes_Arc_Crack
global Key_Data_Format
global num_HC HC_Coor Enriched_Node_Type_HC Node_HC_elem POS_HC Elem_Type_HC
global Killed_Elements num_Killed_Ele Ellipse_Hole_Coor num_Ellipse_Hole
global Title_Font Key_Figure_Control_Widget

disp('    > Plotting deformation....') 

scale = Key_PLOT(2,6);

% Get the new coordinates of all nodes.
New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);

% Get the maximum and minimum value of the new coordinates of all nodes.
Min_X_Coor_New = min(min(New_Node_Coor(:,1)));
Max_X_Coor_New = max(max(New_Node_Coor(:,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(:,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(:,2)));

xi = zeros(4,Num_Elem); yi = xi;
xi_o = zeros(4,Num_Elem); yi_o = xi_o;
xi_1 =[];yi_1 =[];
xi_2 =[];yi_2 =[];
xi_3 =[];yi_3 =[];
xi_4 =[];yi_4 =[];
xi_5 =[];yi_5 =[];
xi_6 =[];yi_6 =[];
xi_7 =[];yi_7 =[];
xi_8 =[];yi_8 =[];
xi_9 =[];yi_9 =[];
xi_10 =[];yi_10 =[];
for iElem = 1:Num_Elem
    NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
	      Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
	% xi(:,iElem) = New_Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
	% yi(:,iElem) = New_Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
	if Elem_Material(iElem)==1
		xi_1(:,iElem) = New_Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
		yi_1(:,iElem) = New_Node_Coor(NN',2);                                     % Initial y-coordinates of nodes		
	elseif Elem_Material(iElem)==2
		xi_2(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_2(:,iElem) = New_Node_Coor(NN',2);   
	elseif Elem_Material(iElem)==3
		xi_3(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_3(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==4
		xi_4(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_4(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==5
		xi_5(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_5(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==6
		xi_6(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_6(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==7
		xi_7(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_7(:,iElem) = New_Node_Coor(NN',2); 	
	elseif Elem_Material(iElem)==8
		xi_8(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_8(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==9
		xi_9(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_9(:,iElem) = New_Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==10
		xi_10(:,iElem) = New_Node_Coor(NN',1);                                     
		yi_10(:,iElem) = New_Node_Coor(NN',2); 		
	end
	if Key_PLOT(2,8) == 1
		xi_o(:,iElem) = Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
		yi_o(:,iElem) = Node_Coor(NN',2);                                     % Initial y-coordinates of nodes	
	end
end

% New figure.
Tools_New_Figure

hold on;
title_defor = ['Deformation (',num2str(Defor_Factor),')'];
title('Deformation','FontName',Title_Font,'FontSize',Size_Font)
%title('Deformation',num2str(Defor_Factor),'FontName',Title_Font,'FontSize',Size_Font)
axis off; 
axis equal;
delta = sqrt(aveg_area_ele);
axis([Min_X_Coor_New-delta  Max_X_Coor_New+delta  Min_Y_Coor_New-delta  Max_Y_Coor_New+delta]); 
    
% Plot deformed mesh.
disp(['      ----- Plotting deformed mesh......'])
patch(xi_1,yi_1,Color_Backgro_Defor_1) 
patch(xi_2,yi_2,Color_Backgro_Defor_2)    
patch(xi_3,yi_3,Color_Backgro_Defor_3)   
patch(xi_4,yi_4,Color_Backgro_Defor_4)  
patch(xi_5,yi_5,Color_Backgro_Defor_5)  
patch(xi_6,yi_6,Color_Backgro_Defor_6)  
patch(xi_7,yi_7,Color_Backgro_Defor_7)   
patch(xi_8,yi_8,Color_Backgro_Defor_7)   
patch(xi_9,yi_9,Color_Backgro_Defor_7)   
patch(xi_10,yi_10,Color_Backgro_Defor_7)   


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%       单元接触状态
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if Key_PLOT(2,13)==1 
    Yes_Contact1=0;  %裂纹面闭合
	Yes_Contact2=0;  %裂纹面由支撑剂支撑
	%如果接触状态文件存在则：
	if exist([Full_Pathname,'.elcs_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting contact status of elements...'])
		ELCS= load([Full_Pathname,'.elcs_',num2str(Num_Step_to_Plot)]);
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
			if ELCS(iElem) == 1    %裂纹面闭合
				xi_Elcs1(:,iElem) = New_Node_Coor(NN',1);                                     
				yi_Elcs1(:,iElem) = New_Node_Coor(NN',2); 	
				Yes_Contact1=1;
			end
			if ELCS(iElem) == 2    %裂纹面由支撑剂支撑
				xi_Elcs2(:,iElem) = New_Node_Coor(NN',1);                                     
				yi_Elcs2(:,iElem) = New_Node_Coor(NN',2); 	
				Yes_Contact2=1;
			end
		end
		if Yes_Contact1==1
			patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])       %用粉红色标记该单元
		end
		if Yes_Contact2==1
			patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])         %用紫色标记该单元
		end
	end
	%如果粘聚裂缝粘聚单元状态文件存在则(2017-04-24)：
	if exist([Full_Pathname,'.elco_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting cohesive crack status of elements...'])
		ELCS= load([Full_Pathname,'.elco_',num2str(Num_Step_to_Plot)]);
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
			if ELCS(iElem) == 1   %粘聚单元
				xi_Elcs1(:,iElem) = New_Node_Coor(NN',1);                                     
				yi_Elcs1(:,iElem) = New_Node_Coor(NN',2); 	
				Yes_Contact1=1;
			end
			% if ELCS(iElem) == 2    %裂纹面由支撑剂支撑
				% xi_Elcs2(:,iElem) = Node_Coor(NN',1);                                     
				% yi_Elcs2(:,iElem) = Node_Coor(NN',2); 	
				% Yes_Contact2=1;
			% end
		end
		if Yes_Contact1==1
			patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])               %用粉红色标记该单元
		end
		if Yes_Contact2==1
			patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])         %用紫色标记该单元
		end
	end
end  

% Plot the node numbers.
if Key_PLOT(2,2) == 1
    disp(['      ----- Plotting node number.....'])
    for iNode = 1:Num_Node
        text(New_Node_Coor(iNode,1),New_Node_Coor(iNode,2),1,num2str(iNode),'FontName',Title_Font,'FontSize',Node_Fontsize,'color',Node_Fontcolor)
    end
end

% Plot the element numbers.
if Key_PLOT(2,3) == 1
    disp(['      ----- Plotting element number.....'])
    for iElem = 1:Num_Elem
        NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
	          Elem_Node(iElem,3) Elem_Node(iElem,4)];
        XN = New_Node_Coor(NN,1);
        YN = New_Node_Coor(NN,2);
        text(mean(XN),mean(YN),1,num2str(iElem),'FontName',Title_Font,'FontSize',Elem_Fontsize,'color',Elem_Fontcolor)
    end
end

% Plot forces.
if Key_PLOT(2,7) == 1 || Key_PLOT(2,7) == 3
    disp(['      ----- Plotting forces of nodes......'])
    Max_x_Force = max(abs(FORCE_Matrix(:,1)));
	Max_y_Force = max(abs(FORCE_Matrix(:,2)));
	Max_Force   = max(Max_x_Force,Max_y_Force);
	
	W = Max_X_Coor_New - Min_X_Coor_New;
	H = Max_Y_Coor_New - Min_Y_Coor_New;
	
	% length of force arrow
    % REMOVE:length_arrow = sqrt(max_area_ele);
	length_arrow = max(W,H)/15.0;          
	
	% Loop through each node.
	for i = 1:Num_Node
	    if FORCE_Matrix(i,3) ~=0           % If the nodes has force load, then:
			c_force_x   = FORCE_Matrix(i,1);
			c_force_y   = FORCE_Matrix(i,2);

			delta_L_x = c_force_x*length_arrow/Max_Force;
			delta_L_y = c_force_y*length_arrow/Max_Force;
			
			StartPoint = [New_Node_Coor(i,1)-delta_L_x   New_Node_Coor(i,2)-delta_L_y     0];
			EndPoint   = [New_Node_Coor(i,1)             New_Node_Coor(i,2)               0];

			line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red')
			
			% The length of the head of the arrow.
			length_arrow_head = length_arrow/3;
			
			% Plot the head of the arrow.
			theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
			theta_1 = pi/2 - theta - pi/3;
			delta_x = -length_arrow_head*cos(theta_1);
			delta_y =  length_arrow_head*sin(theta_1);
			line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red');
			theta_2 = 3*pi/2 - theta + pi/3;
			delta_x = -length_arrow_head*cos(theta_2);
			delta_y =  length_arrow_head*sin(theta_2);
			line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red');
		end	
	end
end

% Plot boundary conditions by darwing triangles.
if Key_PLOT(2,7) == 2 || Key_PLOT(2,7) == 3
    disp(['      ----- Plotting boundary conditions......'])
    W = Max_X_Coor_New - Min_X_Coor_New;
	H = Max_Y_Coor_New - Min_Y_Coor_New;
	% Length of the edge of the boundary triangle.
	% REMOVE:l_edge = max(W,H)/20.0; 
    l_edge = 0.5*sqrt(aveg_area_ele);
	
	% Loop through boundary_x nodes.
	for i = 1:size(Boundary_X_Matrix,1)
	    c_node   = Boundary_X_Matrix(i);
		x_c_node = 	New_Node_Coor(c_node,1);
		y_c_node =  New_Node_Coor(c_node,2);
		% The x location of the boundary triangle.	
	    tri_X = [x_c_node   x_c_node-sqrt(3)/2*l_edge  x_c_node-sqrt(3)/2*l_edge  x_c_node]; 
		% The y location of the boundary triangle.
		tri_Y = [y_c_node   y_c_node+1/2*l_edge        y_c_node-1/2*l_edge        y_c_node]; 
		% Plot the boundary triangle.
        plot(tri_X,tri_Y,'color','black')                                                                  
	end
	% Loop through boundary_y nodes.
	for i = 1:size(Boundary_Y_Matrix,1)
	    c_node   = Boundary_Y_Matrix(i);
		x_c_node = 	New_Node_Coor(c_node,1);
		y_c_node =  New_Node_Coor(c_node,2);
		% The x location of the boundary triangle.	
	    tri_X = [x_c_node   x_c_node-1/2*l_edge         x_c_node+1/2*l_edge        x_c_node]; 
		% The y location of the boundary triangle.
		tri_Y = [y_c_node   y_c_node-sqrt(3)/2*l_edge   y_c_node-sqrt(3)/2*l_edge  y_c_node]; 
		% Plot the boundary triangle.
        plot(tri_X,tri_Y,'color','black')                                                                  
	end
end

% Plot enriched nodes
if isempty(Post_Enriched_Nodes) ~= 1
    if Key_PLOT(2,14)==1
		disp(['      ----- Plotting enriched nodes......'])
		for i =1 :size(Post_Enriched_Nodes,2)
			for j =1:Num_Node
				x_node = New_Node_Coor(j,1);                                          
				y_node = New_Node_Coor(j,2);                                          
				if Post_Enriched_Nodes(j,i)==1     % Tip nodes
					plot(x_node,y_node,'bs')
					% plot(x_node,y_node,'bs','MarkerSize',6)
				elseif Post_Enriched_Nodes(j,i)==2 % Heaviside nodes
					plot(x_node,y_node,'bo')
					% plot(x_node,y_node,'bo','MarkerSize',6)
				elseif Post_Enriched_Nodes(j,i)==3 % Junction nodes
					plot(x_node,y_node,'ko','MarkerSize',10,'Color','blue')
				elseif Post_Enriched_Nodes(j,i)==6 % Junction nodes of crack and hole
					plot(x_node,y_node,'ko','MarkerSize',11,'Color','black')
				elseif Post_Enriched_Nodes(j,i)==4 % Cross nodes
					plot(x_node,y_node,'^','MarkerSize',12,'Color','blue')
				end
			end
		end
	end
end

% Plot enriched nodes of circle holes.
if num_Hole~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of hole...'])
	for i=1:num_Hole
	    for j =1:Num_Node
		    if Enriched_Node_Type_Hl(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',5,'Color','r')
				% plot(x_node,y_node,'bo','MarkerSize',6,'Color',[0/255,0/255,0/255])
			end
		end
		
	end
end

% Plot enriched nodes of ellipse holes.
if num_Ellipse_Hole~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of ellipse hole...'])
	for i=1:num_Ellipse_Hole
	    for j =1:Num_Node
		    if Enriched_Node_Type_Hl(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',5,'Color','r')
				% plot(x_node,y_node,'bo','MarkerSize',6,'Color',[0/255,0/255,0/255])
			end
		end
		
	end
end

% 绘制圆形夹杂的增强节点.
if num_Circ_Inclusion~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of circle inclusions...'])
	for i=1:num_Circ_Inclusion
	    for j =1:Num_Node
		    if Enriched_Node_Type_Incl(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',4,'Color',[255/255,99/255,71/255])
				% plot(x_node,y_node,'bo','MarkerSize',6,'Color',[0/255,0/255,0/255])
			end
		end
		
	end
end

% 绘制多边形夹杂的增强节点.
if num_Poly_Inclusion~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of polygon inclusions...'])
	for i=1:num_Poly_Inclusion
	    for j =1:Num_Node
		    if Enriched_Node_Type_Incl(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',6,'Color',[255/255,99/255,71/255])
			end
		end
		
	end
end

% Plot enriched nodes of crosses.
if num_Cross~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of cross...'])
	for i=1:num_Cross
	    for j =1:Num_Node
		    if Enriched_Node_Type_Cross(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',9,'Color','red')
			end
		end
		
	end
end

% Plot enriched nodes of crosses(五角星).
if num_HC~=0 && Key_PLOT(2,14) == 1
    disp(['      ----- Plotting enriched nodes of HC cross...'])
	for i=1:num_HC
	    for j =1:Num_Node
		    if Enriched_Node_Type_HC(j,i) ==1
		        x_node = New_Node_Coor(j,1);                                          
	            y_node = New_Node_Coor(j,2);   
				plot(x_node,y_node,'p','MarkerSize',18,'Color','m')
			end
		end
		
	end
end

% Plot cracks line if necessary
if Key_PLOT(2,5) == 1
    disp(['      ----- Plotting crack lines......'])
	if num_Crack(isub)~=0
		for i = 1:num_Crack(isub)
			nPt = size(Crack_X{i},2);
			for i_Seg = 1:nPt-1
				x = [Crack_X{i}(i_Seg) Crack_X{i}(i_Seg+1)];
				y = [Crack_Y{i}(i_Seg) Crack_Y{i}(i_Seg+1)];
				for j =1:2
					% Get the local coordinates of the points of the crack. 
					[Kesi,Yita] = Cal_KesiYita_by_Coors(x(j),y(j));
					% Get the element number which contains the points of the crack. 
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(j),y(j));
					% Calculate the displacement of the points of the crack. 
					N1  = Elem_Node(c_Elem_Num,1);                                                  
					N2  = Elem_Node(c_Elem_Num,2);                                                  
					N3  = Elem_Node(c_Elem_Num,3);                                                  
					N4  = Elem_Node(c_Elem_Num,4);                                                
					U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3)...
						 DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
					% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
					[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
					dis_x(j) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
					dis_y(j) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
				end
				last_x = [ x(1)+dis_x(1)*scale x(2)+dis_x(2)*scale];
				last_y = [ y(1)+dis_y(1)*scale y(2)+dis_y(2)*scale];
				%----------------------------
				%如果是弧形裂缝,2017-07-18
				%----------------------------
				%Arc_Crack_Coor:x,y,r,Radian_Start,Radian_End,Radian,Point_Start_x,Point_Start_y,Point_End_x,Point_End_y
			    if abs(sum(Arc_Crack_Coor(i,i_Seg,1:11)))>=1.0e-10 
					c_R=Arc_Crack_Coor(i,i_Seg,4);
					c_Direcction  =Arc_Crack_Coor(i,i_Seg,3);
					c_Radian_Start=Arc_Crack_Coor(i,i_Seg,5)*pi/180;
					c_Radian_End  =Arc_Crack_Coor(i,i_Seg,6)*pi/180;
					%**********************
					%   如果是逆时针圆弧
					%**********************
					if c_Direcction >0.5
						%若结束角度大于起始角度,则直接绘制圆弧
						if c_Radian_End>=c_Radian_Start 
							c_alpha=c_Radian_Start:pi/100:c_Radian_End;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(k),c_y(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(k),c_y(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x(k) = c_x(k)+dis_c_x(k)*scale;last_c_y(k) = c_y(k)+dis_c_y(k)*scale;
							end
							plot(last_c_x,last_c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						%若结束角度小于起始角度,则分成两部分绘制圆弧
						else
							%第1部分:Radian_Start到360
							c_alpha1=c_Radian_Start:pi/100:2*pi;
							c_x=c_R*cos(c_alpha1)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha1)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(k),c_y(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(k),c_y(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x(k) = c_x(k)+dis_c_x(k)*scale;last_c_y(k) = c_y(k)+dis_c_y(k)*scale;
							end
							plot(last_c_x,last_c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
							%第2部分:0到Radian_End
							c_alpha_2=0.0:pi/100:c_Radian_End;
							c_x_2=c_R*cos(c_alpha_2)+Arc_Crack_Coor(i,i_Seg,1);
							c_y_2=c_R*sin(c_alpha_2)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x_2,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x_2(k),c_y_2(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x_2(k),c_y_2(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x_2(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y_2(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x_2(k) = c_x_2(k)+dis_c_x_2(k)*scale;last_c_y_2(k) = c_y_2(k)+dis_c_y_2(k)*scale;
							end
							plot(last_c_x_2,last_c_y_2,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						end
					%**********************
					%   如果是顺时针圆弧
					%**********************
					elseif c_Direcction < (-0.5)
						%若结束角度大于起始角度,则分成两部分绘制圆弧
						if c_Radian_End>=c_Radian_Start 
							%第1部分:Radian_End到360
							c_alpha1=c_Radian_End:pi/50:2*pi;
							c_x=c_R*cos(c_alpha1)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha1)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(k),c_y(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(k),c_y(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x(k) = c_x(k)+dis_c_x(k)*scale;last_c_y(k) = c_y(k)+dis_c_y(k)*scale;
							end
							plot(last_c_x,last_c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
							%第2部分:0到Radian_Start
							c_alpha_2=0.0:pi/50:c_Radian_Start;
							c_x_2=c_R*cos(c_alpha_2)+Arc_Crack_Coor(i,i_Seg,1);
							c_y_2=c_R*sin(c_alpha_2)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x_2,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x_2(k),c_y_2(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x_2(k),c_y_2(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x_2(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y_2(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x_2(k) = c_x_2(k)+dis_c_x_2(k)*scale;last_c_y_2(k) = c_y_2(k)+dis_c_y_2(k)*scale;
							end
							plot(last_c_x_2,last_c_y_2,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						%若结束角度小于起始角度,则直接绘制圆弧
						else
                            c_alpha=c_Radian_End:pi/50:c_Radian_Start;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							for k=1:size(c_x,2)
								[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(k),c_y(k));
								[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(k),c_y(k));
								N1  = Elem_Node(c_Elem_Num,1);N2  = Elem_Node(c_Elem_Num,2);                                                  
								N3  = Elem_Node(c_Elem_Num,3);N4  = Elem_Node(c_Elem_Num,4);                                                
								U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3) DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
								[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
								dis_c_x(k) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
								dis_c_y(k) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
								last_c_x(k) = c_x(k)+dis_c_x(k)*scale;last_c_y(k) = c_y(k)+dis_c_y(k)*scale;
							end
							plot(last_c_x,last_c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						end
					end
				%如果是直线裂缝
				elseif abs(sum(Arc_Crack_Coor(i,i_Seg,1:11)))<1.0e-10 
                    plot(last_x,last_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   					
				end
				
			end

		end	
	end
end

% Plot holes.
if num_Hole ~=0
    disp(['      ----- Plotting hole...'])
	for i = 1:num_Hole
	    Coor_x  = Hole_Coor(i,1);
		Coor_y  = Hole_Coor(i,2);
		c_R  = Hole_Coor(i,3);
		num_fineness = 200;
		for j = 1:num_fineness+1
		    alpha = 2.0*pi/num_fineness*(j-1);
	     	x(j) = Coor_x + c_R*cos(alpha);
		    y(j) = Coor_y + c_R*sin(alpha);
			[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(j),y(j));
            %Sometimes the element does not exist, for example, only part of the hole locates inside the model	
			if c_Elem_Num~=0
			    [Kesi,Yita] = Cal_KesiYita_by_Coors(x(j),y(j));
                [dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
				% 测试代码											 
				% if j==96						
                    % disp(['      x(j),y(j):',num2str(x(j)*1000.0),num2str(y(j)*1000.0)]);
			    	% disp(['      c_Elem_Num:',num2str(c_Elem_Num)]);
                    % [dis_x(j),dis_y(j)] = Cal_Anypoint_Disp2(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
															 % ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  % Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 	
				% end
                %绘制编号txt
				% ttt_x =  x(j) + dis_x(j)*Key_PLOT(2,6);
				% ttt_y =  y(j) + dis_y(j)*Key_PLOT(2,6);
				% text(ttt_x+0.03*delta,ttt_y,1,num2str(j),...
		              % 'FontName',Title_Font,'FontSize',Node_Fontsize,'color',Node_Fontcolor)	
					  
			end
		end
		x_new = x + dis_x*Key_PLOT(2,6);
		y_new = y + dis_y*Key_PLOT(2,6);
		patch(x_new,y_new,'white','edgecolor','black','LineWidth',0.1)	
	end	
end


% Plot ellipse holes, deformation has not been considered (2020-08-10).
if num_Ellipse_Hole ~=0
    disp(['      ----- Plotting ellipse hole...'])
	for i_Hole = 1:num_Ellipse_Hole
	    Coor_x  = Ellipse_Hole_Coor(i_Hole,1);
		Coor_y  = Ellipse_Hole_Coor(i_Hole,2);
		a_Hole  = Ellipse_Hole_Coor(i_Hole,3);
		b_Hole  = Ellipse_Hole_Coor(i_Hole,4);
		c_Hole  = sqrt(a_Hole^2-b_Hole^2);
		theta_Hole  = Ellipse_Hole_Coor(i_Hole,5);
		% Plot the ellipse.
		hold on;
		syms x y;
		theta_Hole = theta_Hole*pi/180.0;
		f=(a_Hole^2-c_Hole^2*cos(theta_Hole)^2)*(x-Coor_x).^2+(a_Hole^2-c_Hole^2*sin(theta_Hole)^2)*(y-Coor_y).^2-c_Hole^2*sin(2*theta_Hole)*(x-Coor_x).*(y-Coor_y)-a_Hole^2*b_Hole^2;
		h=ezplot(f,[-1.5*a_Hole+Coor_x,1.5*a_Hole+Coor_x,-1.5*a_Hole+Coor_y,1.5*a_Hole+Coor_y]);
		h.Color = 'black';
		title('Deformation','FontName',Title_Font,'FontSize',Size_Font)		
	end	
end

% 绘制圆形夹杂.
if num_Circ_Inclusion ~=0
    disp(['      ----- Plotting circle inclusion...'])
	for i = 1:num_Circ_Inclusion
	    Coor_x  = Circ_Inclu_Coor(i,1);
		Coor_y  = Circ_Inclu_Coor(i,2);
		c_R  = Circ_Inclu_Coor(i,3);
		num_fineness = 100;
		for j = 1:num_fineness+1
		    alpha = 2*pi/num_fineness*(j-1);
	     	x(j) = Coor_x + c_R*cos(alpha);
		    y(j) = Coor_y + c_R*sin(alpha);
			[Kesi,Yita] = Cal_KesiYita_by_Coors(x(j),y(j));
			[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(j),y(j));
            [dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
		end
		x_new = x + dis_x*Key_PLOT(2,6);
		y_new = y + dis_y*Key_PLOT(2,6);
		% plot(x_new,y_new,'-')
		patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
	end	
end

% 绘制多边形夹杂,added on 2016-10-04
if num_Poly_Inclusion ~=0
    disp(['      ----- Plotting poly inclusion...'])
	for i = 1:num_Poly_Inclusion
		nEdge = size(Poly_Incl_Coor_x{i},2);
		Num_Diversion = 5;    %多边形的每条边拆分成5份
		k = 0;
		for iEdge = 1:nEdge
		    %获取边线的起点和终点
			if iEdge==nEdge
			    Line_Edge(1,1) = Poly_Incl_Coor_x{i}(iEdge); %边线的起点
				Line_Edge(1,2) = Poly_Incl_Coor_y{i}(iEdge);
				Line_Edge(2,1) = Poly_Incl_Coor_x{i}(1);     %边线的终点
				Line_Edge(2,2) = Poly_Incl_Coor_y{i}(1);
			else
			    Line_Edge(1,1) = Poly_Incl_Coor_x{i}(iEdge);   %边线的起点
				Line_Edge(1,2) = Poly_Incl_Coor_y{i}(iEdge);
				Line_Edge(2,1) = Poly_Incl_Coor_x{i}(iEdge+1); %边线的终点
				Line_Edge(2,2) = Poly_Incl_Coor_y{i}(iEdge+1);
			end
			% 等分点.
			a_x = Line_Edge(1,1);
			a_y = Line_Edge(1,2);
			b_x = Line_Edge(2,1);
			b_y = Line_Edge(2,2);
			%计算边线起点的位移
			k =k+1;
			x(k) = a_x;
		    y(k) = a_y;
			[Kesi,Yita] = Cal_KesiYita_by_Coors(x(k),y(k));
			[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(k),y(k));
            [dis_x(k),dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
			%计算等分点的位移
			for j =  1:Num_Diversion-1
			    k=k+1;
				x(k) = (j*b_x+(Num_Diversion-j)*a_x)/Num_Diversion;
				y(k) = (j*b_y+(Num_Diversion-j)*a_y)/Num_Diversion;
			    [Kesi,Yita] = Cal_KesiYita_by_Coors(x(k),y(k));
			    [c_Elem_Num] = Cal_Ele_Num_by_Coors(x(k),y(k));
                [dis_x(k),dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
			end												      
		end
		x_new = x + dis_x*Key_PLOT(2,6);
		y_new = y + dis_y*Key_PLOT(2,6);
		% plot(x_new,y_new,'-')
		patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
	end	
end

% Plot undeformed outline.
if Key_PLOT(2,8)==1
    disp(['      ----- Plotting undeformed mesh......'])
	%patch(xi_o,yi_o,'r','facealpha',0.01,'EdgeColor','r')  %透明度'facealpha'
	line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
		 [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1.5,'Color',Color_outline_Udefor) 		 
end  

% 绘制天然裂缝.
if Key_PLOT(2,12) == 1
	disp(['      ----- Plotting natural crack line...'])
	if isempty(Na_Crack_X)==0
		for tt_i = 1:num_Na_Crack
			nPt = size(Na_Crack_X{tt_i},2);
			for iPt = 2:nPt
				x = [Na_Crack_X{tt_i}(iPt-1) Na_Crack_X{tt_i}(iPt)];
				y = [Na_Crack_Y{tt_i}(iPt-1) Na_Crack_Y{tt_i}(iPt)]; 
				for jj =1:2
					% Get the local coordinates of the points of the crack. 
					[Kesi,Yita] = Cal_KesiYita_by_Coors(x(jj),y(jj));
					% Get the element number which contains the points of the crack. 
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(jj),y(jj));
					% Calculate the displacement of the points of the crack. 
					N1  = Elem_Node(c_Elem_Num,1);                                                  
					N2  = Elem_Node(c_Elem_Num,2);                                                  
					N3  = Elem_Node(c_Elem_Num,3);                                                  
					N4  = Elem_Node(c_Elem_Num,4);                                                
					U = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3)...
						 DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
					% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
					[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,[],[]);
					dis_x(jj) = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);  
					dis_y(jj) = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
				end
				
				last_x = [ x(1)+dis_x(1)*scale x(2)+dis_x(2)*scale];
				last_y = [ y(1)+dis_y(1)*scale y(2)+dis_y(2)*scale];
				
				plot(last_x,last_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   
				% plot(x,y,'--','Color',Color_Crack,'LineWidth',Width_Crack)  
			end
		end	
	end
end

%杀死的单元用白色
if num_Killed_Ele>0
	for iElem = 1:Num_Elem
		NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
			  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
		if ismember(iElem,Killed_Elements)
			xi_killed(:) = New_Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
			yi_killed(:) = New_Node_Coor(NN',2);  
			patch(xi_killed,yi_killed,'white') 
		end
	end 
end

% Plot Gauss points.
if Key_PLOT(2,4) == 1
    disp(['      ----- Plotting Gauss points......'])
	
	% 读取Gauss点的坐标
	if Key_Data_Format==1 
		Gauss_Coor = load([Full_Pathname,'.gcor_',num2str(isub)]);
	elseif Key_Data_Format==2  %Binary
		c_file = fopen([Full_Pathname,'.gcor_',num2str(isub)],'rb');
		[cc_tem_Gauss_cor,cc_count]   = fread(c_file,inf,'double');
		fclose(c_file);
		%转换成Matlab中的数据格式
		for ccc_i=1:cc_count/2
			Gauss_Coor(ccc_i,1) = ccc_i;
			Gauss_Coor(ccc_i,2) = cc_tem_Gauss_cor(ccc_i*2-1);
			Gauss_Coor(ccc_i,3) = cc_tem_Gauss_cor(ccc_i*2);
		end
	end	
	
	
	
	% 读取Gauss点位移.
	if Key_Data_Format==1 
        Gauss_Disp = load([Full_Pathname,'.disg_',num2str(isub)]);
	elseif Key_Data_Format==2  %Binary
		c_file = fopen([Full_Pathname,'.disg_',num2str(isub)],'rb');
		[cc_tem_Gauss_dis,cc_count]   = fread(c_file,inf,'double');
		fclose(c_file);
		%转换成Matlab中的数据格式
		for ccc_i=1:cc_count/2
			Gauss_Disp(ccc_i,1) = ccc_i;
			Gauss_Disp(ccc_i,2) = cc_tem_Gauss_dis(ccc_i*2-1);
			Gauss_Disp(ccc_i,3) = cc_tem_Gauss_dis(ccc_i*2);
		end
	end
			
	plot(Gauss_Coor(:,2)+scale*Gauss_Disp(:,2),Gauss_Coor(:,3)+scale*Gauss_Disp(:,3),...
	                            'bo','MarkerSize',1,'Color','black')
	clear Gauss_Coor;
	clear Gauss_Disp;
end

% Plot crushed elements.
if isempty(Crushed_element)==0
	xi_crushed =[];    yi_crushed =[];
	size(Crushed_element);
	for iElem = 1:Num_Elem
	    if Crushed_element(iElem)==1
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
			xi_crushed(:,iElem) = New_Node_Coor(NN',1);                                   % Deformed x-coordinates of nodes
			yi_crushed(:,iElem) = New_Node_Coor(NN',2);                                   % Deformed y-coordinates of nodes 
        end		
	end
    patch(xi_crushed,yi_crushed,Color_Crushed_ele) 
end

% 绘制单元应力状态(是否σ1-σ3>Tol),用于描述诱导裂缝体积及诱导效果
if Key_PLOT(2,11) == 1
	%如果单元应力状态文件存在：
	if exist([Full_Pathname,'.elss_',num2str(isub)], 'file') ==2 
	    disp(['      > Plotting induced volume......'])    
        %读取文件		
		Cracked_element=load([Full_Pathname,'.elss_',num2str(isub)]);
		xi_cracked =[];    yi_cracked =[];
		for iElem = 1:Num_Elem
			if Cracked_element(iElem)==1
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
				xi_cracked(:,iElem) = New_Node_Coor(NN',1);                                   % Deformed x-coordinates of nodes
				yi_cracked(:,iElem) = New_Node_Coor(NN',2);                                   % Deformed y-coordinates of nodes 
			end		
		end
		patch(xi_cracked,yi_cracked,'red') 
    end
end
% Plot Weibull.
% xi_crushed =[];    yi_crushed =[];
% for iElem = 1:Num_Elem
	% if any(Ele_Weibull_Par_TS(1,:)==iElem)
		% NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
			  % Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
		% xi_crushed(:,iElem) = New_Node_Coor(NN',1);                                   % Deformed x-coordinates of nodes
		% yi_crushed(:,iElem) = New_Node_Coor(NN',2);                                   % Deformed y-coordinates of nodes 
	% end		
% end
% patch(xi_crushed,yi_crushed,Color_Crushed_ele); 

% Plot shaped cracks if necessary
if Key_PLOT(2,5) == 2 & num_Crack(isub)~=0
    disp(['      ----- Plotting shaped cracks......'])
	Plot_Shaped_Cracks(Shaped_Crack_Points)
end

%绘制支撑剂的圆球
if Key_PLOT(2,10)==1 
	%如果支撑剂直径和坐标文件存在：
	if exist([Full_Pathname,'.epcr_',num2str(isub)], 'file') ==2 
	    disp(['      > Plotting proppant......'])    
        %读取文件		
		c_D_Coors=load([Full_Pathname,'.epcr_',num2str(isub)]);
        num_Proppant = size(c_D_Coors,1);
        for i=1:num_Proppant
			%绘制实心圆
			alpha=0:pi/20:2*pi;
			c_Elem_Num = c_D_Coors(i,1);
			R = c_D_Coors(i,2)/2*Key_PLOT(2,6);
			old_Coor_x = c_D_Coors(i,3);
			old_Coor_y = c_D_Coors(i,4);
			omega = c_D_Coors(i,5);  %对应裂纹片段的倾角
			l_elem= sqrt(aveg_area_ele);
			offset_delta = 0.001*l_elem;
			%原支撑剂所在点的上裂纹偏置微小量之后的点
			old_Coor_x_Up  = old_Coor_x - offset_delta*sin(omega);
			old_Coor_y_Up  = old_Coor_y + offset_delta*cos(omega);
			%原支撑剂所在点的上裂纹偏置微小量之后下偏置点
			old_Coor_x_Low = old_Coor_x + offset_delta*sin(omega);
			old_Coor_y_Low = old_Coor_y - offset_delta*cos(omega);
			%计算支撑剂中心坐标变形后的坐标
			[Kesi_Up,Yita_Up] = Cal_KesiYita_by_Coors(old_Coor_x_Up,old_Coor_y_Up);
            [dis_x_Up,dis_y_Up] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi_Up,Yita_Up...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
			[Kesi_Low,Yita_Low] = Cal_KesiYita_by_Coors(old_Coor_x_Low,old_Coor_y_Low);
            [dis_x_Low,dis_y_Low] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi_Low,Yita_Low...
															 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
															  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
			%真实的支撑剂所在点的位移											  
			dis_x = (dis_x_Up + dis_x_Low)/2;			
            dis_y = (dis_y_Up + dis_y_Low)/2;				
			Coor_x = old_Coor_x  + dis_x*Key_PLOT(2,6);
			Coor_y = old_Coor_y  + dis_y*Key_PLOT(2,6);
			
			x = Coor_x + R*cos(alpha);
			y = Coor_y + R*sin(alpha);
			plot(x,y,'-')
			%axis equal
			fill(x,y,[128/255,138/255,135/255])
		end
	end
end

%绘制破裂区
if Key_PLOT(2,15)==1 
	%如果定义了破裂区
	if Yes_has_FZ ==1
	    disp(['      ----- Plotting fracture zone......'])  
		for i_line = 1:4
		    if i_line ==1
				x = [frac_zone_min_x,frac_zone_max_x];
				y = [frac_zone_min_y,frac_zone_min_y];
		    elseif i_line ==2
				x = [frac_zone_max_x,frac_zone_max_x];
				y = [frac_zone_min_y,frac_zone_max_y];
		    elseif i_line ==3
				x = [frac_zone_max_x,frac_zone_min_x];
				y = [frac_zone_max_y,frac_zone_max_y];
		    elseif i_line ==4
				x = [frac_zone_min_x,frac_zone_min_x];
				y = [frac_zone_max_y,frac_zone_min_y];
			end
			for j =1:2
				%%% Get the local coordinates of the points of the crack. 
				[cKesi,cYita] = Cal_KesiYita_by_Coors(x(j),y(j));
				%%% Get the element number which contains the points of the crack. 
				[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(j),y(j));
				%%% Calculate the displacement of the points of the crack. 
				N1  = Elem_Node(c_Elem_Num,1);                                                  
				N2  = Elem_Node(c_Elem_Num,2);                                                  
				N3  = Elem_Node(c_Elem_Num,3);                                                  
				N4  = Elem_Node(c_Elem_Num,4);                                                
				cU = [DISP(N1,2) DISP(N1,3) DISP(N2,2) DISP(N2,3)...
					  DISP(N3,2) DISP(N3,3) DISP(N4,2) DISP(N4,3)];
				%%% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
				[N,~,~,~]  = Cal_N_dNdkesi_J_detJ(cKesi,cYita,[],[]);
				cdis_x(j) = cU(1)*N(1,1) + cU(3)*N(1,3) + cU(5)*N(1,5) + cU(7)*N(1,7);  
				cdis_y(j) = cU(2)*N(1,1) + cU(4)*N(1,3) + cU(6)*N(1,5) + cU(8)*N(1,7);  
			end
			
			last_x = [ x(1)+cdis_x(1)*scale x(2)+cdis_x(2)*scale];
			last_y = [ y(1)+cdis_y(1)*scale y(2)+cdis_y(2)*scale];
			
			plot(last_x,last_y,'w','LineWidth',1.5,'Color','green')   
		end
	end
end

% Active Figure control widget (2021-08-01)
% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
% Press q to exit.
% Press r (or double-click) to reset to the initial.
if Key_Figure_Control_Widget==1
    fcw(gca);
end

% Save pictures.
Save_Picture(c_figure,Full_Pathname,'defm')