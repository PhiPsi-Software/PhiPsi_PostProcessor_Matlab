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

function Plot_Mesh3D_Slice(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
% This function plots the 3D slice of the mesh.

global Node_Coor Elem_Node Key_POST_HF
global Num_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Min_Z_Coor Max_Z_Coor
global Key_PLOT aveg_area_ele
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global num_Crack num_of_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Elem_Material Num_Step_to_Plot DISP Num_Foc_z
global Num_Foc_x Num_Foc_y Foc_x Foc_y Num_Foc_z Foc_z FORCE_Matrix
global Key_Data_Format
global Crack_node_X Crack_node_Y Crack_node_Z
global Crack_Ele_1 Crack_Ele_2 Crack_Ele_3
global Model_Outline Model_OutArea
global Num_Bou_x Num_Bou_y Num_Bou_z Bou_x Bou_y Bou_z
global Crack_node_local_X Crack_node_local_Y Crack_node_local_Z
global Crack_Node_in_Ele
global Crack3D_Vector_S1_X Crack3D_Vector_S1_Y Crack3D_Vector_S1_Z
global Crack3D_CalP_X Crack3D_CalP_Y Crack3D_CalP_Z
global Crack3D_Vertex_Vector_X_X 
global Crack3D_Vertex_Vector_X_Y  
global Crack3D_Vertex_Vector_X_Z 
global Crack3D_Vertex_Vector_Y_X  
global Crack3D_Vertex_Vector_Y_Y  
global Crack3D_Vertex_Vector_Y_Z  
global Crack3D_Vertex_Vector_Z_X  
global Crack3D_Vertex_Vector_Z_Y  
global Crack3D_Vertex_Vector_Z_Z 	
global Crack3D_Meshed_Outline
global Num_Accuracy_Contour Key_Contour_Metd
global volume_ele max_volume_ele min_volume_ele aveg_volume_ele 
global Num_Contourlevel Key_Flipped_Gray Color_Contourlevel
global Color_Mesh
global Crack3D_CalP_Aperture
global Cracks_FluidEle_Vector_3D_X
global Cracks_FluidEle_Vector_3D_Y
global Cracks_FluidEle_Vector_3D_Z
global Max_Aperture_of_each_Crack
global Title_Font Key_Figure_Control_Widget

disp('    > Plotting slice of mesh....') 

%需要绘制网格Key_PLOT(1,9) = 1
if Key_PLOT(1,9) ==0
    disp('      ***************************************************************') 
	disp('      WARNING :: Key_PLOT(1,9) is set to 1 for 3D mesh slice ploting!') 
	disp('                 Key_PLOT(1,9) = 1                                   ') 
    disp('      ***************************************************************') 
    Key_PLOT(1,9) = 1;
end

Slice_Type    = Key_PLOT(1, 1)-1;
Location_Plot = Key_PLOT(1,15);


% Get the maximum and minimum value of the new coordinates of all nodes.
Min_X_Coor = min(min(Node_Coor(1:Num_Node,1)));
Max_X_Coor = max(max(Node_Coor(1:Num_Node,1)));
Min_Y_Coor = min(min(Node_Coor(1:Num_Node,2)));
Max_Y_Coor = max(max(Node_Coor(1:Num_Node,2)));
Min_Z_Coor = min(min(Node_Coor(1:Num_Node,3)));
Max_Z_Coor = max(max(Node_Coor(1:Num_Node,3)));

c_X_Length = Max_X_Coor-Min_X_Coor;
c_Y_Length = Max_Y_Coor-Min_Y_Coor;
c_Z_Length = Max_Z_Coor-Min_Z_Coor;

% 根据切片类型不同进行设置待绘制的坐标
if Slice_Type==1   %x轴切片
   Ploted_Nodes = find(Node_Coor(1:Num_Node,1)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,2);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,3);  
elseif Slice_Type==2   %y轴切片
   Ploted_Nodes = find(Node_Coor(1:Num_Node,2)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,1);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,3);
elseif Slice_Type==3   %z轴切片
   Ploted_Nodes = find(Node_Coor(1:Num_Node,3)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,1);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = Node_Coor(Ploted_Nodes,2);  
end


% New figure.
Tools_New_Figure
hold on;


% Plot the deformed mesh.
if Key_PLOT(1,9)==1
	disp('      Plot deformed mesh....') 	
	for iElem = 1:Num_Elem
		NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
			  Elem_Node(iElem,3) Elem_Node(iElem,4) ...
			  Elem_Node(iElem,5) Elem_Node(iElem,6) ...
			  Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
			  
		%6个面循环
		for i_face =1:6
			if i_face==1
				NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
						Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)];    
			elseif i_face==2
				NN_4 = [Elem_Node(iElem,5) Elem_Node(iElem,6) ...
						Elem_Node(iElem,7) Elem_Node(iElem,8) Elem_Node(iElem,5)];     
			elseif i_face==3
				NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
						Elem_Node(iElem,6) Elem_Node(iElem,5) Elem_Node(iElem,1)];  
			elseif i_face==4
				NN_4 = [Elem_Node(iElem,2) Elem_Node(iElem,6) ...
						Elem_Node(iElem,7) Elem_Node(iElem,3) Elem_Node(iElem,2)];  
			elseif i_face==5
				NN_4 = [Elem_Node(iElem,3) Elem_Node(iElem,4) ...
						Elem_Node(iElem,8) Elem_Node(iElem,7) Elem_Node(iElem,3)];     
			elseif i_face==6
				NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,4) ...
						Elem_Node(iElem,8) Elem_Node(iElem,5) Elem_Node(iElem,1)];   
			end
			
			%假如当前面的4个节点在Slice上		
			if ismember(NN_4,Ploted_Nodes)
				if Slice_Type==1       %x轴切片
					xi = Node_Coor(NN_4',2);                                         
					yi = Node_Coor(NN_4',3);      
				elseif Slice_Type==2   %y轴切片
					xi = Node_Coor(NN_4',1);                                         
					yi = Node_Coor(NN_4',3);   
				elseif Slice_Type==3   %z轴切片
					xi = Node_Coor(NN_4',1);                                         
					yi = Node_Coor(NN_4',2);   
				end                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
				patch(xi,yi,Color_Backgro_Mesh_1) 
			end	
		end			
	end 
	if  Slice_Type==1
		title(['Mesh of x-axis slice (x= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font) 
	elseif  Slice_Type==2
		title(['Mesh of y-axis slice (y= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)
	elseif Slice_Type==3 
		title(['Mesh of z-axis slice (z= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)			
	end	
end
% Plot crack.
if Key_PLOT(1,5)~=0
  if isempty(Crack_X)==0
	disp('      Plot the crack....') 	
		for i_Crack = 1:num_Crack(isub)
			nnode = size(Crack_node_X{i_Crack},2);
			% for i_nnode=1:nnode
				% c_node_x = [Crack_node_X{i_Crack}(i_nnode)];
				% c_node_y = [Crack_node_Y{i_Crack}(i_nnode)];
				% c_node_z = [Crack_node_Z{i_Crack}(i_nnode)];
			% end
			
			nele = size(Crack_Ele_1{i_Crack},2);
			for i_nele=1:nele
				
 
				cr_node_1 = Crack_Ele_1{i_Crack}(i_nele);
				cr_node_2 = Crack_Ele_2{i_Crack}(i_nele);
				cr_node_3 = Crack_Ele_3{i_Crack}(i_nele);
				cr_el_line1_P1 = [Crack_node_X{i_Crack}(cr_node_1),Crack_node_Y{i_Crack}(cr_node_1),Crack_node_Z{i_Crack}(cr_node_1)];
				cr_el_line1_P2 = [Crack_node_X{i_Crack}(cr_node_2),Crack_node_Y{i_Crack}(cr_node_2),Crack_node_Z{i_Crack}(cr_node_2)];
				cr_el_line2_P1 = [Crack_node_X{i_Crack}(cr_node_2),Crack_node_Y{i_Crack}(cr_node_2),Crack_node_Z{i_Crack}(cr_node_2)];
				cr_el_line2_P2 = [Crack_node_X{i_Crack}(cr_node_3),Crack_node_Y{i_Crack}(cr_node_3),Crack_node_Z{i_Crack}(cr_node_3)];
				cr_el_line3_P1 = [Crack_node_X{i_Crack}(cr_node_3),Crack_node_Y{i_Crack}(cr_node_3),Crack_node_Z{i_Crack}(cr_node_3)];
				cr_el_line3_P2 = [Crack_node_X{i_Crack}(cr_node_1),Crack_node_Y{i_Crack}(cr_node_1),Crack_node_Z{i_Crack}(cr_node_1)];		
				%ooooooooooooooooooooooooooooooooooooooo
				%                x轴切片
				%ooooooooooooooooooooooooooooooooooooooo
				if  Slice_Type==1 
					num_intersection = 0;
					%裂缝面三角形单元的第1条边
					if  (Location_Plot >=cr_el_line1_P1(1) & Location_Plot <= cr_el_line1_P2(1)) || (Location_Plot >=cr_el_line1_P2(1) & Location_Plot <= cr_el_line1_P1(1))
						num_intersection = num_intersection +1;
						%........................................................................
						%计算交点坐标
						%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
						%........................................................................
						x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
						x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
						inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
						inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
						
					end		
					%裂缝面三角形单元的第2条边
					if  (Location_Plot >=cr_el_line2_P1(1) & Location_Plot <= cr_el_line2_P2(1)) || (Location_Plot >=cr_el_line2_P2(1) & Location_Plot <= cr_el_line2_P1(1))
						num_intersection = num_intersection +1;
						x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
						x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
						inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
						inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
					end	
					%裂缝面三角形单元的第3条边
					if  (Location_Plot >=cr_el_line3_P1(1) & Location_Plot <= cr_el_line3_P2(1)) || (Location_Plot >=cr_el_line3_P2(1) & Location_Plot <= cr_el_line3_P1(1))
						num_intersection = num_intersection +1;
						x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
						x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
						inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
						inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
					end	
				%ooooooooooooooooooooooooooooooooooooooo
				%                y轴切片
				%ooooooooooooooooooooooooooooooooooooooo
				elseif  Slice_Type==2 %y轴切片
					num_intersection = 0;
					%裂缝面三角形单元的第1条边
					if  (Location_Plot >=cr_el_line1_P1(2) & Location_Plot <= cr_el_line1_P2(2)) || (Location_Plot >=cr_el_line1_P2(2) & Location_Plot <= cr_el_line1_P1(2))
						num_intersection = num_intersection +1;
						%........................................................................
						%计算交点坐标
						%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
						%........................................................................
						x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
						x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
						inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
						inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_z];
						
					end		
					%裂缝面三角形单元的第2条边
					if  (Location_Plot >=cr_el_line2_P1(2) & Location_Plot <= cr_el_line2_P2(2)) || (Location_Plot >=cr_el_line2_P2(2) & Location_Plot <= cr_el_line2_P1(2))
						num_intersection = num_intersection +1;
						x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
						x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
						inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
						inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_z];
					end	
					%裂缝面三角形单元的第3条边
					if  (Location_Plot >=cr_el_line3_P1(2) & Location_Plot <= cr_el_line3_P2(2)) || (Location_Plot >=cr_el_line3_P2(2) & Location_Plot <= cr_el_line3_P1(2))
						num_intersection = num_intersection +1;
						x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
						x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
						inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
						inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_z];
					end	
				%ooooooooooooooooooooooooooooooooooooooo
				%                z轴切片
				%ooooooooooooooooooooooooooooooooooooooo						
				elseif Slice_Type==3  %z轴切片
					num_intersection = 0;
					%裂缝面三角形单元的第1条边
					if  (Location_Plot >=cr_el_line1_P1(3) & Location_Plot <= cr_el_line1_P2(3)) || (Location_Plot >=cr_el_line1_P2(3) & Location_Plot <= cr_el_line1_P1(3))
						num_intersection = num_intersection +1;
						%........................................................................
						%计算交点坐标
						%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
						%........................................................................
						x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
						x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
						inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
						inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
						
					end		
					%裂缝面三角形单元的第2条边
					if  (Location_Plot >=cr_el_line2_P1(3) & Location_Plot <= cr_el_line2_P2(3)) || (Location_Plot >=cr_el_line2_P2(3) & Location_Plot <= cr_el_line2_P1(3))
						num_intersection = num_intersection +1;
						x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
						x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
						inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
						inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
					end	
					%裂缝面三角形单元的第3条边
					if  (Location_Plot >=cr_el_line3_P1(3) & Location_Plot <= cr_el_line3_P2(3)) || (Location_Plot >=cr_el_line3_P2(3) & Location_Plot <= cr_el_line3_P1(3))
						num_intersection = num_intersection +1;
						x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
						x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
						inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
						inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
						intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
					end						
				end
				%%/////////////////////////////////////////////////////
				%  如果有两个交点，绘制交点(暂时不支持变形后的裂缝)
				%//////////////////////////////////////////////////////
				if num_intersection==2
					plot(intersection_Points(1:2,1),intersection_Points(1:2,2),'LineWidth',2.0,'Color','black')
				end						
			end					
		end
  end
end


%绘制坐标轴
% h = Tools_mArrow3([0 0 0],[(c_X_Length+c_Y_Length+c_Z_Length)/15 0 0],'color','red','facealpha',1.0);
% ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
% h = Tools_mArrow3([0 0 0],[0 (c_X_Length+c_Y_Length+c_Z_Length)/15 0],'color','green','facealpha',1.0);
% ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
% h = Tools_mArrow3([0 0 0],[0 0 (c_X_Length+c_Y_Length+c_Z_Length)/15],'color','blue','facealpha',1.0);
% ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');


%绘制裂缝开度向量
if Key_PLOT(1,7)==4
	if isempty(Crack3D_CalP_Aperture)==0 && isempty(Cracks_FluidEle_Vector_3D_X)==0 && isempty(Crack3D_CalP_X)==0
		for i = 1:num_Crack(isub)
			max_length = (c_X_Length+c_Y_Length+c_Z_Length)/35;
			Max_Aperture =  max(Max_Aperture_of_each_Crack);
			nCalP = size(Crack3D_CalP_X{i},2);
			for j=1:nCalP
				c_CalP_x = Crack3D_CalP_X{i}(j);
				c_CalP_y = Crack3D_CalP_Y{i}(j);
				c_CalP_z = Crack3D_CalP_Z{i}(j);
				c_CalP_V_x = Cracks_FluidEle_Vector_3D_X{i}(j);
				c_CalP_V_y = Cracks_FluidEle_Vector_3D_Y{i}(j);
				c_CalP_V_z = Cracks_FluidEle_Vector_3D_Z{i}(j);
				c_Calp_Apeture = Crack3D_CalP_Aperture{i}(j);
				length_arrow = c_Calp_Apeture/Max_Aperture*max_length;
				c_delta_x = length_arrow*c_CalP_V_x;
				c_delta_y = length_arrow*c_CalP_V_y;
				c_delta_z = length_arrow*c_CalP_V_z;
				
				if Slice_Type==1  %x轴切片
				    if c_CalP_x ==Location_Plot
					    Arrow_Line_Width = 1.2;  %箭头线宽
					    %------------
						% 上箭头
						%------------
						StartPoint = [c_CalP_y           c_CalP_z          ];
						EndPoint   = [c_CalP_y+c_delta_y c_CalP_z+c_delta_z];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);	
					    %------------
						% 下箭头
						%------------
						StartPoint = [c_CalP_y           c_CalP_z          ];
						EndPoint   = [c_CalP_y-c_delta_y c_CalP_z-c_delta_z];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);						
					end
                elseif Slice_Type==2  %y轴切片
				    if c_CalP_y ==Location_Plot
					    Arrow_Line_Width = 1.2;  %箭头线宽
					    %------------
						% 上箭头
						%------------
						StartPoint = [c_CalP_x           c_CalP_z          ];
						EndPoint   = [c_CalP_x+c_delta_x c_CalP_z+c_delta_z];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);	
					    %------------
						% 下箭头
						%------------
						StartPoint = [c_CalP_x           c_CalP_z          ];
						EndPoint   = [c_CalP_x-c_delta_x c_CalP_z-c_delta_z];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);						
					end	
                elseif Slice_Type==3  %z轴切片
				    if c_CalP_z ==Location_Plot
					    Arrow_Line_Width = 1.2;  %箭头线宽
					    %------------
						% 上箭头
						%------------
						StartPoint = [c_CalP_x           c_CalP_y          ];
						EndPoint   = [c_CalP_x+c_delta_x c_CalP_y+c_delta_y];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);	
					    %------------
						% 下箭头
						%------------
						StartPoint = [c_CalP_x           c_CalP_y          ];
						EndPoint   = [c_CalP_x-c_delta_x c_CalP_y-c_delta_y];
						line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',Arrow_Line_Width)
						% The length of the head of the arrow.
						length_arrow_head = length_arrow/3.5;
						
						% Plot the head of the arrow.
						theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
						theta_1 = pi/2 - theta - pi/3;
						delta_x = -length_arrow_head*cos(theta_1);
						delta_y =  length_arrow_head*sin(theta_1);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);
						theta_2 = 3*pi/2 - theta + pi/3;
						delta_x = -length_arrow_head*cos(theta_2);
						delta_y =  length_arrow_head*sin(theta_2);
						line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',Arrow_Line_Width);						
					end										
				end

			end
		end
	end
end



% The range of the plot.
min_x = Min_X_Coor; max_x = Max_X_Coor;
min_y = Min_Y_Coor; max_y = Max_Y_Coor;	
min_z = Min_Z_Coor; max_z = Max_Z_Coor;

if Slice_Type==1       %x轴切片
	axis([min_y max_y min_z max_z]);	      
elseif Slice_Type==2   %y轴切片
	axis([min_x max_x min_z max_z]);	 
elseif Slice_Type==3   %z轴切片
	axis([min_x max_x min_y max_y]);	
end     

axis equal; 
%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
% colorbar('FontName',Title_Font,'FontSize',Size_Font);
set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')



% Active Figure control widget (2021-08-01)
% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
% Press q to exit.
% Press r (or double-click) to reset to the initial.
if Key_Figure_Control_Widget==1
    fcw(gca);
end


% Save pictures.
Save_Picture(c_figure,Full_Pathname,'mehs')



