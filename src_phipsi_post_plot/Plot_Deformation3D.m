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

function Plot_Deformation3D(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
%绘制三维变形图.

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
global Full_Pathname
global Tip_Enriched_Ele_BaseLine Tip_Enriched_Ele_BaseLine_Vector_x 
global Tip_Enriched_Ele_BaseLine_Vector_y Tip_Enriched_Ele_BaseLine_Vector_z
global Tip_Enriched_Node_Ref_Element
global Cracks_CalP_UpDis_3D_X Cracks_CalP_UpDis_3D_Y Cracks_CalP_UpDis_3D_Z
global Cracks_CalP_LowDis_3D_X Cracks_CalP_LowDis_3D_Y Cracks_CalP_LowDis_3D_Z
global Crack3D_FluEl_Aperture Crack3D_Fluid_Ele_Num Crack3D_Fluid_Ele_Nodes
global Cracks_FluNode_Vector_3D_X Cracks_FluNode_Vector_3D_Y Cracks_FluNode_Vector_3D_Z
global Title_Font Key_Figure_Control_Widget

% global Tri_BCD
disp(['      ----- Plotting deformed mesh......'])
xi_1 =[];yi_1 =[];zi_1 =[];
xi_2 =[];yi_2 =[];zi_2 =[];
xi_3 =[];yi_3 =[];zi_3 =[];
xi_4 =[];yi_4 =[];zi_4 =[];
xi_5 =[];yi_5 =[];zi_5 =[];
xi_6 =[];yi_6 =[];zi_6 =[];
xi_7 =[];yi_7 =[];zi_7 =[];
xi_8 =[];yi_8 =[];zi_8 =[];
xi_9 =[];yi_9 =[];zi_9 =[];
xi_10 =[];yi_10 =[];zi_10 =[];

scale = Key_PLOT(2,6);

% Get the new coordinates of all nodes.
New_Node_Coor(1:Num_Node,1) = Node_Coor(1:Num_Node,1) + scale*DISP(1:Num_Node,2);
New_Node_Coor(1:Num_Node,2) = Node_Coor(1:Num_Node,2) + scale*DISP(1:Num_Node,3);
New_Node_Coor(1:Num_Node,3) = Node_Coor(1:Num_Node,3) + scale*DISP(1:Num_Node,4);


% Get the maximum and minimum value of the new coordinates of all nodes.
Min_X_Coor_New = min(min(New_Node_Coor(1:Num_Node,1)));
Max_X_Coor_New = max(max(New_Node_Coor(1:Num_Node,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(1:Num_Node,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(1:Num_Node,2)));
Min_Z_Coor_New = min(min(New_Node_Coor(1:Num_Node,3)));
Max_Z_Coor_New = max(max(New_Node_Coor(1:Num_Node,3)));

c_X_Length = Max_X_Coor_New-Min_X_Coor_New;
c_Y_Length = Max_Y_Coor_New-Min_Y_Coor_New;
c_Z_Length = Max_Z_Coor_New-Min_Z_Coor_New;
	     
% New figure.
Tools_New_Figure
hold on;
title('Deformation','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
delta = sqrt(aveg_area_ele);
axis([Min_X_Coor_New-delta Max_X_Coor_New+delta ...
      Min_Y_Coor_New-delta Max_Y_Coor_New+delta ...
	  Min_Z_Coor_New-delta Max_Z_Coor_New+delta]);
view(-37.5,30)  %视角
% view(0,90)        %xy视角
% view(180,0)        %xz视角
% view(90,0)        %yz视角
% view(-7,-15)  %视角


%绘制裂尖应力计算点及计算球内的高斯点
if Key_PLOT(2,12)==1 |  Key_PLOT(2,12)==2   
    %绘制裂尖应力计算点
	if exist([Full_Pathname,'.vspt_',num2str(isub)], 'file') ==2
		test= load([Full_Pathname,'.vspt_',num2str(isub)]);
		if isempty(test)==0 %不为空
		    plot3(test(:,1),test(:,2),test(:,3),'k.','MarkerSize',15.0,'color','m') 
			num_Ball = size(test,1);
			%绘制透明圆球
			for i_Ball = 1:num_Ball
				r_ball = test(i_Ball,4);
				x_cor = test(i_Ball,1);
				y_cor = test(i_Ball,2);
				z_cor = test(i_Ball,3);
				[x_ball,y_ball,z_ball]=sphere(30);
				surf(r_ball*x_ball +x_cor,r_ball*y_ball +y_cor,r_ball*z_ball +z_cor,'FaceColor','red','LineStyle','none','FaceAlpha',0.2);
			end
		end
	end
	% 计算球内的高斯点
	if Key_PLOT(2,12)==2 
		if exist([Full_Pathname,'.vspg_',num2str(isub)], 'file') ==2
			test= load([Full_Pathname,'.vspg_',num2str(isub)]);
			num_Total_Gauss = size(test,1);
			if isempty(test)==0%不为空
			    plot3(test(:,1),test(:,2),test(:,3),'k.','MarkerSize',5.0) 
			end
		end	
	end
end

%绘制单元面
Color_3D_ele_face = [189/255,252/255,201/255];
FaceAlpha_3D_ele_face = 0.8;
if Key_PLOT(2,3) == 1
  c_plot_count = 0;
  Ploted_Ele_Face =[0, 0, 0, 0];
  c_x_1=[];c_y_1=[];c_z_1=[];
  c_x_2=[];c_y_2=[];c_z_2=[];
  c_x_3=[];c_y_3=[];c_z_3=[];
  c_x_4=[];c_y_4=[];c_z_4=[];
  c_x_5=[];c_y_5=[];c_z_5=[];
  c_x_6=[];c_y_6=[];c_z_6=[];
  count_1 = 0;
  count_2 = 0;
  count_3 = 0;
  count_4 = 0;
  count_5 = 0;
  count_6 = 0;
  for iElem = 1:Num_Elem
    NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
	      Elem_Node(iElem,3) Elem_Node(iElem,4) ...
		  Elem_Node(iElem,5) Elem_Node(iElem,6) ...
		  Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
	  %第1个面
	  if not(ismember(sort([NN(1),NN(2),NN(3),NN(4)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(1),NN(2),NN(3),NN(4)]);	
		count_1 = count_1 +1;
	    c_x_1(1:4,count_1) = [New_Node_Coor(NN(1),1),New_Node_Coor(NN(2),1),New_Node_Coor(NN(3),1),New_Node_Coor(NN(4),1)];
	    c_y_1(1:4,count_1) = [New_Node_Coor(NN(1),2),New_Node_Coor(NN(2),2),New_Node_Coor(NN(3),2),New_Node_Coor(NN(4),2)];
	    c_z_1(1:4,count_1) = [New_Node_Coor(NN(1),3),New_Node_Coor(NN(2),3),New_Node_Coor(NN(3),3),New_Node_Coor(NN(4),3)];		
	  end	  
	  %第2个面
	  if not(ismember(sort([NN(5),NN(6),NN(7),NN(8)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(5),NN(6),NN(7),NN(8)]);	
		count_2 = count_2 +1;
		c_x_2(1:4,count_2) = [New_Node_Coor(NN(5),1),New_Node_Coor(NN(6),1),New_Node_Coor(NN(7),1),New_Node_Coor(NN(8),1)];
		c_y_2(1:4,count_2) = [New_Node_Coor(NN(5),2),New_Node_Coor(NN(6),2),New_Node_Coor(NN(7),2),New_Node_Coor(NN(8),2)];
		c_z_2(1:4,count_2) = [New_Node_Coor(NN(5),3),New_Node_Coor(NN(6),3),New_Node_Coor(NN(7),3),New_Node_Coor(NN(8),3)];			
	  end	
	  %第3个面
	  if not(ismember(sort([NN(1),NN(2),NN(6),NN(5)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(1),NN(2),NN(6),NN(5)]);	
		count_3 = count_3 +1;
		c_x_3(1:4,count_3) = [New_Node_Coor(NN(1),1),New_Node_Coor(NN(2),1),New_Node_Coor(NN(6),1),New_Node_Coor(NN(5),1)];
		c_y_3(1:4,count_3) = [New_Node_Coor(NN(1),2),New_Node_Coor(NN(2),2),New_Node_Coor(NN(6),2),New_Node_Coor(NN(5),2)];
		c_z_3(1:4,count_3) = [New_Node_Coor(NN(1),3),New_Node_Coor(NN(2),3),New_Node_Coor(NN(6),3),New_Node_Coor(NN(5),3)];			
	  end	
	  %第4个面
	  if not(ismember(sort([NN(2),NN(3),NN(7),NN(6)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(2),NN(3),NN(7),NN(6)]);	
		count_4 = count_4 +1;
	    c_x_4(1:4,count_4) = [New_Node_Coor(NN(2),1),New_Node_Coor(NN(3),1),New_Node_Coor(NN(7),1),New_Node_Coor(NN(6),1)];
	    c_y_4(1:4,count_4) = [New_Node_Coor(NN(2),2),New_Node_Coor(NN(3),2),New_Node_Coor(NN(7),2),New_Node_Coor(NN(6),2)];
	    c_z_4(1:4,count_4) = [New_Node_Coor(NN(2),3),New_Node_Coor(NN(3),3),New_Node_Coor(NN(7),3),New_Node_Coor(NN(6),3)];			
	  end
	  %第5个面
	  if not(ismember(sort([NN(7),NN(8),NN(4),NN(3)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(7),NN(8),NN(4),NN(3)]);	
		count_5 = count_5 +1;
        c_x_5(1:4,count_5) = [New_Node_Coor(NN(7),1),New_Node_Coor(NN(8),1),New_Node_Coor(NN(4),1),New_Node_Coor(NN(3),1)];
        c_y_5(1:4,count_5) = [New_Node_Coor(NN(7),2),New_Node_Coor(NN(8),2),New_Node_Coor(NN(4),2),New_Node_Coor(NN(3),2)];
        c_z_5(1:4,count_5) = [New_Node_Coor(NN(7),3),New_Node_Coor(NN(8),3),New_Node_Coor(NN(4),3),New_Node_Coor(NN(3),3)];			
	  end
	  %第6个面
	  if not(ismember(sort([NN(5),NN(1),NN(4),NN(8)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(5),NN(1),NN(4),NN(8)]);	
		count_6 = count_6 +1;
	    c_x_6(1:4,count_6) = [New_Node_Coor(NN(5),1),New_Node_Coor(NN(1),1),New_Node_Coor(NN(4),1),New_Node_Coor(NN(8),1)];
	    c_y_6(1:4,count_6) = [New_Node_Coor(NN(5),2),New_Node_Coor(NN(1),2),New_Node_Coor(NN(4),2),New_Node_Coor(NN(8),2)];
	    c_z_6(1:4,count_6) = [New_Node_Coor(NN(5),3),New_Node_Coor(NN(1),3),New_Node_Coor(NN(4),3),New_Node_Coor(NN(8),3)];		
	  end
  end 
  patch(c_x_1(1:4,1:count_1),c_y_1(1:4,1:count_1),c_z_1(1:4,1:count_1),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)    %2021-08-02
  patch(c_x_2(1:4,1:count_2),c_y_2(1:4,1:count_2),c_z_2(1:4,1:count_2),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_3(1:4,1:count_3),c_y_3(1:4,1:count_3),c_z_3(1:4,1:count_3),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_4(1:4,1:count_4),c_y_4(1:4,1:count_4),c_z_4(1:4,1:count_4),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_5(1:4,1:count_5),c_y_5(1:4,1:count_5),c_z_5(1:4,1:count_5),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_6(1:4,1:count_6),c_y_6(1:4,1:count_6),c_z_6(1:4,1:count_6),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)    
end




%绘制变形后的网格
Line_width =0.1;
if Key_PLOT(2,1)==1
    if Key_PLOT(2,9)==1 &&  Key_PLOT(2,3) ~= 1
		c_plot_count = 0;
		Ploted_Ele_lines =[0 0];
		to_be_plot_count = 0;
		to_be_plot_x = [];to_be_plot_y = [];to_be_plot_z = [];  %2021-08-02
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) Elem_Node(iElem,3) Elem_Node(iElem,4) ...
				  Elem_Node(iElem,5) Elem_Node(iElem,6) Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
			for i=1:3
				if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	
					c_plot_count = c_plot_count + 1;
					Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
					to_be_plot_count = to_be_plot_count +1;
					to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+1),1)];
					to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+1),2)];
					to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+1),3)];		
				end			
			end
			for i=5:7
				%检查该边线是否已经绘制
				if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	  
					c_plot_count = c_plot_count + 1;
					Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
					to_be_plot_count = to_be_plot_count +1;
					to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+1),1)];
					to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+1),2)];
					to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+1),3)];		
				end					  
			end
			for i=1:4
				%检查该边线是否已经绘制
				if not(ismember(sort([NN(i),NN(i+4)]),Ploted_Ele_lines,'rows')) 		
					c_plot_count = c_plot_count + 1;
					Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+4)]);
					to_be_plot_count = to_be_plot_count +1;
					to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+4),1)];
					to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+4),2)];
					to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+4),3)];						
				end				  
			end	
			%检查该边线是否已经绘制
			if not(ismember(sort([NN(1),NN(4)]),Ploted_Ele_lines,'rows')) 		
				c_plot_count = c_plot_count + 1;
				Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(1),NN(4)]);		
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),1) New_Node_Coor(NN(4),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),2) New_Node_Coor(NN(4),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),3) New_Node_Coor(NN(4),3)];				
			end
			if not(ismember(sort([NN(5),NN(8)]),Ploted_Ele_lines,'rows')) 	
				c_plot_count = c_plot_count + 1;
				Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(5),NN(8)]);	
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),1) New_Node_Coor(NN(8),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),2) New_Node_Coor(NN(8),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),3) New_Node_Coor(NN(8),3)];					
			end
		end 	
		plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',Line_width,'Color','green')	%灰色				
	end
	
    if Key_PLOT(2,9)==2  | Key_PLOT(2,9)==4 %增强单元的网格
		c_plot_count = 0;
		Ploted_Ele_lines =[0 0];	
		to_be_plot_count = 0;
		to_be_plot_x = [];to_be_plot_y = [];to_be_plot_z = [];  %2021-08-02		
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4) ...
				  Elem_Node(iElem,5) Elem_Node(iElem,6) ...
				  Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
  			%如果该单元含有增强节点	
			if isempty(Post_Enriched_Nodes)==0
				if sum(sum(Post_Enriched_Nodes(NN,1:size(Post_Enriched_Nodes,2))))~=0
					for i=1:3
						%检查该边线是否已经绘制(防止重复绘制)
						if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 		
							c_plot_count = c_plot_count + 1;
							Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
							c_plot_count = c_plot_count + 1;
							Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
							to_be_plot_count = to_be_plot_count +1;
							to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+1),1)];
							to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+1),2)];
							to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+1),3)];								
						end			
					end
					for i=5:7
						%检查该边线是否已经绘制
						if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 		
							c_plot_count = c_plot_count + 1;
							Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
							to_be_plot_count = to_be_plot_count +1;
							to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+1),1)];
							to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+1),2)];
							to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+1),3)];									
						end					  
					end
					for i=1:4
						%检查该边线是否已经绘制
						if not(ismember(sort([NN(i),NN(i+4)]),Ploted_Ele_lines,'rows')) 		
							c_plot_count = c_plot_count + 1;
							Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+4)]);
							to_be_plot_count = to_be_plot_count +1;
							to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),1) New_Node_Coor(NN(i+4),1)];
							to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),2) New_Node_Coor(NN(i+4),2)];
							to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(i),3) New_Node_Coor(NN(i+4),3)];	
						end				  
					end	
					%检查该边线是否已经绘制
					if not(ismember(sort([NN(1),NN(4)]),Ploted_Ele_lines,'rows')) 		
						c_plot_count = c_plot_count + 1;
						Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(1),NN(4)]);		
						to_be_plot_count = to_be_plot_count +1;
						to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),1) New_Node_Coor(NN(4),1)];
						to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),2) New_Node_Coor(NN(4),2)];
						to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(1),3) New_Node_Coor(NN(4),3)];							
					end
					if not(ismember(sort([NN(5),NN(8)]),Ploted_Ele_lines,'rows')) 	
						c_plot_count = c_plot_count + 1;
						Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(5),NN(8)]);			
						to_be_plot_count = to_be_plot_count +1;
						to_be_plot_x(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),1) New_Node_Coor(NN(8),1)];
						to_be_plot_y(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),2) New_Node_Coor(NN(8),2)];
						to_be_plot_z(to_be_plot_count,1:2) = [New_Node_Coor(NN(5),3) New_Node_Coor(NN(8),3)];							
					end					  
				end		
			end
		end 
        plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',Line_width,'Color',[.5 .5 .5])	%灰色			

		
		for i=1:size(Model_Outline,1)
			plot3([New_Node_Coor(Model_Outline(i,1),1),New_Node_Coor(Model_Outline(i,2),1)],...
				  [New_Node_Coor(Model_Outline(i,1),2),New_Node_Coor(Model_Outline(i,2),2)],...
				  [New_Node_Coor(Model_Outline(i,1),3),New_Node_Coor(Model_Outline(i,2),3)],'LineWidth',Line_width,'Color','black')	    
		end			
	end	
    if Key_PLOT(2,9)==3  | Key_PLOT(2,9)==4  %模型的表面网格	
		for i=1:size(Model_OutArea,1)
			plot3([New_Node_Coor(Model_OutArea(i,1),1),New_Node_Coor(Model_OutArea(i,2),1)],...
				  [New_Node_Coor(Model_OutArea(i,1),2),New_Node_Coor(Model_OutArea(i,2),2)],...
				  [New_Node_Coor(Model_OutArea(i,1),3),New_Node_Coor(Model_OutArea(i,2),3)],'LineWidth',Line_width,'Color',[.6 .5 .4])	 			  
		end	
		for i=1:size(Model_Outline,1)
			plot3([New_Node_Coor(Model_Outline(i,1),1),New_Node_Coor(Model_Outline(i,2),1)],...
				  [New_Node_Coor(Model_Outline(i,1),2),New_Node_Coor(Model_Outline(i,2),2)],...
				  [New_Node_Coor(Model_Outline(i,1),3),New_Node_Coor(Model_Outline(i,2),3)],'LineWidth',Line_width,'Color','black')	    
		end		
	end		
    if Key_PLOT(2,9) ==0 %仅绘制模型边框(Outlines)
		for i=1:size(Model_Outline,1)
			plot3([New_Node_Coor(Model_Outline(i,1),1),New_Node_Coor(Model_Outline(i,2),1)],...
				  [New_Node_Coor(Model_Outline(i,1),2),New_Node_Coor(Model_Outline(i,2),2)],...
				  [New_Node_Coor(Model_Outline(i,1),3),New_Node_Coor(Model_Outline(i,2),3)],'LineWidth',Line_width,'Color','black')	    
		end
	end
    if Key_PLOT(2,8) ==1 %变形前的网格边界
		for i=1:size(Model_Outline,1)
			plot3([Node_Coor(Model_Outline(i,1),1),Node_Coor(Model_Outline(i,2),1)],...
				  [Node_Coor(Model_Outline(i,1),2),Node_Coor(Model_Outline(i,2),2)],...
				  [Node_Coor(Model_Outline(i,1),3),Node_Coor(Model_Outline(i,2),3)],'--','LineWidth',Line_width,'Color','black')	    
		end
	end	
end

%绘制变形前的网格
Line_width =0.1;
if Key_PLOT(2,8)==1
    %不允许仅绘制变形前的网格而不绘制变形后的网格(2021-08-02)
	if Key_PLOT(2,9)~=1
		disp(' >> Error :: Key_PLOT(2,8)==1 but Key_PLOT(2,9)~=1, not allow!') 
		Error_Message
	end
	
	c_plot_count = 0;
	Ploted_Ele_lines =[0 0];
	to_be_plot_count = 0;
	to_be_plot_x = [];to_be_plot_y = [];to_be_plot_z = [];  %2021-08-02
	to_be_plot_x_1 = [];to_be_plot_y_1 = [];to_be_plot_z_1 = [];  %2021-08-02
	to_be_plot_x_2 = [];to_be_plot_y_2 = [];to_be_plot_z_2 = [];
	p_1 = [];
	p_2 = [];
	lines =[];
	for iElem = 1:Num_Elem
		NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) Elem_Node(iElem,3) Elem_Node(iElem,4) ...
			  Elem_Node(iElem,5) Elem_Node(iElem,6) Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
		for i=1:3
			if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	
				c_plot_count = c_plot_count + 1;
				Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+1),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+1),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+1),3)];		
			end			
		end
		for i=5:7
			%检查该边线是否已经绘制
			if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	  
				c_plot_count = c_plot_count + 1;
				Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+1),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+1),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+1),3)];		
				
			end					  
		end
		for i=1:4
			%检查该边线是否已经绘制
			if not(ismember(sort([NN(i),NN(i+4)]),Ploted_Ele_lines,'rows')) 		
				c_plot_count = c_plot_count + 1;
				Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+4)]);
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+4),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+4),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+4),3)];						
			end				  
		end	
		%检查该边线是否已经绘制
		if not(ismember(sort([NN(1),NN(4)]),Ploted_Ele_lines,'rows')) 		
			c_plot_count = c_plot_count + 1;
			Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(1),NN(4)]);		
			to_be_plot_count = to_be_plot_count +1;
			to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(1),1) Node_Coor(NN(4),1)];
			to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(1),2) Node_Coor(NN(4),2)];
			to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(1),3) Node_Coor(NN(4),3)];				
		end
		if not(ismember(sort([NN(5),NN(8)]),Ploted_Ele_lines,'rows')) 	
			c_plot_count = c_plot_count + 1;
			Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(5),NN(8)]);	
			to_be_plot_count = to_be_plot_count +1;
			to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(5),1) Node_Coor(NN(8),1)];
			to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(5),2) Node_Coor(NN(8),2)];
			to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(5),3) Node_Coor(NN(8),3)];					
		end
	end 	
	plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',Line_width,'Color',[0.5,0.5,0.5])	%灰色			
end
if Key_PLOT(2,8)==2	
	%仅绘制模型边框(Outlines)
	for i=1:size(Model_Outline,1)
			plot3([Node_Coor(Model_Outline(i,1),1),Node_Coor(Model_Outline(i,2),1)],...
				  [Node_Coor(Model_Outline(i,1),2),Node_Coor(Model_Outline(i,2),2)],...
				  [Node_Coor(Model_Outline(i,1),3),Node_Coor(Model_Outline(i,2),3)],'--','LineWidth',Line_width,'Color',[.5 .5 .5])	%灰色	    
	end
end

%绘制坐标轴
Arrow_length = (c_X_Length+c_Y_Length+c_Z_Length)/15;
h = Tools_mArrow3([0 0 0],[Arrow_length 0 0],'color','red','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 Arrow_length 0],'color','green','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 0 Arrow_length],'color','blue','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');



%绘制载荷
if Key_PLOT(2,7) == 1 || Key_PLOT(2,7) == 3
    disp(['      ----- Plotting forces of nodes......'])
    Max_x_Force = max(abs(FORCE_Matrix(:,1)));
	Max_y_Force = max(abs(FORCE_Matrix(:,2)));
	Max_z_Force = max(abs(FORCE_Matrix(:,3)));
	Max_Force   = max(Max_x_Force,max(Max_y_Force,Max_z_Force));
	
	% length of force arrow
    % REMOVE:length_arrow = sqrt(max_area_ele);
	length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/30.0;          
	
	% Loop through each node.
	for i = 1:Num_Node
	    if FORCE_Matrix(i,1) ~=0  | FORCE_Matrix(i,2) ~=0 | FORCE_Matrix(i,3) ~=0         % If the nodes has force load, then:
			c_force_x   = FORCE_Matrix(i,1);
			c_force_y   = FORCE_Matrix(i,2);
            c_force_z   = FORCE_Matrix(i,3);
			delta_L_x = c_force_x*length_arrow/Max_Force;
			delta_L_y = c_force_y*length_arrow/Max_Force;
			delta_L_z = c_force_z*length_arrow/Max_Force;
			
			StartPoint = [New_Node_Coor(i,1)-delta_L_x   New_Node_Coor(i,2)-delta_L_y     New_Node_Coor(i,3)-delta_L_z];
			EndPoint   = [New_Node_Coor(i,1)             New_Node_Coor(i,2)               New_Node_Coor(i,3)          ];
			%plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',2.0,'Color',[153/255,51/255,250/255]);
			%plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',1.5,'Color','red');
			h = Tools_mArrow3(StartPoint,EndPoint,'color','red','stemWidth',length_arrow/35.0,'tipWidth',length_arrow/15.0,'facealpha',1.0);	
		end	
	end
end

%绘制边界条件
if Key_PLOT(2,7) == 2 || Key_PLOT(2,7) == 3
    disp(['      ----- Plotting boundary conditions......'])
	length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/80.0;          
	for i = 1:size(Bou_x,1)
		delta_L_x = length_arrow/2;
		StartPoint = [New_Node_Coor(Bou_x(i),1)-delta_L_x   New_Node_Coor(Bou_x(i),2)     New_Node_Coor(Bou_x(i),3)];
		EndPoint   = [New_Node_Coor(Bou_x(i),1)+0   New_Node_Coor(Bou_x(i),2)     New_Node_Coor(Bou_x(i),3)];
		plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',1.0,'Color','black');
	end
	for i = 1:size(Bou_y,1)
		delta_L_y = length_arrow/2;
		StartPoint = [New_Node_Coor(Bou_y(i),1)   New_Node_Coor(Bou_y(i),2)-delta_L_y     New_Node_Coor(Bou_y(i),3)];
		EndPoint   = [New_Node_Coor(Bou_y(i),1)   New_Node_Coor(Bou_y(i),2)+0     New_Node_Coor(Bou_y(i),3)];
		plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',1.0,'Color','black');
	end
	for i = 1:size(Bou_z,1)
		delta_L_z = length_arrow/2;
		StartPoint = [New_Node_Coor(Bou_z(i),1)   New_Node_Coor(Bou_z(i),2)     New_Node_Coor(Bou_z(i),3)-delta_L_z];
		EndPoint   = [New_Node_Coor(Bou_z(i),1)   New_Node_Coor(Bou_z(i),2)     New_Node_Coor(Bou_z(i),3)+0];
		plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',1.0,'Color','black');
	end	
end
% Plot nodes
if isempty(Post_Enriched_Nodes) ~= 1
    length_min = min(c_X_Length,min(c_Y_Length,c_Z_Length))/5.0; 
	r = length_min/20;
    if Key_PLOT(2,2)==1
		disp(['      ----- Plotting nodes......'])
		for i =1 :size(Post_Enriched_Nodes,2)
			x_node =[];y_node =[];z_node =[];
			count = 0;		
			for j =1:Num_Node
			    count = count +1;
				x_node(count) = New_Node_Coor(j,1);                                          
				y_node(count) = New_Node_Coor(j,2);  
                z_node(count) = New_Node_Coor(j,3);  				
			end
		end
		plot3(x_node,y_node,z_node,'k.','MarkerSize',9.0)    % MarkerSize 表示点的大小,黑色点				
	end
end


% Plot enriched nodes
if isempty(Post_Enriched_Nodes) ~= 1
    length_min = min(c_X_Length,min(c_Y_Length,c_Z_Length))/5.0; 
	r = length_min/10;
    if Key_PLOT(2,14)==1
		disp(['      ----- Plotting enriched nodes......'])
		x_node_tip =[];y_node_tip =[];z_node_tip =[];
		count_tip =0;
		x_node_H=[];y_node_H =[];z_node_H =[];
		count_H =0;		
		for i =1 :size(Post_Enriched_Nodes,2)
			for j =1:Num_Node
				x_node = New_Node_Coor(j,1);                                          
				y_node = New_Node_Coor(j,2);  
                z_node = New_Node_Coor(j,3);  				
				if Post_Enriched_Nodes(j,i)==1     % Tip nodes
					count_tip = count_tip +1;
					x_node_tip(count_tip) = x_node;
					y_node_tip(count_tip) = y_node;
					z_node_tip(count_tip) = z_node;
				elseif Post_Enriched_Nodes(j,i)==2 % Heaviside nodes
					count_H = count_H +1;
					x_node_H(count_H) = x_node;
					y_node_H(count_H) = y_node;
					z_node_H(count_H) = z_node;
				elseif Post_Enriched_Nodes(j,i)==3 % Junction nodes
				    % TO BE DONE
				end
			end
		end
		plot3(x_node_tip,y_node_tip,z_node_tip,'r.','MarkerSize',13.0)    % MarkerSize 表示点的大小，b.表示绿色的点
		plot3(x_node_H,y_node_H,z_node_H,'b.','MarkerSize',9.0)    % MarkerSize 表示点的大小，b.表示绿色的点		
	end
end


%绘制增强节点的位移向量(节点的实际位移，不是增强自由度位移),2020-01-05
if Key_PLOT(2,12)==3
	% length of displacement arrow
	length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/25.0;  
	Max_Value    = max(max(abs(DISP(1:Num_Node,2:4))));
    if isempty(Post_Enriched_Nodes) ~= 1
		disp(['      ----- Plotting displacement vector of enriched nodes......'])
		for i =1 :size(Post_Enriched_Nodes,2)
            for j =1:Num_Node
			    if Post_Enriched_Nodes(j,i)~=0     %增强节点
					x_node = New_Node_Coor(j,1);                                          
					y_node = New_Node_Coor(j,2);  
					z_node = New_Node_Coor(j,3);  				
					c_dx   = DISP(j,2);
					c_dy   = DISP(j,3);
					c_dz   = DISP(j,4);
					delta_L_x = c_dx*length_arrow/Max_Value;
					delta_L_y = c_dy*length_arrow/Max_Value;
					delta_L_z = c_dz*length_arrow/Max_Value;
				
					StartPoint = [x_node,    y_node,     z_node];
					EndPoint   = [x_node+delta_L_x,    y_node+delta_L_y,     z_node+delta_L_z];
					% plot3([StartPoint(1),EndPoint(1)],[StartPoint(2),EndPoint(2)], [StartPoint(3),EndPoint(3)],'LineWidth',1.0,'Color','black');
					h = Tools_mArrow3(StartPoint,EndPoint,'color','black','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
				end
			end		
        end		
	end
end
				
%绘制裂缝面(模型之外的离散裂缝面节点不变形)
if Key_PLOT(2,5) == 1
    disp(['      ----- Plotting crack surface...'])
	if isempty(Crack_X)==0
		for i = 1:num_Crack(isub)
			% nPt = size(Crack_X{i},2);
			%--------------
			% option 1
			%--------------
			%目前一个裂缝只能由4个点构成
			% c_x = [Crack_X{i}(1:4)];
			% c_y = [Crack_Y{i}(1:4)];
			% c_z = [Crack_Z{i}(1:4)];
			% fill3(c_x,c_y,c_z,'r','FaceAlpha',0.8,'FaceLighting','gouraud')
			%--------------
			% option 2
			%--------------
			%绘制裂缝面离散点
			nnode = size(Crack_node_X{i},2);
			for j=1:nnode
				c_node_x = [Crack_node_X{i}(j)];
				c_node_y = [Crack_node_Y{i}(j)];
				c_node_z = [Crack_node_Z{i}(j)];
				% 获得离散裂缝面节点所在模型单元号
				c_Ele = Crack_Node_in_Ele{i}(j);
				if c_Ele==0   %离散裂缝的节点在模型外
				    plot3(c_node_x,c_node_y,c_node_z,'c.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点
					% ts = text(c_node_x,c_node_y,c_node_z,num2str(j),'Color','blue','FontSize',12,'FontName','Consolas','FontAngle','italic');
					Crack_node_X_new{i}(j) = c_node_x;
					Crack_node_Y_new{i}(j) = c_node_y;
					Crack_node_Z_new{i}(j) = c_node_z;
				else          %离散裂缝的节点在模型内
				    % Get the local coordinates of the points of the crack. 
					Kesi = Crack_node_local_X{i}(j); 
					Yita = Crack_node_local_Y{i}(j); 
					Zeta = Crack_node_local_Z{i}(j); 
					NN = [Elem_Node(c_Ele,1) Elem_Node(c_Ele,2) ...
						  Elem_Node(c_Ele,3) Elem_Node(c_Ele,4) ...
						  Elem_Node(c_Ele,5) Elem_Node(c_Ele,6) ...
						  Elem_Node(c_Ele,7) Elem_Node(c_Ele,8)]; 		
					% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
				    [N]  = Cal_N_3D(Kesi,Yita,Zeta);
				    dis_x(j) =     DISP(NN(1),2)*N(1) + DISP(NN(2),2)*N(2) + DISP(NN(3),2)*N(3) + DISP(NN(4),2)*N(4) ...
                                 + DISP(NN(5),2)*N(5) + DISP(NN(6),2)*N(6) + DISP(NN(7),2)*N(7) + DISP(NN(8),2)*N(8); 				
				    dis_y(j) =     DISP(NN(1),3)*N(1) + DISP(NN(2),3)*N(2) + DISP(NN(3),3)*N(3) + DISP(NN(4),3)*N(4) ...
                                 + DISP(NN(5),3)*N(5) + DISP(NN(6),3)*N(6) + DISP(NN(7),3)*N(7) + DISP(NN(8),3)*N(8);	
				    dis_z(j) =     DISP(NN(1),4)*N(1) + DISP(NN(2),4)*N(2) + DISP(NN(3),4)*N(3) + DISP(NN(4),4)*N(4) ...
                                 + DISP(NN(5),4)*N(5) + DISP(NN(6),4)*N(6) + DISP(NN(7),4)*N(7) + DISP(NN(8),4)*N(8); 
					if isnan(dis_x(j)) ==1
					    dis_x(j) = 0.0;
						% disp(['            WARNING :: dis_x of node ',num2str(j),' is NAN!']) 
					end
					if isnan(dis_y(j)) ==1
					    dis_y(j) = 0.0;
						% disp(['            WARNING :: dis_y of node ',num2str(j),' is NAN!']) 
					end
					if isnan(dis_z(j)) ==1
					    dis_z(j) = 0.0;	
                        % disp(['            WARNING :: dis_z of node ',num2str(j),' is NAN!']) 						
                    end						
                    last_c_node_x = c_node_x +	dis_x(j)*scale;	
					last_c_node_y = c_node_y +	dis_y(j)*scale;
					last_c_node_z = c_node_z +	dis_z(j)*scale;	
					Crack_node_X_new{i}(j) = last_c_node_x;
					Crack_node_Y_new{i}(j) = last_c_node_y;
					Crack_node_Z_new{i}(j) = last_c_node_z;					
                    %绘制离散裂缝点		
                    if Key_PLOT(2,5) == 2					
        			    plot3(last_c_node_x,last_c_node_y,last_c_node_z,'c.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点 
					end
					%显示离散裂缝点编号(测试用)
					% ts = text(c_node_x,c_node_y,c_node_z,num2str(j),'Color','blue','FontSize',17,'FontName','Consolas','FontAngle','italic');
					
					%绘制离散裂缝面边界节点的主应力方向(测试用)
					if Key_PLOT(2,11) == 1
						length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/60.0;  
						StartPoint = [last_c_node_x last_c_node_y last_c_node_z];
						EndPoint   = [last_c_node_x+length_arrow*Crack3D_Vector_S1_X{i}(j) last_c_node_y+length_arrow*Crack3D_Vector_S1_Y{i}(j) last_c_node_z+length_arrow*Crack3D_Vector_S1_Z{i}(j)];
						h = Tools_mArrow3(StartPoint,EndPoint,'color','m','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
						EndPoint   = [last_c_node_x-length_arrow*Crack3D_Vector_S1_X{i}(j) last_c_node_y-length_arrow*Crack3D_Vector_S1_Y{i}(j) last_c_node_z-length_arrow*Crack3D_Vector_S1_Z{i}(j)];	
						h = Tools_mArrow3(StartPoint,EndPoint,'color','m','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);		
                    end					
				end
			end
			%绘制单元面
			nele = size(Crack_Ele_1{i},2);
			c_x =[];c_y =[];c_z =[];
			for j=1:nele
				c_x(j,1:3) = [Crack_node_X_new{i}(Crack_Ele_1{i}(j)),Crack_node_X_new{i}(Crack_Ele_2{i}(j)),Crack_node_X_new{i}(Crack_Ele_3{i}(j))];
				c_y(j,1:3) = [Crack_node_Y_new{i}(Crack_Ele_1{i}(j)),Crack_node_Y_new{i}(Crack_Ele_2{i}(j)),Crack_node_Y_new{i}(Crack_Ele_3{i}(j))];
				c_z(j,1:3) = [Crack_node_Z_new{i}(Crack_Ele_1{i}(j)),Crack_node_Z_new{i}(Crack_Ele_2{i}(j)),Crack_node_Z_new{i}(Crack_Ele_3{i}(j))];					
			end		
            patch(c_x',c_y',c_z',[1,1,0],'FaceAlpha',0.5,'FaceLighting','gouraud')				
		end	
	end
end

%绘制裂缝体,2021-02-11.
if Key_PLOT(2,5) == 2
    disp(['      ----- Plotting crack volume...'])
	if isempty(Crack_X)==0
		for i_C = 1:num_Crack(isub)
			if i_C==1
				c_clor = 'green';
			elseif i_C==2
				c_clor = 'red';
			elseif i_C==3
				c_clor = 'blue';
			elseif i_C==4	
				c_clor = 'cyan';
			elseif i_C==5	
				c_clor = 'magenta';
			elseif i_C==6	
				c_clor = 'yellow';
			end	
			if isempty(Crack3D_FluEl_Aperture)==0
				nfluid_Ele = Crack3D_Fluid_Ele_Num{i_C};
				for j=1:nfluid_Ele
					c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i_C,j,1:3);
					old_coor_x =Crack3D_CalP_X{i_C}(c_fluid_nodes);
					old_coor_y =Crack3D_CalP_Y{i_C}(c_fluid_nodes);
					old_coor_z =Crack3D_CalP_Z{i_C}(c_fluid_nodes);
					Up_dis_V_x  = Cracks_CalP_UpDis_3D_X{i_C}(c_fluid_nodes);
					Up_dis_V_y  = Cracks_CalP_UpDis_3D_Y{i_C}(c_fluid_nodes);
					Up_dis_V_z  = Cracks_CalP_UpDis_3D_Z{i_C}(c_fluid_nodes);
					new_coor_x = old_coor_x + scale*Up_dis_V_x;
					new_coor_y = old_coor_y + scale*Up_dis_V_y;
					new_coor_z = old_coor_z + scale*Up_dis_V_z;
					%填充上表面
                    fill3(new_coor_x,new_coor_y,new_coor_z,c_clor,'FaceAlpha',0.5,'FaceLighting','gouraud')				
					Low_dis_V_x  = Cracks_CalP_LowDis_3D_X{i_C}(c_fluid_nodes);
					Low_dis_V_y  = Cracks_CalP_LowDis_3D_Y{i_C}(c_fluid_nodes);
					Low_dis_V_z  = Cracks_CalP_LowDis_3D_Z{i_C}(c_fluid_nodes);
					new_coor_x = old_coor_x + scale*Low_dis_V_x;
					new_coor_y = old_coor_y + scale*Low_dis_V_y;
					new_coor_z = old_coor_z + scale*Low_dis_V_z;
					%填充下表面
                    fill3(new_coor_x,new_coor_y,new_coor_z,c_clor,'FaceAlpha',0.5,'FaceLighting','gouraud')						
				end
			end
		end	
	end
end


%绘制离散裂缝面边界节点的局部坐标系
if Key_PLOT(2,10) == 1
    for i = 1:num_Crack(isub)
	    for j=1:size(Crack3D_Meshed_Outline{i},2)
		    Vector_x_x = Crack3D_Vertex_Vector_X_X{i}(j);
			Vector_x_y = Crack3D_Vertex_Vector_X_Y{i}(j);
			Vector_x_z = Crack3D_Vertex_Vector_X_Z{i}(j);
		    Vector_y_x = Crack3D_Vertex_Vector_Y_X{i}(j);
			Vector_y_y = Crack3D_Vertex_Vector_Y_Y{i}(j);
			Vector_y_z = Crack3D_Vertex_Vector_Y_Z{i}(j);
		    Vector_z_x = Crack3D_Vertex_Vector_Z_X{i}(j);
			Vector_z_y = Crack3D_Vertex_Vector_Z_Y{i}(j);
			Vector_z_z = Crack3D_Vertex_Vector_Z_Z{i}(j);	
            Crack_Node = Crack3D_Meshed_Outline{i}(j);	
			c_node_x = [Crack_node_X{i}(Crack_Node)];
			c_node_y = [Crack_node_Y{i}(Crack_Node)];
			c_node_z = [Crack_node_Z{i}(Crack_Node)];			
			% 获得离散裂缝面节点所在模型单元号
			c_Ele = Crack_Node_in_Ele{i}(Crack_Node);
			% 单元号不能为0
			if  c_Ele ~= 0
				% Get the local coordinates of the points of the crack. 
				Kesi = Crack_node_local_X{i}(Crack_Node); 
				Yita = Crack_node_local_Y{i}(Crack_Node); 
				Zeta = Crack_node_local_Z{i}(Crack_Node); 
				NN = [Elem_Node(c_Ele,1) Elem_Node(c_Ele,2) ...
					  Elem_Node(c_Ele,3) Elem_Node(c_Ele,4) ...
					  Elem_Node(c_Ele,5) Elem_Node(c_Ele,6) ...
					  Elem_Node(c_Ele,7) Elem_Node(c_Ele,8)]; 		
				% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
				[N]  = Cal_N_3D(Kesi,Yita,Zeta);
				c_dis_x =     DISP(NN(1),2)*N(1) + DISP(NN(2),2)*N(2) + DISP(NN(3),2)*N(3) + DISP(NN(4),2)*N(4) ...
							 + DISP(NN(5),2)*N(5) + DISP(NN(6),2)*N(6) + DISP(NN(7),2)*N(7) + DISP(NN(8),2)*N(8); 				
				c_dis_y =     DISP(NN(1),3)*N(1) + DISP(NN(2),3)*N(2) + DISP(NN(3),3)*N(3) + DISP(NN(4),3)*N(4) ...
							 + DISP(NN(5),3)*N(5) + DISP(NN(6),3)*N(6) + DISP(NN(7),3)*N(7) + DISP(NN(8),3)*N(8);	
				c_dis_z =     DISP(NN(1),4)*N(1) + DISP(NN(2),4)*N(2) + DISP(NN(3),4)*N(3) + DISP(NN(4),4)*N(4) ...
							 + DISP(NN(5),4)*N(5) + DISP(NN(6),4)*N(6) + DISP(NN(7),4)*N(7) + DISP(NN(8),4)*N(8); 
				last_c_node_x = c_node_x +	c_dis_x*scale;	
				last_c_node_y = c_node_y +	c_dis_y*scale;
				last_c_node_z = c_node_z +	c_dis_z*scale;	
				% Crack_node_X_new{i}(j) = last_c_node_x;
				% Crack_node_Y_new{i}(j) = last_c_node_y;
				% Crack_node_Z_new{i}(j) = last_c_node_z;	
				
				%显示离散裂缝面边界节点编号(测试用)
				% ts = text(last_c_node_x,last_c_node_y,last_c_node_z,num2str(j),'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');
					
				% 以下以彩色绘制
				length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/80.0;  
				StartPoint = [last_c_node_x last_c_node_y last_c_node_z];
				EndPoint   = [last_c_node_x+length_arrow*Vector_x_x last_c_node_y+length_arrow*Vector_x_y last_c_node_z+length_arrow*Vector_x_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','r','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
				ts = text(last_c_node_x+1.1*length_arrow*Vector_x_x,last_c_node_y+1.1*length_arrow*Vector_x_y,last_c_node_z+1.1*length_arrow*Vector_x_z,  ...
																"x",'Color','red','FontSize',12,'FontName','Consolas','FontAngle','italic');											
				EndPoint   = [last_c_node_x+length_arrow*Vector_y_x last_c_node_y+length_arrow*Vector_y_y last_c_node_z+length_arrow*Vector_y_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','g','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
				ts = text(last_c_node_x+1.1*length_arrow*Vector_y_x,last_c_node_y+1.1*length_arrow*Vector_y_y,last_c_node_z+1.1*length_arrow*Vector_y_z,  ...
																"y",'Color','g','FontSize',12,'FontName','Consolas','FontAngle','italic');						
				EndPoint   = [last_c_node_x+length_arrow*Vector_z_x last_c_node_y+length_arrow*Vector_z_y last_c_node_z+length_arrow*Vector_z_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','b','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);		
				ts = text(last_c_node_x+1.1*length_arrow*Vector_z_x,last_c_node_y+1.1*length_arrow*Vector_z_y,last_c_node_z+1.1*length_arrow*Vector_z_z,  ...
																"z",'Color','b','FontSize',12,'FontName','Consolas','FontAngle','italic');		
				% 以下以黑绘制
				% length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/80.0;  
				% StartPoint = [last_c_node_x last_c_node_y last_c_node_z];
				% EndPoint   = [last_c_node_x+length_arrow*Vector_x_x last_c_node_y+length_arrow*Vector_x_y last_c_node_z+length_arrow*Vector_x_z];
				% h = Tools_mArrow3(StartPoint,EndPoint,'color','black','stemWidth',0.020,'facealpha',1.0);
				% ts = text(last_c_node_x+1.1*length_arrow*Vector_x_x,last_c_node_y+1.1*length_arrow*Vector_x_y,last_c_node_z+1.1*length_arrow*Vector_x_z,  ...
																% "x",'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');											
				% EndPoint   = [last_c_node_x+length_arrow*Vector_y_x last_c_node_y+length_arrow*Vector_y_y last_c_node_z+length_arrow*Vector_y_z];
				% h = Tools_mArrow3(StartPoint,EndPoint,'color','black','stemWidth',0.020,'facealpha',1.0);
				% ts = text(last_c_node_x+1.1*length_arrow*Vector_y_x,last_c_node_y+1.1*length_arrow*Vector_y_y,last_c_node_z+1.1*length_arrow*Vector_y_z,  ...
																% "y",'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');					
				% EndPoint   = [last_c_node_x+length_arrow*Vector_z_x last_c_node_y+length_arrow*Vector_z_y last_c_node_z+length_arrow*Vector_z_z];
				% h = Tools_mArrow3(StartPoint,EndPoint,'color','black','stemWidth',0.020,'facealpha',1.0);		
				% ts = text(last_c_node_x+1.1*length_arrow*Vector_z_x,last_c_node_y+1.1*length_arrow*Vector_z_y,last_c_node_z+1.1*length_arrow*Vector_z_z,  ...
																% "z",'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');	
				% EndPoint   = [last_c_node_x-length_arrow*Crack3D_Vector_S1_X{i}(j) last_c_node_y-length_arrow*Crack3D_Vector_S1_Y{i}(j) last_c_node_z-length_arrow*Crack3D_Vector_S1_Z{i}(j)];	
				% h = Tools_mArrow3(StartPoint,EndPoint,'color','m','stemWidth',0.025,'facealpha',1.0);		
			end
		end
	end
end	

%裂尖增强单元的基准线(baseline)和基准线的局部坐标系,2020-01-05
if Key_PLOT(2,10) == 2 |  Key_PLOT(2,10) == 3
	for iElem = 1:Num_Elem
	    if sum(abs(Tip_Enriched_Ele_BaseLine(iElem,1:6)))>1.0D-10
			Point_A = Tip_Enriched_Ele_BaseLine(iElem,1:3);
			Point_B = Tip_Enriched_Ele_BaseLine(iElem,4:6);
			plot3([Point_A(1),Point_B(1)],[Point_A(2),Point_B(2)],[Point_A(3),Point_B(3)],'LineWidth',1.5,'Color','m')	
			plot3(Point_A(1),Point_A(2),Point_A(3),'g.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点 
			plot3(Point_B(1),Point_B(2),Point_B(3),'g.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点 
			%基准线AB的中心点
			mid_point_x =  (Point_A(1)+Point_B(1))/2.0;
			mid_point_y =  (Point_A(2)+Point_B(2))/2.0;
			mid_point_z =  (Point_A(3)+Point_B(3))/2.0;
			plot3(mid_point_x,mid_point_y,mid_point_z,'*','MarkerSize',15.0,'Color','m');    % MarkerSize 表示点的大小,青色点
			%绘制基准线的局部坐标系
			Vector_x_x = Tip_Enriched_Ele_BaseLine_Vector_x(iElem,1);
			Vector_x_y = Tip_Enriched_Ele_BaseLine_Vector_x(iElem,2);
			Vector_x_z = Tip_Enriched_Ele_BaseLine_Vector_x(iElem,3);
			Vector_y_x = Tip_Enriched_Ele_BaseLine_Vector_y(iElem,1);
			Vector_y_y = Tip_Enriched_Ele_BaseLine_Vector_y(iElem,2);
			Vector_y_z = Tip_Enriched_Ele_BaseLine_Vector_y(iElem,3);
			Vector_z_x = Tip_Enriched_Ele_BaseLine_Vector_z(iElem,1);
			Vector_z_y = Tip_Enriched_Ele_BaseLine_Vector_z(iElem,2);
			Vector_z_z = Tip_Enriched_Ele_BaseLine_Vector_z(iElem,3);		
			
			if Key_PLOT(2,10) == 3
				length_arrow = (c_X_Length+c_Y_Length+c_Z_Length)/60.0;  %局部坐标系的大小
				StartPoint = [mid_point_x,mid_point_y,mid_point_z];
				EndPoint   = [mid_point_x+length_arrow*Vector_x_x mid_point_y+length_arrow*Vector_x_y mid_point_z+length_arrow*Vector_x_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','r','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
				ts = text(mid_point_x+1.1*length_arrow*Vector_x_x,mid_point_y+1.1*length_arrow*Vector_x_y,mid_point_z+1.1*length_arrow*Vector_x_z,  ...
																"x",'Color','red','FontSize',12,'FontName','Consolas','FontAngle','italic');											
				EndPoint   = [mid_point_x+length_arrow*Vector_y_x mid_point_y+length_arrow*Vector_y_y mid_point_z+length_arrow*Vector_y_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','g','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);
				ts = text(mid_point_x+1.1*length_arrow*Vector_y_x,mid_point_y+1.1*length_arrow*Vector_y_y,mid_point_z+1.1*length_arrow*Vector_y_z,  ...
																"y",'Color','g','FontSize',12,'FontName','Consolas','FontAngle','italic');						
				EndPoint   = [mid_point_x+length_arrow*Vector_z_x mid_point_y+length_arrow*Vector_z_y mid_point_z+length_arrow*Vector_z_z];
				h = Tools_mArrow3(StartPoint,EndPoint,'color','b','stemWidth',length_arrow/25.0,'tipWidth',length_arrow/10.0,'facealpha',1.0);		
				ts = text(mid_point_x+1.1*length_arrow*Vector_z_x,mid_point_y+1.1*length_arrow*Vector_z_y,mid_point_z+1.1*length_arrow*Vector_z_z,  ...
																"z",'Color','b','FontSize',12,'FontName','Consolas','FontAngle','italic');		
			end	
        end		
	end
end	

%裂尖增强节点对应的单元(参考单元),2020-02-10
if Key_PLOT(2,10) == 4 
    if sum(abs(Tip_Enriched_Node_Ref_Element))>1.0D-10
	    for i_Node=1:Num_Node
		    c_Elem  = Tip_Enriched_Node_Ref_Element(i_Node,1);
			c_Crack = Tip_Enriched_Node_Ref_Element(i_Node,2);
			if c_Elem >=1
				NN = [Elem_Node(c_Elem,1) Elem_Node(c_Elem,2) ...
					  Elem_Node(c_Elem,3) Elem_Node(c_Elem,4) ...
					  Elem_Node(c_Elem,5) Elem_Node(c_Elem,6) ...
					  Elem_Node(c_Elem,7) Elem_Node(c_Elem,8)];                             % Nodes for current element
				%b(bule) c(cyan) g(green) k(key) m(magenta) r(red) w(white) y(yellow)
				if c_Crack==1
				    c_clor = 'red';
				elseif c_Crack==2
				    c_clor = 'black';
				elseif c_Crack==3
				    c_clor = 'green';
                elseif c_Crack==4	
				    c_clor = 'cyan';
                elseif c_Crack==5	
				    c_clor = 'magenta';
                elseif c_Crack==6	
				    c_clor = 'yellow';
                end				
				for i=1:3
					plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+1),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+1),2)],[Node_Coor(NN(i),3),Node_Coor(NN(i+1),3)],'LineWidth',Line_width,'Color',c_clor)	
				end
				for i=5:7
					plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+1),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+1),2)],[Node_Coor(NN(i),3),Node_Coor(NN(i+1),3)],'LineWidth',Line_width,'Color',c_clor)	
				end
				for i=1:4
					plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+4),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+4),2)],[Node_Coor(NN(i),3),Node_Coor(NN(i+4),3)],'LineWidth',Line_width,'Color',c_clor)	
				end	
				plot3([Node_Coor(NN(1),1),Node_Coor(NN(4),1)],[Node_Coor(NN(1),2),Node_Coor(NN(4),2)],[Node_Coor(NN(1),3),Node_Coor(NN(4),3)],'LineWidth',Line_width,'Color',c_clor)		
				plot3([Node_Coor(NN(5),1),Node_Coor(NN(8),1)],[Node_Coor(NN(5),2),Node_Coor(NN(8),2)],[Node_Coor(NN(5),3),Node_Coor(NN(8),3)],'LineWidth',Line_width,'Color',c_clor)
			end
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
