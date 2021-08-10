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

function Plot_Mesh3D(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
% This function plots the initial geometry,三维.

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
global Elem_Material Num_Step_to_Plot
global Crack_node_X Crack_node_Y Crack_node_Z
global Crack_Ele_1 Crack_Ele_2 Crack_Ele_3
global Model_Outline
global Crack_node_local_X Crack_node_local_Y Crack_node_local_Z
global Crack_Node_in_Ele
global Crack3D_CalP_X Crack3D_CalP_Y Crack3D_CalP_Z
global Crack3D_Vertex_Vector_X_X Crack3D_Vertex_Vector_X_Y Crack3D_Vertex_Vector_X_Z 
global Crack3D_Vertex_Vector_Y_X Crack3D_Vertex_Vector_Y_Y Crack3D_Vertex_Vector_Y_Z  
global Crack3D_Vertex_Vector_Z_X Crack3D_Vertex_Vector_Z_Y Crack3D_Vertex_Vector_Z_Z 	
global Crack3D_Meshed_Outline Crack3D_Fluid_Ele_Num
global Crack3D_Fluid_Ele_Nodes Crack3D_CalP_Aperture
global Cracks_FluNode_Vector_3D_X Cracks_FluNode_Vector_3D_Y Cracks_FluNode_Vector_3D_Z
global Max_Aperture_of_each_Crack Crack3D_FluEl_Aperture Min_Aperture_of_each_Crack
global Tip_Enriched_Ele_BaseLine Tip_Enriched_Ele_BaseLine_Vector_x 
global Tip_Enriched_Ele_BaseLine_Vector_y Tip_Enriched_Ele_BaseLine_Vector_z
global FluidEle_GaussNor_3D_X FluidEle_GaussNor_3D_Y FluidEle_GaussNor_3D_Z
global FluidEle_LCS_VectorX_X FluidEle_LCS_VectorX_Y FluidEle_LCS_VectorX_Z
global FluidEle_LCS_VectorY_X FluidEle_LCS_VectorY_Y FluidEle_LCS_VectorY_Z
global FluidEle_LCS_VectorZ_X FluidEle_LCS_VectorZ_Y FluidEle_LCS_VectorZ_Z
global FluidEle_Contact_Force_X FluidEle_Contact_Force_Y FluidEle_Contact_Force_Z
global Title_Font Key_Figure_Control_Widget
global G_X_NODES G_Y_NODES G_Z_NODES G_NN G_X_Min G_X_Max G_Y_Min G_Y_Max G_Z_Min G_Z_Max

% global Tri_BCD
disp(['      ----- Plotting undeformed mesh......'])
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

c_X_Length = Max_X_Coor-Min_X_Coor;
c_Y_Length = Max_Y_Coor-Min_Y_Coor;
c_Z_Length = Max_Z_Coor-Min_Z_Coor;

%===================  
%   New figure
%===================  
Tools_New_Figure
hold on;
title('Finite Element Mesh','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
delta = sqrt(aveg_area_ele);
axis([Min_X_Coor-delta Max_X_Coor+delta ...
      Min_Y_Coor-delta Max_Y_Coor+delta ...
	  Min_Z_Coor-delta Max_Z_Coor+delta]);
view(-37.5,30)  %视角
% view(0,90)        %xy视角
% view(180,0)        %xz视角
% view(90,0)        %yz视角
% view(-7,-15)  %视角

%===================  
%    绘制坐标轴
%===================  
Arrow_length = (c_X_Length+c_Y_Length+c_Z_Length)/15;
h = Tools_mArrow3([0 0 0],[Arrow_length 0 0],'color','red','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 Arrow_length 0],'color','green','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 0 Arrow_length],'color','blue','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');

%===================  	  
%   绘制单元网格
%===================  
Line_width =0.1;
if Key_PLOT(1,9) ==1  && Key_PLOT(1,3)~=1
    c_plot_count = 0;
	Ploted_Ele_lines =[0 0];
	to_be_plot_count = 0;
	to_be_plot_x = [];to_be_plot_y = [];to_be_plot_z = [];  %2021-08-02
    for iElem = 1:Num_Elem
	    % iElem
		NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) Elem_Node(iElem,3) Elem_Node(iElem,4) ...
			  Elem_Node(iElem,5) Elem_Node(iElem,6) Elem_Node(iElem,7) Elem_Node(iElem,8)];      % Nodes for current element
	    % if iElem ==92
		    % NN
		% end
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
		    %%%%%%检查该边线是否已经绘制
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
		    %%%%%%%检查该边线是否已经绘制
		    if not(ismember(sort([NN(i),NN(i+4)]),Ploted_Ele_lines,'rows')) 		
			    c_plot_count = c_plot_count + 1;
			    Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+4)]);
				to_be_plot_count = to_be_plot_count +1;
				to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+4),1)];
				to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+4),2)];
				to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+4),3)];						
            end				  
		end	
		%%%%%%%%%检查该边线是否已经绘制
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
	plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',Line_width,'Color',[.5 .5 .5])	%灰色		
	
elseif Key_PLOT(1,9) ==0 %仅绘制模型边框(Outlines)
    for i=1:size(Model_Outline,1)
			plot3([Node_Coor(Model_Outline(i,1),1),Node_Coor(Model_Outline(i,2),1)],...
				  [Node_Coor(Model_Outline(i,1),2),Node_Coor(Model_Outline(i,2),2)],...
				  [Node_Coor(Model_Outline(i,1),3),Node_Coor(Model_Outline(i,2),3)],'LineWidth',Line_width,'Color','black')    
	end	
end
if Key_PLOT(1,10) >=1%绘制给定的单元及其节点号
	Plot_Elem      = Key_PLOT(1,10);
	if (Key_PLOT(1,11) > Key_PLOT(1,10)) &  (Key_PLOT(1,11)<=Num_Elem )
	    Plot_Elem_End  = Key_PLOT(1,11);
	else
	    Plot_Elem_End  = Plot_Elem;
	end
	if Plot_Elem <=Num_Elem 
		for iElem = Plot_Elem:Plot_Elem_End
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4) ...
				  Elem_Node(iElem,5) Elem_Node(iElem,6) ...
				  Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
			%节点编号
			for i=1:8	  
				ts = text(Node_Coor(NN(i),1),Node_Coor(NN(i),2),Node_Coor(NN(i),3),num2str(NN(i)),'Color','blue','FontSize',10,'FontName','Consolas','FontAngle','italic');
			end
			for i=1:3
				plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+1),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+1),2)],...
					  [Node_Coor(NN(i),3),Node_Coor(NN(i+1),3)],'LineWidth',1.0,'Color','m')	
			end
			for i=5:7
				plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+1),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+1),2)],...
					  [Node_Coor(NN(i),3),Node_Coor(NN(i+1),3)],'LineWidth',1.0,'Color','m')	
			end
			for i=1:4
				plot3([Node_Coor(NN(i),1),Node_Coor(NN(i+4),1)],[Node_Coor(NN(i),2),Node_Coor(NN(i+4),2)],...
					  [Node_Coor(NN(i),3),Node_Coor(NN(i+4),3)],'LineWidth',1.0,'Color','m')	
			end	
			plot3([Node_Coor(NN(1),1),Node_Coor(NN(4),1)],[Node_Coor(NN(1),2),Node_Coor(NN(4),2)],...
				  [Node_Coor(NN(1),3),Node_Coor(NN(4),3)],'LineWidth',1.0,'Color','m')		
				  
			plot3([Node_Coor(NN(5),1),Node_Coor(NN(8),1)],[Node_Coor(NN(5),2),Node_Coor(NN(8),2)],...
				  [Node_Coor(NN(5),3),Node_Coor(NN(8),3)],'LineWidth',1.0,'Color','m')
		end
	end
end
%===================  
%   绘制单元面
%===================  
Color_3D_ele_face = [189/255,252/255,201/255];
FaceAlpha_3D_ele_face = 0.8;
if Key_PLOT(1,3) == 1
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
	    c_x_1(1:4,count_1) = [Node_Coor(NN(1),1),Node_Coor(NN(2),1),Node_Coor(NN(3),1),Node_Coor(NN(4),1)];
	    c_y_1(1:4,count_1) = [Node_Coor(NN(1),2),Node_Coor(NN(2),2),Node_Coor(NN(3),2),Node_Coor(NN(4),2)];
	    c_z_1(1:4,count_1) = [Node_Coor(NN(1),3),Node_Coor(NN(2),3),Node_Coor(NN(3),3),Node_Coor(NN(4),3)];		
	  end	  
	  %第2个面
	  if not(ismember(sort([NN(5),NN(6),NN(7),NN(8)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(5),NN(6),NN(7),NN(8)]);	
		count_2 = count_2 +1;
		c_x_2(1:4,count_2) = [Node_Coor(NN(5),1),Node_Coor(NN(6),1),Node_Coor(NN(7),1),Node_Coor(NN(8),1)];
		c_y_2(1:4,count_2) = [Node_Coor(NN(5),2),Node_Coor(NN(6),2),Node_Coor(NN(7),2),Node_Coor(NN(8),2)];
		c_z_2(1:4,count_2) = [Node_Coor(NN(5),3),Node_Coor(NN(6),3),Node_Coor(NN(7),3),Node_Coor(NN(8),3)];			
	  end	
	  %第3个面
	  if not(ismember(sort([NN(1),NN(2),NN(6),NN(5)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(1),NN(2),NN(6),NN(5)]);	
		count_3 = count_3 +1;
		c_x_3(1:4,count_3) = [Node_Coor(NN(1),1),Node_Coor(NN(2),1),Node_Coor(NN(6),1),Node_Coor(NN(5),1)];
		c_y_3(1:4,count_3) = [Node_Coor(NN(1),2),Node_Coor(NN(2),2),Node_Coor(NN(6),2),Node_Coor(NN(5),2)];
		c_z_3(1:4,count_3) = [Node_Coor(NN(1),3),Node_Coor(NN(2),3),Node_Coor(NN(6),3),Node_Coor(NN(5),3)];			
	  end	
	  %第4个面
	  if not(ismember(sort([NN(2),NN(3),NN(7),NN(6)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(2),NN(3),NN(7),NN(6)]);	
		count_4 = count_4 +1;
	    c_x_4(1:4,count_4) = [Node_Coor(NN(2),1),Node_Coor(NN(3),1),Node_Coor(NN(7),1),Node_Coor(NN(6),1)];
	    c_y_4(1:4,count_4) = [Node_Coor(NN(2),2),Node_Coor(NN(3),2),Node_Coor(NN(7),2),Node_Coor(NN(6),2)];
	    c_z_4(1:4,count_4) = [Node_Coor(NN(2),3),Node_Coor(NN(3),3),Node_Coor(NN(7),3),Node_Coor(NN(6),3)];			
	  end
	  %第5个面
	  if not(ismember(sort([NN(7),NN(8),NN(4),NN(3)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(7),NN(8),NN(4),NN(3)]);	
		count_5 = count_5 +1;
        c_x_5(1:4,count_5) = [Node_Coor(NN(7),1),Node_Coor(NN(8),1),Node_Coor(NN(4),1),Node_Coor(NN(3),1)];
        c_y_5(1:4,count_5) = [Node_Coor(NN(7),2),Node_Coor(NN(8),2),Node_Coor(NN(4),2),Node_Coor(NN(3),2)];
        c_z_5(1:4,count_5) = [Node_Coor(NN(7),3),Node_Coor(NN(8),3),Node_Coor(NN(4),3),Node_Coor(NN(3),3)];			
	  end
	  %第6个面
	  if not(ismember(sort([NN(5),NN(1),NN(4),NN(8)]),Ploted_Ele_Face,'rows')) %检查该面是否已经绘制(防止重复绘制)
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_Face(c_plot_count,1:4) = sort([NN(5),NN(1),NN(4),NN(8)]);	
		count_6 = count_6 +1;
	    c_x_6(1:4,count_6) = [Node_Coor(NN(5),1),Node_Coor(NN(1),1),Node_Coor(NN(4),1),Node_Coor(NN(8),1)];
	    c_y_6(1:4,count_6) = [Node_Coor(NN(5),2),Node_Coor(NN(1),2),Node_Coor(NN(4),2),Node_Coor(NN(8),2)];
	    c_z_6(1:4,count_6) = [Node_Coor(NN(5),3),Node_Coor(NN(1),3),Node_Coor(NN(4),3),Node_Coor(NN(8),3)];		
	  end
  end 
  patch(c_x_1(1:4,1:count_1),c_y_1(1:4,1:count_1),c_z_1(1:4,1:count_1),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)    %2021-08-02
  patch(c_x_2(1:4,1:count_2),c_y_2(1:4,1:count_2),c_z_2(1:4,1:count_2),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_3(1:4,1:count_3),c_y_3(1:4,1:count_3),c_z_3(1:4,1:count_3),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_4(1:4,1:count_4),c_y_4(1:4,1:count_4),c_z_4(1:4,1:count_4),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_5(1:4,1:count_5),c_y_5(1:4,1:count_5),c_z_5(1:4,1:count_5),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)
  patch(c_x_6(1:4,1:count_6),c_y_6(1:4,1:count_6),c_z_6(1:4,1:count_6),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face)     
end

%===========================  
%  绘制裂缝面及其他附属量
%===========================  
if Key_PLOT(1,5) >= 1 || (Key_PLOT(1,5) == 0  && (Key_PLOT(1,6) + Key_PLOT(1,7))>0)
    disp(['      ----- Plotting crack surface...'])
	if isempty(Crack_X)==0
		for i = 1:num_Crack(isub)
			if i==1
				c_clor = 'green';
			elseif i==2
				c_clor = 'red';
			elseif i==3
				c_clor = 'blue';
			elseif i==4	
				c_clor = 'cyan';
			elseif i==5	
				c_clor = 'magenta';
			elseif i==6	
				c_clor = 'yellow';
			end
			nPt = size(Crack_X{i},2);
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
			%ooooooooooooooooooooooooooo
			%    绘制裂缝面离散点
			%ooooooooooooooooooooooooooo
			if Key_PLOT(1,5) >= 2
				nnode = size(Crack_node_X{i},2);
				for j=1:nnode
					c_node_x = [Crack_node_X{i}(j)];
					c_node_y = [Crack_node_Y{i}(j)];
					c_node_z = [Crack_node_Z{i}(j)];
					plot3(c_node_x,c_node_y,c_node_z,'c.','MarkerSize',16.0)    % MarkerSize 表示点的大小,黑色点
					%裂缝面离散点编号
					if Key_PLOT(1,5) >= 3
					    ts = text(c_node_x,c_node_y,c_node_z,num2str(j),'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');
					end
				end
			end
			%ooooooooooooooooooooooooooo
			%        绘制单元面
			%ooooooooooooooooooooooooooo
			if Key_PLOT(1,5) >= 1
				nele = size(Crack_Ele_1{i},2);
				c_x =[];c_y =[];c_z =[];
				for j=1:nele
					c_x(j,1:3) = [Crack_node_X{i}(Crack_Ele_1{i}(j)),Crack_node_X{i}(Crack_Ele_2{i}(j)),Crack_node_X{i}(Crack_Ele_3{i}(j))];
					c_y(j,1:3) = [Crack_node_Y{i}(Crack_Ele_1{i}(j)),Crack_node_Y{i}(Crack_Ele_2{i}(j)),Crack_node_Y{i}(Crack_Ele_3{i}(j))];
					c_z(j,1:3) = [Crack_node_Z{i}(Crack_Ele_1{i}(j)),Crack_node_Z{i}(Crack_Ele_2{i}(j)),Crack_node_Z{i}(Crack_Ele_3{i}(j))];					
				end		
				patch(c_x',c_y',c_z',[1,1,0],'FaceAlpha',0.5,'FaceLighting','gouraud')			    %橘黄色,不透明,黄色边线							
			end
			%绘制单元面
			% nele = size(Crack_Ele_1{i},2);
			% c_x =[];c_y =[];c_z =[];
			% for j=1:nele
				% c_x(j,1:3) = [Crack_node_X_new{i}(Crack_Ele_1{i}(j)),Crack_node_X_new{i}(Crack_Ele_2{i}(j)),Crack_node_X_new{i}(Crack_Ele_3{i}(j))];
				% c_y(j,1:3) = [Crack_node_Y_new{i}(Crack_Ele_1{i}(j)),Crack_node_Y_new{i}(Crack_Ele_2{i}(j)),Crack_node_Y_new{i}(Crack_Ele_3{i}(j))];
				% c_z(j,1:3) = [Crack_node_Z_new{i}(Crack_Ele_1{i}(j)),Crack_node_Z_new{i}(Crack_Ele_2{i}(j)),Crack_node_Z_new{i}(Crack_Ele_3{i}(j))];					
			% end		
            % patch(c_x',c_y',c_z',[1,1,0],'FaceAlpha',0.5,'FaceLighting','gouraud')	
			
			%ooooooooooooooooooooooooooo
			%绘制流体单元计算点(流体节点)
			%ooooooooooooooooooooooooooo
			if Key_PLOT(1,6) == 1  || Key_PLOT(1,6) == 2
			    if isempty(Crack3D_CalP_X)==0
					nCalP = size(Crack3D_CalP_X{i},2);
					for j=1:nCalP
						c_CalP_x = Crack3D_CalP_X{i}(j);
						c_CalP_y = Crack3D_CalP_Y{i}(j);
						c_CalP_z = Crack3D_CalP_Z{i}(j);
						plot3(c_CalP_x,c_CalP_y,c_CalP_z,'k.','MarkerSize',6.0)    % MarkerSize 表示点的大小,黑色点;k表示黑色.
						%计算点编号
						if  Key_PLOT(1,6) == 2
						    ts = text(c_CalP_x,c_CalP_y,c_CalP_z,num2str(j),'Color','black','FontSize',8,'FontName','Consolas','FontAngle','italic');
						end
					end
					nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
					for j=1:nfluid_Ele
						c_num_fluid_nodes = 3;
						c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
						c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
						c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
						c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
						% fill3(c_x,c_y,c_z,c_clor,'FaceAlpha',0.1,'FaceLighting','gouraud')	
						patch(c_x,c_y,c_z,c_clor,'FaceAlpha',0.1,'FaceLighting','gouraud')	
						% fill3(c_x,c_y,c_z,'b')	
						% fill3(c_x,c_y,c_z,[1*double(j/nfluid_Ele),1*(1-double(j/nfluid_Ele)),1*double(j/nfluid_Ele)])			
						% fill3(c_x,c_y,c_z,[rand(1),rand(1),rand(1)])	%颜色控制
						% 流体单元编号
						if Key_PLOT(1,6) == 4 || Key_PLOT(1,6) == 5
							flu_ele_ave_x = sum(c_x)/c_num_fluid_nodes;
							flu_ele_ave_y = sum(c_y)/c_num_fluid_nodes;
							flu_ele_ave_z = sum(c_z)/c_num_fluid_nodes;
							ts = text(flu_ele_ave_x,flu_ele_ave_y,flu_ele_ave_z,num2str(j),'Color','black','FontSize',11,'FontName','Consolas','FontAngle','italic'); % 流体单元编号
						end					
                    end					
				end
			end
			%ooooooooooooooooooooooooooooooooooooooo
			% 流体单元高斯点外法线向量,2020-01-22
			%ooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,7) == 21
			    nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
				for j=1:nfluid_Ele
					c_num_fluid_nodes = 3;
					c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
					c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
					c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
					c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
					flu_ele_ave_x = sum(c_x)/c_num_fluid_nodes;
					flu_ele_ave_y = sum(c_y)/c_num_fluid_nodes;
					flu_ele_ave_z = sum(c_z)/c_num_fluid_nodes;
					%绘制流体单元Gauss点
					plot3(flu_ele_ave_x,flu_ele_ave_y,flu_ele_ave_z,'black.','MarkerSize',8.0)    % MarkerSize 表示点的大小,黑色点	
					%获得法向量
					c_FlEl_Nor_x = FluidEle_GaussNor_3D_X{i}(j);
					c_FlEl_Nor_y = FluidEle_GaussNor_3D_Y{i}(j);
					c_FlEl_Nor_z = FluidEle_GaussNor_3D_Z{i}(j);
					c_Length = (c_X_Length+c_Y_Length+c_Z_Length)/50;  %法向量的长度
					c_delta_x = c_Length*c_FlEl_Nor_x;
					c_delta_y = c_Length*c_FlEl_Nor_y;
					c_delta_z = c_Length*c_FlEl_Nor_z;
					h = Tools_mArrow3([flu_ele_ave_x flu_ele_ave_y flu_ele_ave_z], ...
					                  [flu_ele_ave_x+c_delta_x flu_ele_ave_y+c_delta_y flu_ele_ave_z+c_delta_z],...
									  'color','blue','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);					
				end					
			end		
			%oooooooooooooooooooooooooooooooooooooooo	
			% 流体单元高斯点局部坐标系,2020-01-22
			%oooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,7) == 22
			    nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
				for j=1:nfluid_Ele
					c_num_fluid_nodes = 3;
					c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
					c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
					c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
					c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
					flu_ele_ave_x = sum(c_x)/c_num_fluid_nodes;
					flu_ele_ave_y = sum(c_y)/c_num_fluid_nodes;
					flu_ele_ave_z = sum(c_z)/c_num_fluid_nodes;
					%绘制流体单元Gauss点
					plot3(flu_ele_ave_x,flu_ele_ave_y,flu_ele_ave_z,'black.','MarkerSize',8.0)    % MarkerSize 表示点的大小,黑色点	
					%坐标系的长度
					c_Length = (c_X_Length+c_Y_Length+c_Z_Length)/60;  
					%局部坐标系X轴
					c_FlEl_LCS_X_X = FluidEle_LCS_VectorX_X{i}(j);
					c_FlEl_LCS_X_Y = FluidEle_LCS_VectorX_Y{i}(j);
					c_FlEl_LCS_X_Z = FluidEle_LCS_VectorX_Z{i}(j);
					c_delta_x = c_Length*c_FlEl_LCS_X_X;
					c_delta_y = c_Length*c_FlEl_LCS_X_Y;
					c_delta_z = c_Length*c_FlEl_LCS_X_Z;
					h = Tools_mArrow3([flu_ele_ave_x flu_ele_ave_y flu_ele_ave_z], ...
					                  [flu_ele_ave_x+c_delta_x flu_ele_ave_y+c_delta_y flu_ele_ave_z+c_delta_z],...
									  'color','red','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);		
					%局部坐标系Y轴
					c_FlEl_LCS_Y_X = FluidEle_LCS_VectorY_X{i}(j);
					c_FlEl_LCS_Y_Y = FluidEle_LCS_VectorY_Y{i}(j);
					c_FlEl_LCS_Y_Z = FluidEle_LCS_VectorY_Z{i}(j);
					c_delta_x = c_Length*c_FlEl_LCS_Y_X;
					c_delta_y = c_Length*c_FlEl_LCS_Y_Y;
					c_delta_z = c_Length*c_FlEl_LCS_Y_Z;
					h = Tools_mArrow3([flu_ele_ave_x flu_ele_ave_y flu_ele_ave_z], ...
					                  [flu_ele_ave_x+c_delta_x flu_ele_ave_y+c_delta_y flu_ele_ave_z+c_delta_z],...
									  'color','green','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);	
					%局部坐标系Z轴
					c_FlEl_LCS_Z_X = FluidEle_LCS_VectorZ_X{i}(j);
					c_FlEl_LCS_Z_Y = FluidEle_LCS_VectorZ_Y{i}(j);
					c_FlEl_LCS_Z_Z = FluidEle_LCS_VectorZ_Z{i}(j);
					c_delta_x = c_Length*c_FlEl_LCS_Z_X;
					c_delta_y = c_Length*c_FlEl_LCS_Z_Y;
					c_delta_z = c_Length*c_FlEl_LCS_Z_Z;
					h = Tools_mArrow3([flu_ele_ave_x flu_ele_ave_y flu_ele_ave_z], ...
					                  [flu_ele_ave_x+c_delta_x flu_ele_ave_y+c_delta_y flu_ele_ave_z+c_delta_z],...
									  'color','blue','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);									  
				end					
			end		
			
			%oooooooooooooooooooooooooooooooooooooooooooo	
			% 流体单元(接触单元)高斯点接触力,2020-01-22
			%oooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,7) == 25
			    if isempty(Crack3D_CalP_X)==0 && isempty(FluidEle_Contact_Force_X)==0
					nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
					max_Force_X = max(abs(FluidEle_Contact_Force_X{i}));
					max_Force_Y = max(abs(FluidEle_Contact_Force_Y{i}));
					max_Force_Z = max(abs(FluidEle_Contact_Force_Z{i}));
					max_Force = sqrt(max_Force_X^2 + max_Force_Y^2 + max_Force_Z^2);
					for j=1:nfluid_Ele
						c_num_fluid_nodes = 3;
						c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
						c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
						c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
						c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
						flu_ele_ave_x = sum(c_x)/c_num_fluid_nodes;
						flu_ele_ave_y = sum(c_y)/c_num_fluid_nodes;
						flu_ele_ave_z = sum(c_z)/c_num_fluid_nodes;
						%绘制流体单元Gauss点
						plot3(flu_ele_ave_x,flu_ele_ave_y,flu_ele_ave_z,'black.','MarkerSize',8.0)    % MarkerSize 表示点的大小,黑色点	
						%坐标系的长度
						c_Length = (c_X_Length+c_Y_Length+c_Z_Length)/25;  
						c_Contact_Force_X = FluidEle_Contact_Force_X{i}(j);
						c_Contact_Force_Y = FluidEle_Contact_Force_Y{i}(j);
						c_Contact_Force_Z = FluidEle_Contact_Force_Z{i}(j);
						c_delta_x = c_Length*c_Contact_Force_X/max_Force;
						c_delta_y = c_Length*c_Contact_Force_Y/max_Force;
						c_delta_z = c_Length*c_Contact_Force_Z/max_Force;
						h = Tools_mArrow3([flu_ele_ave_x+c_delta_x flu_ele_ave_y+c_delta_y flu_ele_ave_z+c_delta_z], ...
										  [flu_ele_ave_x flu_ele_ave_y flu_ele_ave_z],...
										  'color','red','stemWidth',c_Length/35.0,'tipWidth',c_Length/15.0,'facealpha',1.0);						
					end	
                end				
			end								
			%oooooooooooooooooooooooooooooooooooooooo					
			% 绘制流体单元计算点开度(矢量表示)
			%oooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,7) == 4 
			    if isempty(Crack3D_CalP_Aperture)==0 && isempty(Cracks_FluNode_Vector_3D_X)==0 && isempty(Crack3D_CalP_X)==0
				    max_length = (c_X_Length+c_Y_Length+c_Z_Length)/30;
					Max_Aperture =  max(abs(Max_Aperture_of_each_Crack));
					nCalP = size(Crack3D_CalP_X{i},2);
					for j=1:nCalP
						c_CalP_x = Crack3D_CalP_X{i}(j);
						c_CalP_y = Crack3D_CalP_Y{i}(j);
						c_CalP_z = Crack3D_CalP_Z{i}(j);
						c_CalP_V_x = Cracks_FluNode_Vector_3D_X{i}(j);
						c_CalP_V_y = Cracks_FluNode_Vector_3D_Y{i}(j);
						c_CalP_V_z = Cracks_FluNode_Vector_3D_Z{i}(j);
						c_Calp_Apeture = Crack3D_CalP_Aperture{i}(j);
						c_Length = c_Calp_Apeture/Max_Aperture*max_length;
						c_delta_x = c_Length*c_CalP_V_x;
						c_delta_y = c_Length*c_CalP_V_y;
						c_delta_z = c_Length*c_CalP_V_z;
						% plot3(c_CalP_x,c_CalP_y,c_CalP_z,'g.','MarkerSize',12.0)    % MarkerSize 表示点的大小,黑色点
						% ts = text(c_CalP_x,c_CalP_y,c_CalP_z,num2str(j),'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic');
						if c_Calp_Apeture >=0 
						    %张开型裂缝箭头朝外
						    %张开型裂缝用黑色
						    h = Tools_mArrow3([c_CalP_x c_CalP_y c_CalP_z],[c_CalP_x+c_delta_x c_CalP_y+c_delta_y c_CalP_z+c_delta_z],...
							'color','black','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);
						    h = Tools_mArrow3([c_CalP_x c_CalP_y c_CalP_z],[c_CalP_x-c_delta_x c_CalP_y-c_delta_y c_CalP_z-c_delta_z],...
							'color','black','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);
						else
						    %闭合型裂缝箭头朝内(指向裂缝面)
						    h = Tools_mArrow3([c_CalP_x+c_delta_x c_CalP_y+c_delta_y c_CalP_z+c_delta_z],[c_CalP_x c_CalP_y c_CalP_z],...
							'color','black','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);
						    h = Tools_mArrow3([c_CalP_x-c_delta_x c_CalP_y-c_delta_y c_CalP_z-c_delta_z],[c_CalP_x c_CalP_y c_CalP_z],...
							'color','black','stemWidth',c_Length/25.0,'tipWidth',c_Length/10.0,'facealpha',1.0);
                        end							
						
					end
				end
			end			
            %oooooooooooooooooooooooooooooooooooooooo
			% 绘制流体单元开度(云图)
			%oooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,7) == 14
			    if isempty(Crack3D_FluEl_Aperture)==0
				    Max_Aper = max(Max_Aperture_of_each_Crack(i));
					Min_Aper = min(Min_Aperture_of_each_Crack(i));
					nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
					for j=1:nfluid_Ele
					    c_Aperture = Crack3D_FluEl_Aperture{i}(j);
						c_num_fluid_nodes = 3;
						c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
						c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
						c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
						c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
						c_Color_Value = (c_Aperture-Min_Aper)/(Max_Aper-Min_Aper);
						
						if c_Color_Value >=0.0 && c_Color_Value <=0.1 
						    patch(c_x,c_y,c_z,[0/255,0/255,253/255])	
                        elseif c_Color_Value >0.1 && c_Color_Value <=0.2 
						    patch(c_x,c_y,c_z,[1/255,173/255,255/255])	
                        elseif c_Color_Value >0.2 && c_Color_Value <=0.3 
						    patch(c_x,c_y,c_z,[0/255,255/255,255/255])	
                        elseif c_Color_Value >0.3 && c_Color_Value <=0.4 
						    patch(c_x,c_y,c_z,[0/255,255/255,173/255])	
                        elseif c_Color_Value >0.4 && c_Color_Value <=0.5 
						    patch(c_x,c_y,c_z,[0/255,255/255,85/255])	
                        elseif c_Color_Value >0.5 && c_Color_Value <=0.6 
						    patch(c_x,c_y,c_z,[84/255,255/255,0/255])	
                        elseif c_Color_Value >0.6 && c_Color_Value <=0.7 
						    patch(c_x,c_y,c_z,[173/255,255/255,0/255])	
                        elseif c_Color_Value >0.7 && c_Color_Value <=0.8 
						    patch(c_x,c_y,c_z,[255/255,255/255,0/255])	
                        elseif c_Color_Value >0.8 && c_Color_Value <=0.9 
						    patch(c_x,c_y,c_z,[255/255,173/255,1/255])	
                        elseif c_Color_Value >0.9 && c_Color_Value <=1.0 
						    patch(c_x,c_y,c_z,[255/255,0/255,0/255])								
                        end						
					end
				end
			end	

			%oooooooooooooooooooooooooooooooooooooooo
			%绘制流体单元(彩色绘制)
			%oooooooooooooooooooooooooooooooooooooooo
			if Key_PLOT(1,6) == 3 || Key_PLOT(1,6) == 4 || Key_PLOT(1,6) == 5
			    if isempty(Crack3D_CalP_X)==0
					nCalP = size(Crack3D_CalP_X{i},2);
					for j=1:nCalP
						c_CalP_x = [Crack3D_CalP_X{i}(j)];
						c_CalP_y = [Crack3D_CalP_Y{i}(j)];
						c_CalP_z = [Crack3D_CalP_Z{i}(j)];
						if Key_PLOT(1,6) ~= 3 
						    plot3(c_CalP_x,c_CalP_y,c_CalP_z,'g.','MarkerSize',10.0)    % MarkerSize 表示点的大小,黑色点
						end
						%ts = text(c_CalP_x,c_CalP_y,c_CalP_z,num2str(j),'Color','blue','FontSize',12,'FontName','Consolas','FontAngle','italic'); % 流体单元节点编号
					end
				end
				nfluid_Ele = Crack3D_Fluid_Ele_Num{i};
				for j=1:nfluid_Ele
				    c_num_fluid_nodes = 3;
				    c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i,j,1:c_num_fluid_nodes);
					c_x =Crack3D_CalP_X{i}(c_fluid_nodes);
					c_y =Crack3D_CalP_Y{i}(c_fluid_nodes);
					c_z =Crack3D_CalP_Z{i}(c_fluid_nodes);
					patch(c_x,c_y,c_z,'blue','FaceAlpha',0.1,'FaceLighting','gouraud')	
					% fill3(c_x,c_y,c_z,'b')	
                    % fill3(c_x,c_y,c_z,[1*double(j/nfluid_Ele),1*(1-double(j/nfluid_Ele)),1*double(j/nfluid_Ele)])			
                    % fill3(c_x,c_y,c_z,[rand(1),rand(1),rand(1)])	%颜色控制
					
                    % 流体单元编号
                    if Key_PLOT(1,6) == 4 || Key_PLOT(1,6) == 5
					    flu_ele_ave_x = sum(c_x)/c_num_fluid_nodes;
						flu_ele_ave_y = sum(c_y)/c_num_fluid_nodes;
						flu_ele_ave_z = sum(c_z)/c_num_fluid_nodes;
						ts = text(flu_ele_ave_x,flu_ele_ave_y,flu_ele_ave_z,num2str(j),'Color','blue','FontSize',10,'FontName','Consolas','FontAngle','italic'); % 流体单元编号
                    end								
				end
				%绘制计算点及编号
				if  Key_PLOT(1,6) == 5
					nCalP = size(Crack3D_CalP_X{i},2);
					for j=1:nCalP
						c_CalP_x = Crack3D_CalP_X{i}(j);
						c_CalP_y = Crack3D_CalP_Y{i}(j);
						c_CalP_z = Crack3D_CalP_Z{i}(j);
						plot3(c_CalP_x,c_CalP_y,c_CalP_z,'g.','MarkerSize',10.0)    % MarkerSize 表示点的大小,黑色点
						ts = text(c_CalP_x,c_CalP_y,c_CalP_z,num2str(j),'Color','red','FontSize',10,'FontName','Consolas','FontAngle','italic');
					end
				end					
			end				
		end	
	end
end
                        
% Plot the node numbers.
if Key_PLOT(1,2) ==2
    disp(['      ----- Plotting node number...'])
	tem_string = 1:1:Num_Node;
    text(Node_Coor(:,1),Node_Coor(:,2),Node_Coor(:,3),string(tem_string),'FontName',Title_Font,'FontSize',9,'color',[3/255,168/255,158/255])	%锰蓝色
end

% Plot the element numbers.
if Key_PLOT(1,3) >=2
    disp(['      ----- Plotting element number...'])
    for iElem = 1:Num_Elem
        NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,5) Elem_Node(iElem,6) Elem_Node(iElem,7) Elem_Node(iElem,8)];
        XN = Node_Coor(NN,1);
        YN = Node_Coor(NN,2);
		ZN = Node_Coor(NN,3);
		% mean_XN(iElem) = mean(XN)+rand()*0.001;
		% mean_YN(iElem) = mean(YN)+rand()*0.001;
		% mean_ZN(iElem) = mean(ZN)+rand()*0.001;
		mean_XN(iElem) = mean(XN);
		mean_YN(iElem) = mean(YN);
		mean_ZN(iElem) = mean(ZN);		
		tem_string = 1:1:Num_Elem;
    end
    text(mean_XN,mean_YN,mean_ZN,string(tem_string),'FontName',Title_Font,'FontSize',10,'color','black')	
	%Plot element center
	if Key_PLOT(1,3) >=3
	    plot3(mean_XN,mean_YN,mean_ZN,'k*','MarkerSize',9.0)    % MarkerSize 表示点的大小,黑色点
	end
end

% Plot nodes
if isempty(Node_Coor) ~= 1
    length_min = min(c_X_Length,min(c_Y_Length,c_Z_Length))/5.0; 
	r = length_min/20;
    if Key_PLOT(1,2)>=1
		disp(['      ----- Plotting nodes......'])
		x_node =[];y_node =[];z_node =[];
		count = 0;
		for j =1:Num_Node
			count = count +1;
			x_node(count) = Node_Coor(j,1);                                          
			y_node(count) = Node_Coor(j,2);  
			z_node(count) = Node_Coor(j,3);  		
		end
		plot3(x_node,y_node,z_node,'k.','MarkerSize',9.0)    % MarkerSize 表示点的大小,黑色点
	end
end


% Plot enriched nodes
if isempty(Post_Enriched_Nodes) ~= 1
    length_min = min(c_X_Length,min(c_Y_Length,c_Z_Length))/5.0; 
	r = length_min/10;
    if Key_PLOT(1,8)==1
		disp(['      ----- Plotting enriched nodes......'])
		count_1 = 0;count_2 = 0;count_3 = 0;
		x_node_1 =[];y_node_1 =[];z_node_1 =[];		
		x_node_2 =[];y_node_2 =[];z_node_2 =[];	
		x_node_3 =[];y_node_3 =[];z_node_3 =[];	
		for i =1 :size(Post_Enriched_Nodes,2)
			for j =1:Num_Node
				x_node = Node_Coor(j,1);                                          
				y_node = Node_Coor(j,2);  
                z_node = Node_Coor(j,3);  				
				if Post_Enriched_Nodes(j,i)==1     % Tip nodes
				    count_1 = count_1 +1;
				    x_node_1(count_1) = Node_Coor(j,1);                                          
				    y_node_1(count_1) = Node_Coor(j,2);  
                    z_node_1(count_1) = Node_Coor(j,3);  						
				elseif Post_Enriched_Nodes(j,i)==2 % Heaviside nodes
				    count_2 = count_2 +1;
				    x_node_2(count_2) = Node_Coor(j,1);                                          
				    y_node_2(count_2) = Node_Coor(j,2);  
                    z_node_2(count_2) = Node_Coor(j,3);  					
				elseif Post_Enriched_Nodes(j,i)==3 % Junction nodes
				    count_3 = count_3 +1;
				    x_node_3(count_3) = Node_Coor(j,1);                                          
				    y_node_3(count_3) = Node_Coor(j,2);  
                    z_node_3(count_3) = Node_Coor(j,3);  					
				end
			end
		end
		plot3(x_node_1,y_node_1,z_node_1,'r.','MarkerSize',14.0)  		
		plot3(x_node_2,y_node_2,z_node_2,'b.','MarkerSize',12.0) 
		plot3(x_node_3,y_node_3,z_node_3,'g.','MarkerSize',14.0) 
	end
end

%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
%裂尖增强单元的基准线(baseline)和基准线的局部坐标系,2020-01-21
%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
% for iElem = 1:Num_Elem
	% if sum(abs(Tip_Enriched_Ele_BaseLine(iElem,1:6)))>1.0D-10
		% Point_A = Tip_Enriched_Ele_BaseLine(iElem,1:3);
		% Point_B = Tip_Enriched_Ele_BaseLine(iElem,4:6);
		% plot3([Point_A(1),Point_B(1)],[Point_A(2),Point_B(2)],[Point_A(3),Point_B(3)],'LineWidth',1.5,'Color','m')	
		% plot3(Point_A(1),Point_A(2),Point_A(3),'r.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点 
		% plot3(Point_B(1),Point_B(2),Point_B(3),'r.','MarkerSize',16.0);    % MarkerSize 表示点的大小,青色点 
        % StartPoint = Point_A;
	    % EndPoint   = Point_B;
		% h = Tools_mArrow3(StartPoint,EndPoint,'color','black','stemWidth',0.018,'facealpha',1.0);
	% end		
% end
		
% Plot Gauss points.
% if Key_PLOT(1,4) == 1
    % disp(['      ----- Plotting Gauss points...'])
    %%%% Read gauss point coordinates file.
	% Gauss_Coor = load([Full_Pathname,'.gcor_',num2str(isub)]);
	% plot(Gauss_Coor(:,2),Gauss_Coor(:,3),'bo','MarkerSize',1,'Color','black')

	% clear Gauss_Coor
% end



% Active Figure control widget (2021-08-01)
% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
% Press q to exit.
% Press r (or double-click) to reset to the initial.
if Key_Figure_Control_Widget==1
    fcw(gca);
end

% Save pictures.
Save_Picture(c_figure,Full_Pathname,'mesh')
