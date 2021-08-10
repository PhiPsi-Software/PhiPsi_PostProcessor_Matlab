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

function Plot_Crack_Contour_3D(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
%Plot 3D crack contour.
%2020-03-12

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
global Crack3D_Node_Aperture 
global Key_Flipped_Gray Color_Contourlevel
global Title_Font Key_Figure_Control_Widget
global Crack3D_Fluid_Ele_Nodes Crack3D_CalP_Aperture Crack3D_Fluid_Ele_Num

% global Tri_BCD
disp(['      ----- Plotting crack contour......'])
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
title('Crack contour','FontName',Title_Font,'FontSize',Size_Font)
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

%----------------------
%     绘制坐标轴
%----------------------
Arrow_length = (c_X_Length+c_Y_Length+c_Z_Length)/15;
h = Tools_mArrow3([0 0 0],[Arrow_length 0 0],'color','red','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 Arrow_length 0],'color','green','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
h = Tools_mArrow3([0 0 0],[0 0 Arrow_length],'color','blue','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');

%----------------------
%   绘制变形后的网格
%----------------------
%绘制变形后的网格
Line_width =0.1;
if Key_PLOT(5,9)==1 
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
	plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',Line_width,'Color','red')	%灰色				
end

if Key_PLOT(5,9)==2  | Key_PLOT(5,9)==4 %增强单元的网格
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

if Key_PLOT(5,9)==3  | Key_PLOT(5,9)==4  %模型的表面网格	
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
if Key_PLOT(5,9) ==0 %仅绘制模型边框(Outlines)
	for i=1:size(Model_Outline,1)
		plot3([New_Node_Coor(Model_Outline(i,1),1),New_Node_Coor(Model_Outline(i,2),1)],...
			  [New_Node_Coor(Model_Outline(i,1),2),New_Node_Coor(Model_Outline(i,2),2)],...
			  [New_Node_Coor(Model_Outline(i,1),3),New_Node_Coor(Model_Outline(i,2),3)],'LineWidth',Line_width,'Color','black')	    
	end
end




%-----------------------------
%Plot crack aperture contour.
%-----------------------------
if Key_PLOT(5,1)==1 
	if Key_PLOT(5,2)==1 %results of fluid elements and nodes
		if isempty(Crack3D_CalP_Aperture)==0  && isempty(Crack3D_CalP_X)==0
			for i_Crack = 1:num_Crack(isub)    
				nfluid_Ele = Crack3D_Fluid_Ele_Num{i_Crack};
				n_fl_nodes_per_elem = 3;   %Every fluid element has 3 fluid nodes.
				% Initialization of the required matrices
				X = zeros(n_fl_nodes_per_elem,nfluid_Ele);
				Y = zeros(n_fl_nodes_per_elem,nfluid_Ele);
				Z = zeros(n_fl_nodes_per_elem,nfluid_Ele);
				for i_fluid_Ele=1:nfluid_Ele
					c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i_Crack,i_fluid_Ele,1:n_fl_nodes_per_elem);
					X(1:n_fl_nodes_per_elem,i_fluid_Ele) =Crack3D_CalP_X{i_Crack}(c_fluid_nodes);
					Y(1:n_fl_nodes_per_elem,i_fluid_Ele) =Crack3D_CalP_Y{i_Crack}(c_fluid_nodes);
					Z(1:n_fl_nodes_per_elem,i_fluid_Ele) =Crack3D_CalP_Z{i_Crack}(c_fluid_nodes);
					c_X = Crack3D_CalP_X{i_Crack}(c_fluid_nodes);
					c_Y = Crack3D_CalP_Y{i_Crack}(c_fluid_nodes);
					c_Z = Crack3D_CalP_Z{i_Crack}(c_fluid_nodes);
					profile(1:n_fl_nodes_per_elem,i_fluid_Ele) =  Crack3D_CalP_Aperture{i_Crack}(c_fluid_nodes)*1000.0 ; % extract component value of the fluid node 				
				end	
				patch(X,Y,Z,profile)
			end
		end
    elseif Key_PLOT(5,2)==2    %results of discrete crack surface
		if isempty(Crack_X)==0
			for i = 1:num_Crack(isub)    
				nel   = size(Crack_Ele_1{i},2);    % number of elements
				nnode = size(Crack_node_X{i},2);   % total number of nodes in system
				nnel  = 3;                         % number of nodes per element
				% Initialization of the required matrices
				X = zeros(nnel,nel);
				Y = zeros(nnel,nel);
				Z = zeros(nnel,nel);
				profile = zeros(nnel,nel) ;
				for iel=1:nel
					nd=[Crack_Ele_1{i}(iel) Crack_Ele_2{i}(iel) Crack_Ele_3{i}(iel)];         % extract connected node for (iel)-th element
					X(1:nnel,iel)=Crack_node_X{i}(nd);    % extract x value of the node
					Y(1:nnel,iel)=Crack_node_Y{i}(nd);    % extract y value of the node
					Z(1:nnel,iel)=Crack_node_Z{i}(nd) ;   % extract z value of the node
					profile(1:nnel,iel) = Crack3D_Node_Aperture{i}(nd)*1000.0 ; % extract component value of the node 		
				end	
				patch(X,Y,Z,profile)
			end
		end
	end
end	

%---------------------------------------
%Plot crack node and element numnber.
%---------------------------------------
% element
if Key_PLOT(5,5)>=2 
	if Key_PLOT(5,2)==1 %results of fluid elements and nodes
		if isempty(Crack3D_CalP_Aperture)==0  && isempty(Crack3D_CalP_X)==0
			for i_Crack = 1:num_Crack(isub)    
				nfluid_Ele = Crack3D_Fluid_Ele_Num{i_Crack};
				for i_fluid_Ele=1:nfluid_Ele
					c_fluid_nodes = Crack3D_Fluid_Ele_Nodes(i_Crack,i_fluid_Ele,1:3);
					c_X = Crack3D_CalP_X{i_Crack}(c_fluid_nodes);
					c_Y = Crack3D_CalP_Y{i_Crack}(c_fluid_nodes);
					c_Z = Crack3D_CalP_Z{i_Crack}(c_fluid_nodes);
                    ts = text(sum(c_X(1:3))/3.0,sum(c_Y(1:3))/3.0,sum(c_Z(1:3))/3.0,num2str(i_fluid_Ele),'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic'); % 流体单元编号						
				end	
			end
		end
    elseif Key_PLOT(5,2)==2    %results of discrete crack surface
		if isempty(Crack_X)==0
			for i = 1:num_Crack(isub)    
				nel   = size(Crack_Ele_1{i},2);    % number of elements
				nnode = size(Crack_node_X{i},2);   % total number of nodes in system
				nnel  = 3;                         % number of nodes per element
				for iel=1:nel
					nd=[Crack_Ele_1{i}(iel) Crack_Ele_2{i}(iel) Crack_Ele_3{i}(iel)];         % extract connected node for (iel)-th element
					c_X(1:nnel)=Crack_node_X{i}(nd);    % extract x value of the node
					c_Y(1:nnel)=Crack_node_Y{i}(nd);    % extract y value of the node
					c_Z(1:nnel)=Crack_node_Z{i}(nd) ;   % extract z value of the node
                    ts = text(sum(c_X(1:3))/3.0,sum(c_Y(1:3))/3.0,sum(c_Z(1:3))/3.0,num2str(iel),'Color','black','FontSize',12,'FontName','Consolas','FontAngle','italic'); % 裂缝面离散单元编号		
				end	
			end
		end
	end
end
% node
if Key_PLOT(5,5)>=1
	if Key_PLOT(5,2)==1 %results of fluid elements and nodes
		if isempty(Crack3D_CalP_Aperture)==0  && isempty(Crack3D_CalP_X)==0
			for i_Crack = 1:num_Crack(isub)    
				nCalP = size(Crack3D_CalP_X{i_Crack},2);
				for j=1:nCalP
					c_CalP_x = Crack3D_CalP_X{i_Crack}(j);
					c_CalP_y = Crack3D_CalP_Y{i_Crack}(j);
					c_CalP_z = Crack3D_CalP_Z{i_Crack}(j);
					plot3(c_CalP_x,c_CalP_y,c_CalP_z,'k.','MarkerSize',10.0)    % MarkerSize 表示点的大小,黑色点;k表示黑色.
					ts = text(c_CalP_x,c_CalP_y,c_CalP_z,num2str(j),'Color','red','FontSize',10,'FontName','Consolas','FontAngle','italic');
				end
			end
		end
    elseif Key_PLOT(5,2)==2    %results of discrete crack surface
		if isempty(Crack_X)==0
		    for i_Crack = 1:num_Crack(isub)    
				nnode = size(Crack_node_X{i_Crack},2);
				for j=1:nnode
					c_node_x = [Crack_node_X{i_Crack}(j)];
					c_node_y = [Crack_node_Y{i_Crack}(j)];
					c_node_z = [Crack_node_Z{i_Crack}(j)];
					plot3(c_node_x,c_node_y,c_node_z,'c.','MarkerSize',10.0)    % MarkerSize 表示点的大小,黑色点
					ts = text(c_node_x,c_node_y,c_node_z,num2str(j),'Color','red','FontSize',10,'FontName','Consolas','FontAngle','italic');
				end
			end
		end
	end
end

%----------------------
% Set colormap.
%----------------------
if Key_Flipped_Gray==0
	colormap(Color_Contourlevel)
elseif Key_Flipped_Gray==1
	colormap(flipud(gray))
end				
					
%----------------------	
% Set colorbar.
%----------------------
cbar = colorbar;
brighten(0.5); 
if Key_PLOT(5,1)==1
    set(get(cbar,'title'),'string','mm');
end
% get the color limits
clim = caxis;
ylim(cbar,[clim(1) clim(2)]);
numpts = 10 ;    % Number of points to be displayed on colorbar
kssv = linspace(clim(1),clim(2),numpts);
set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
for i = 1:numpts
	imep = num2str(kssv(i),'%+3.2E');
	vasu(i) = {imep} ;
end
set(cbar,'YTickLabel',vasu(1:numpts),'FontName',Title_Font,'FontSize',Size_Font);


%------------------------------------------------------	
% Active Figure control widget (2021-08-01)
% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
% Press q to exit.
% Press r (or double-click) to reset to the initial.
%------------------------------------------------------	
if Key_Figure_Control_Widget==1
    fcw(gca);
end

%----------------------
% Save pictures.
%----------------------
Save_Picture(c_figure,Full_Pathname,'ccon')


