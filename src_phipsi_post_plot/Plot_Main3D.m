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

function Plot_Main3D(POST_Substep)
% This function plots mesh and deformed mesh for 3D problems.

global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y Num_Foc_z Foc_z
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush DISP  FORCE_Matrix
global Key_Data_Format
global Crack_node_X Crack_node_Y Crack_node_Z
global Crack_Ele_1 Crack_Ele_2 Crack_Ele_3
global Model_Outline Model_OutArea
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
global Crack3D_Fluid_Ele_Num
global Crack3D_Fluid_Ele_Nodes
global Crack3D_CalP_Aperture
global Cracks_FluNode_Vector_3D_X Cracks_FluNode_Vector_3D_Y Cracks_FluNode_Vector_3D_Z
global Max_Aperture_of_each_Crack Crack3D_FluEl_Aperture
global Min_Aperture_of_each_Crack Stress_Matrix
global Tip_Enriched_Ele_BaseLine Tip_Enriched_Ele_BaseLine_Vector_x 
global Tip_Enriched_Ele_BaseLine_Vector_y Tip_Enriched_Ele_BaseLine_Vector_z
global FluidEle_GaussNor_3D_X FluidEle_GaussNor_3D_Y FluidEle_GaussNor_3D_Z
global FluidEle_LCS_VectorX_X FluidEle_LCS_VectorX_Y FluidEle_LCS_VectorX_Z
global FluidEle_LCS_VectorY_X FluidEle_LCS_VectorY_Y FluidEle_LCS_VectorY_Z
global FluidEle_LCS_VectorZ_X FluidEle_LCS_VectorZ_Y FluidEle_LCS_VectorZ_Z
global FluidEle_Contact_Force_X FluidEle_Contact_Force_Y FluidEle_Contact_Force_Z
global Crack3D_Node_Aperture
global Tip_Enriched_Node_Ref_Element
global Cracks_CalP_UpDis_3D_X Cracks_CalP_UpDis_3D_Y Cracks_CalP_UpDis_3D_Z
global Cracks_CalP_LowDis_3D_X Cracks_CalP_LowDis_3D_Y Cracks_CalP_LowDis_3D_Z
global Title_Font Key_Figure_Control_Widget
global Only_Plot_Mesh Node_Coor Elem_Node Num_Elem
global G_X_NODES G_Y_NODES G_Z_NODES G_NN G_X_Min G_X_Max G_Y_Min G_Y_Max G_Z_Min G_Z_Max

scale = Key_PLOT(2,6);

% Check Key_Figure_Control_Widget (2021-8-02)
if Key_Figure_Control_Widget==1
	disp('    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
	disp('    ATTENTION :: Figure Control Widget is active!') 
	disp('    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
end	

% Get the number of the substep for post process.
if isa(POST_Substep,'double') == 0
    if length(POST_Substep) == 4
		if  lower(POST_Substep)  ==  'last'
			if Key_Dynamic ==1
			    if Output_Freq  ==0
				    POST_Substep = Real_Iteras;
				else
				    POST_Substep = num_Output_Sub;
				end
			elseif Key_Dynamic ==0
			    if Output_Freq==0
				    POST_Substep = Real_Sub;
				else
			    	POST_Substep = num_Output_Sub;
				end
			end
		else
		    disp('    Error :: Unrecognized parameters of *POST_Substep.') 
		    Error_Message
		end
	elseif length(POST_Substep) == 5
		if  lower(POST_Substep) == 'first'
			POST_Substep = 1;
		else 
			disp('    Error :: Unrecognized parameters of *POST_Substep.') 
		    Error_Message
		end
	else
	    disp('    Error :: Unrecognized parameters of *POST_Substep.') 
		Error_Message
	end
end


%若存在，则读取局部加密后的网格节点坐标及单元编号, 2021-06-30.
if exist([Full_Pathname,'.nodr_',num2str(POST_Substep)], 'file') ==2
    disp('    > Reading refined mesh....') 
    Node_Coor   = load([Full_Pathname,'.nodr_',num2str(POST_Substep)]);
	Elem_Info   = load([Full_Pathname,'.eler_',num2str(POST_Substep)]);
	Elem_Node = Elem_Info(:,1:8);
	Elem_Material = Elem_Info(:,9);       % Material number of elements
	num_of_Material = max(Elem_Material); % Number of materials
	% Get the numbers of nodes, elements, boundaries and forces.
	Num_Node = size(Node_Coor,1);
	Num_Elem = size(Elem_Node,1);
	% Initial Element_Info Matrix.
	Elemet_Info = zeros(Num_Elem,10);
	for i=1:Num_Elem
		Elemet_Info(i,1) = i;
	end
	% Initial Node_Info Matrix.
	Node_Info   = zeros(Num_Node,10);
	for i=1:Num_Node
		Node_Info(i,1) = i;
	end
	% Initial "Yes_Checked" Matrix, flag for initiation check.
	Yes_Checked = zeros(Num_Elem,1);	
	% Get the max and min value of node coordinates.
	Max_X_Coor = max(max(Node_Coor(:,1)));
	Min_X_Coor = min(min(Node_Coor(:,1)));
	Max_Y_Coor = max(max(Node_Coor(:,2)));
	Min_Y_Coor = min(min(Node_Coor(:,2)));
	Max_Z_Coor = max(max(Node_Coor(:,3)));
	Min_Z_Coor = min(min(Node_Coor(:,3)));	
	
	% Get the maximum & minimum & average volume.
	% Num_Elem
	for i=1:Num_Elem
		N1  = Elem_Node(i,1);                                                  % Node 1 for current element
		N2  = Elem_Node(i,2);                                                  % Node 2 for current element
		N3  = Elem_Node(i,3);                                                  % Node 3 for current element
		N4  = Elem_Node(i,4);                                                  % Node 4 for current element
		N5  = Elem_Node(i,5);                                                  % Node 1 for current element
		N6  = Elem_Node(i,6);                                                  % Node 2 for current element
		N7  = Elem_Node(i,7);                                                  % Node 3 for current element
		N8  = Elem_Node(i,8);                                                  % Node 4 for current element	
		NN = [N1 N2 N3 N4 N5 N6 N7 N8];                                        % 8 nodes of the current element  
		X_NODES = Node_Coor(NN,1);                                             % Coordinates of the four nodes of the current element
		Y_NODES = Node_Coor(NN,2);
		Z_NODES = Node_Coor(NN,3);
		%--------------------------------
		%Get the volume of each element.
		%--------------------------------
		shp = alphaShape(X_NODES,Y_NODES,Z_NODES); %Matlab内置函数,根据坐标生成几何体
		volume_ele(i) = volume(shp);               %Matlab内置函数,计算生成的几何体的体积
	end

	max_volume_ele = max(volume_ele);
	min_volume_ele = min(volume_ele);
	aveg_volume_ele = mean(volume_ele);	
	
	% Get G_NN,G_X_NODES,G_Y_NODES,G_X_Min,G_X_Max,G_Y_Min,G_Y_Max,G means Global.
	G_NN      = zeros(8,Num_Elem);
	G_X_NODES = zeros(8,Num_Elem);
	G_Y_NODES = zeros(8,Num_Elem);
	G_Z_NODES = zeros(8,Num_Elem);
	G_X_Min   = zeros(1,Num_Elem);
	G_X_Max   = zeros(1,Num_Elem);
	G_Y_Min   = zeros(1,Num_Elem);
	G_Y_Max   = zeros(1,Num_Elem);
	G_Z_Min   = zeros(1,Num_Elem);
	G_Z_Max   = zeros(1,Num_Elem);

	G_NN = [Elem_Node(:,1)';Elem_Node(:,2)';Elem_Node(:,3)';Elem_Node(:,4)';
			Elem_Node(:,5)';Elem_Node(:,6)';Elem_Node(:,7)';Elem_Node(:,8)';];

	G_X_NODES = [Node_Coor(G_NN(1,:),1)';Node_Coor(G_NN(2,:),1)';Node_Coor(G_NN(3,:),1)';Node_Coor(G_NN(4,:),1)';
				 Node_Coor(G_NN(5,:),1)';Node_Coor(G_NN(6,:),1)';Node_Coor(G_NN(7,:),1)';Node_Coor(G_NN(8,:),1)'];
	G_Y_NODES = [Node_Coor(G_NN(1,:),2)';Node_Coor(G_NN(2,:),2)';Node_Coor(G_NN(3,:),2)';Node_Coor(G_NN(4,:),2)';
				 Node_Coor(G_NN(5,:),2)';Node_Coor(G_NN(6,:),2)';Node_Coor(G_NN(7,:),2)';Node_Coor(G_NN(8,:),2)'];
	G_Z_NODES = [Node_Coor(G_NN(1,:),3)';Node_Coor(G_NN(2,:),3)';Node_Coor(G_NN(3,:),3)';Node_Coor(G_NN(4,:),3)';
				 Node_Coor(G_NN(5,:),3)';Node_Coor(G_NN(6,:),3)';Node_Coor(G_NN(7,:),3)';Node_Coor(G_NN(8,:),3)'];
	G_X_Min   = min(G_X_NODES,[],1);
	G_X_Max   = max(G_X_NODES,[],1);
	G_Y_Min   = min(G_Y_NODES,[],1);
	G_Y_Max   = max(G_Y_NODES,[],1);
	G_Z_Min   = min(G_Z_NODES,[],1);
	G_Z_Max   = max(G_Z_NODES,[],1);
end



% Read nodal displacement file.
disp('    > Reading nodal displacement file....') 
% DISP= load([Full_Pathname,'.disn_',num2str(POST_Substep)]);
if exist([Full_Pathname,'.disn_',num2str(POST_Substep)], 'file') ==2 
	if Key_Data_Format==1 
		DISP   = load([Full_Pathname,'.disn_',num2str(POST_Substep)]);
	elseif Key_Data_Format==2  %Binary
		c_file = fopen([Full_Pathname,'.disn_',num2str(POST_Substep)],'rb');
		[cc_DISP,cc_count]   = fread(c_file,inf,'double');
		fclose(c_file);
		%转换成Matlab中的数据格式
		for ccc_i=1:cc_count/3
			DISP(ccc_i,1) = ccc_i;
			DISP(ccc_i,2) = cc_DISP(ccc_i*3-2);
			DISP(ccc_i,3) = cc_DISP(ccc_i*3-1);
			DISP(ccc_i,4) = cc_DISP(ccc_i*3);
		end
	end
end

% Read nodal stress file.
disp('    > Reading nodal stress file....') 
%读取节点应力
if Key_Data_Format==1 
	if exist([Full_Pathname,'.strn_',num2str(POST_Substep)], 'file') ==2 
		Stress_Matrix = load([Full_Pathname,'.strn_',num2str(POST_Substep)]);
	else
		Stress_Matrix =[];
	end
elseif Key_Data_Format==2  %Binary
	if exist([Full_Pathname,'.strn_',num2str(POST_Substep)], 'file') ==2 
		c_file = fopen([Full_Pathname,'.strn_',num2str(POST_Substep)],'rb');
		[cc_Stress_Matrix,cc_count]   = fread(c_file,inf,'double');
		fclose(c_file);
		%转换成Matlab中的数据格式
		for ccc_i=1:cc_count/7
			Stress_Matrix(ccc_i,1) = ccc_i;
			Stress_Matrix(ccc_i,2) = cc_Stress_Matrix(ccc_i*7-6);
			Stress_Matrix(ccc_i,3) = cc_Stress_Matrix(ccc_i*7-5);
			Stress_Matrix(ccc_i,4) = cc_Stress_Matrix(ccc_i*7-4);			
			Stress_Matrix(ccc_i,5) = cc_Stress_Matrix(ccc_i*7-3);
			Stress_Matrix(ccc_i,6) = cc_Stress_Matrix(ccc_i*7-2);
			Stress_Matrix(ccc_i,7) = cc_Stress_Matrix(ccc_i*7-1);
			Stress_Matrix(ccc_i,8) = cc_Stress_Matrix(ccc_i*7);
		end
	else
		Stress_Matrix =[];
	end		
end

% Read force file.
disp('    > Reading force file....') 
Force_X_Matrix = load([Full_Pathname,'.focx']);
Force_Y_Matrix = load([Full_Pathname,'.focy']);
Force_Z_Matrix = load([Full_Pathname,'.focz']);

% Read Boundary file.
disp('    > Reading Boundary file....') 
Boundary_X_Matrix = load([Full_Pathname,'.boux']);
Boundary_Y_Matrix = load([Full_Pathname,'.bouy']);
Boundary_Z_Matrix = load([Full_Pathname,'.bouz']);

% Read outl file.
disp('    > Reading outl (model outlines) file....') 
Model_Outline = load([Full_Pathname,'.outl']);
disp('    > Reading outa (model outareas) file....') 
Model_OutArea = load([Full_Pathname,'.outa']);
% Read coordinates and other info of cracks if cracks exist.
if num_Crack(POST_Substep) ~= 0
	disp('    > Reading coordinates of cracks....') 
	file_Crack_X = fopen([Full_Pathname,'.crax_',num2str(POST_Substep)]);
	file_Crack_Y = fopen([Full_Pathname,'.cray_',num2str(POST_Substep)]);
	file_Crack_Z = fopen([Full_Pathname,'.craz_',num2str(POST_Substep)]);
	Crack_X = cell(num_Crack(POST_Substep));
	Crack_Y = cell(num_Crack(POST_Substep));
	Crack_Z = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack_X{i} = str2num(fgetl(file_Crack_X));
        Crack_Y{i} = str2num(fgetl(file_Crack_Y));
		Crack_Z{i} = str2num(fgetl(file_Crack_Z));
	end
	fclose(file_Crack_X);
	fclose(file_Crack_Y);
	fclose(file_Crack_Z);
	
	% disp('    > Reading ennd file....') 
	Enriched_Node_Type = load([Full_Pathname,'.ennd_',num2str(POST_Substep)]);
	% disp('    > Reading posi file....') 
	POS = load([Full_Pathname,'.posi_',num2str(POST_Substep)]);
	% disp('    > Reading elty file....');
	Elem_Type = load([Full_Pathname,'.elty_',num2str(POST_Substep)]);
	%读取裂缝面节点坐标
	file_Crack_node_X = fopen([Full_Pathname,'.cnox_',num2str(POST_Substep)]);
	file_Crack_node_Y = fopen([Full_Pathname,'.cnoy_',num2str(POST_Substep)]);
	file_Crack_node_Z = fopen([Full_Pathname,'.cnoz_',num2str(POST_Substep)]);
	Crack_node_X = cell(num_Crack(POST_Substep));
	Crack_node_Y = cell(num_Crack(POST_Substep));
	Crack_node_Z = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack_node_X{i} = str2num(fgetl(file_Crack_node_X));
        Crack_node_Y{i} = str2num(fgetl(file_Crack_node_Y));
		Crack_node_Z{i} = str2num(fgetl(file_Crack_node_Z));
	end
	fclose(file_Crack_node_X);
	fclose(file_Crack_node_Y);
	fclose(file_Crack_node_Z);
	%读取裂缝面节点所在单元局部坐标
	file_Crack_node_local_X = fopen([Full_Pathname,'.cnlx_',num2str(POST_Substep)]);
	file_Crack_node_local_Y = fopen([Full_Pathname,'.cnly_',num2str(POST_Substep)]);
	file_Crack_node_local_Z = fopen([Full_Pathname,'.cnlz_',num2str(POST_Substep)]);
	Crack_node_local_X = cell(num_Crack(POST_Substep));
	Crack_node_local_Y = cell(num_Crack(POST_Substep));
	Crack_node_local_Z = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack_node_local_X{i} = str2num(fgetl(file_Crack_node_local_X));
        Crack_node_local_Y{i} = str2num(fgetl(file_Crack_node_local_Y));
		Crack_node_local_Z{i} = str2num(fgetl(file_Crack_node_local_Z));
	end
	fclose(file_Crack_node_local_X);
	fclose(file_Crack_node_local_Y);
	fclose(file_Crack_node_local_Z);	
	%读取裂尖增强单元的基准线
	if exist([Full_Pathname,'.blab_',num2str(POST_Substep)], 'file') ==2 
		Tip_Enriched_Ele_BaseLine = load([Full_Pathname,'.blab_',num2str(POST_Substep)]);
	else
		Tip_Enriched_Ele_BaseLine =[];
	end
	%读取裂尖节点对应的增强单元号(参考单元号)及裂缝号,2021-02-10
	if exist([Full_Pathname,'.tere_',num2str(POST_Substep)], 'file') ==2 
		Tip_Enriched_Node_Ref_Element = load([Full_Pathname,'.tere_',num2str(POST_Substep)]);
	else
		Tip_Enriched_Node_Ref_Element =[];
	end	
    %读取裂尖增强单元的基准线局部坐标系的x轴	
	if exist([Full_Pathname,'.blvx_',num2str(POST_Substep)], 'file') ==2 
		Tip_Enriched_Ele_BaseLine_Vector_x = load([Full_Pathname,'.blvx_',num2str(POST_Substep)]);
	else
		Tip_Enriched_Ele_BaseLine_Vector_x =[];
	end	
    %读取裂尖增强单元的基准线局部坐标系的y轴	
	if exist([Full_Pathname,'.blvy_',num2str(POST_Substep)], 'file') ==2 
		Tip_Enriched_Ele_BaseLine_Vector_y = load([Full_Pathname,'.blvy_',num2str(POST_Substep)]);
	else
		Tip_Enriched_Ele_BaseLine_Vector_y =[];
	end	
    %读取裂尖增强单元的基准线局部坐标系的z轴	
	if exist([Full_Pathname,'.blvz_',num2str(POST_Substep)], 'file') ==2 
		Tip_Enriched_Ele_BaseLine_Vector_z = load([Full_Pathname,'.blvz_',num2str(POST_Substep)]);
	else
		Tip_Enriched_Ele_BaseLine_Vector_z =[];
	end		
	%读取裂缝面离散单元
	file_Crack_Ele_1 = fopen([Full_Pathname,'.cms1_',num2str(POST_Substep)]);
	file_Crack_Ele_2 = fopen([Full_Pathname,'.cms2_',num2str(POST_Substep)]);
	file_Crack_Ele_3 = fopen([Full_Pathname,'.cms3_',num2str(POST_Substep)]);
	Crack_Ele_1 = cell(num_Crack(POST_Substep));
	Crack_Ele_2 = cell(num_Crack(POST_Substep));
	Crack_Ele_3 = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack_Ele_1{i} = str2num(fgetl(file_Crack_Ele_1));
        Crack_Ele_2{i} = str2num(fgetl(file_Crack_Ele_2));
		Crack_Ele_3{i} = str2num(fgetl(file_Crack_Ele_3));
	end
	fclose(file_Crack_Ele_1);
	fclose(file_Crack_Ele_2);
	fclose(file_Crack_Ele_3);	
	%读取裂缝面离散单元节点所在的单元号（模型单元号）
	file_Crack_Node_in_Ele = fopen([Full_Pathname,'.cmse_',num2str(POST_Substep)]);
	Crack_Node_in_Ele = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack_Node_in_Ele{i} = str2num(fgetl(file_Crack_Node_in_Ele));
	end
	fclose(file_Crack_Node_in_Ele);	
	
	%读取离散裂缝面外边界节点
	file_Crack3D_Meshed_Outline= fopen([Full_Pathname,'.cmso_',num2str(POST_Substep)]);
	Crack3D_Meshed_Outline = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack3D_Meshed_Outline{i} = str2num(fgetl(file_Crack3D_Meshed_Outline));
	end
	fclose(file_Crack3D_Meshed_Outline);
	
	%读取裂缝面前缘点最小主应力方向
	if exist([Full_Pathname,'.cndx_',num2str(POST_Substep)], 'file') ==2 
		file_Crack3D_Vector_S1_X = fopen([Full_Pathname,'.cndx_',num2str(POST_Substep)]);
		file_Crack3D_Vector_S1_Y = fopen([Full_Pathname,'.cndy_',num2str(POST_Substep)]);
		file_Crack3D_Vector_S1_Z = fopen([Full_Pathname,'.cndz_',num2str(POST_Substep)]);
		Crack3D_Vector_S1_X = cell(num_Crack(POST_Substep));
		Crack3D_Vector_S1_Y = cell(num_Crack(POST_Substep));
		Crack3D_Vector_S1_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Crack3D_Vector_S1_X{i} = str2num(fgetl(file_Crack3D_Vector_S1_X));
			Crack3D_Vector_S1_Y{i} = str2num(fgetl(file_Crack3D_Vector_S1_Y));
			Crack3D_Vector_S1_Z{i} = str2num(fgetl(file_Crack3D_Vector_S1_Z));
		end
		fclose(file_Crack3D_Vector_S1_X);
		fclose(file_Crack3D_Vector_S1_Y);
		fclose(file_Crack3D_Vector_S1_Z);
	end
	
	%读取裂缝面边界节点的局部坐标系
	file_Crack3D_Vertex_Vector_X_X = fopen([Full_Pathname,'.cvxx_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_X_Y = fopen([Full_Pathname,'.cvxy_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_X_Z = fopen([Full_Pathname,'.cvxz_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Y_X = fopen([Full_Pathname,'.cvyx_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Y_Y = fopen([Full_Pathname,'.cvyy_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Y_Z = fopen([Full_Pathname,'.cvyz_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Z_X = fopen([Full_Pathname,'.cvzx_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Z_Y = fopen([Full_Pathname,'.cvzy_',num2str(POST_Substep)]);
	file_Crack3D_Vertex_Vector_Z_Z = fopen([Full_Pathname,'.cvzz_',num2str(POST_Substep)]);	
	Crack3D_Vertex_Vector_X_X = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_X_Y = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_X_Z = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Y_X = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Y_Y = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Y_Z = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Z_X = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Z_Y = cell(num_Crack(POST_Substep));
	Crack3D_Vertex_Vector_Z_Z = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack3D_Vertex_Vector_X_X{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_X_X));
        Crack3D_Vertex_Vector_X_Y{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_X_Y));
		Crack3D_Vertex_Vector_X_Z{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_X_Z));
		Crack3D_Vertex_Vector_Y_X{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Y_X));
        Crack3D_Vertex_Vector_Y_Y{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Y_Y));
		Crack3D_Vertex_Vector_Y_Z{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Y_Z));
		Crack3D_Vertex_Vector_Z_X{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Z_X));
        Crack3D_Vertex_Vector_Z_Y{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Z_Y));
		Crack3D_Vertex_Vector_Z_Z{i} = str2num(fgetl(file_Crack3D_Vertex_Vector_Z_Z));		
	end
	fclose(file_Crack3D_Vertex_Vector_X_X);  
	fclose(file_Crack3D_Vertex_Vector_X_Y);
	fclose(file_Crack3D_Vertex_Vector_X_Z);  
	fclose(file_Crack3D_Vertex_Vector_Y_X);  
	fclose(file_Crack3D_Vertex_Vector_Y_Y); 
	fclose(file_Crack3D_Vertex_Vector_Y_Z); 
	fclose(file_Crack3D_Vertex_Vector_Z_X); 
	fclose(file_Crack3D_Vertex_Vector_Z_Y);  
	fclose(file_Crack3D_Vertex_Vector_Z_Z);  
	
	%读取裂缝面流体单元数目cpfn文件
	file_cpfn = fopen([Full_Pathname,'.cpfn_',num2str(POST_Substep)]);
	Crack3D_Fluid_Ele_Num = cell(num_Crack(POST_Substep));
	for i=1:num_Crack(POST_Substep) 
		Crack3D_Fluid_Ele_Num{i} = str2num(fgetl(file_cpfn));
	end
	fclose(file_cpfn);
	
	%读取裂缝面流体单元计算点编号文件cpno
	% file_cpno = fopen([Full_Pathname,'.cpno_',num2str(POST_Substep)]);
	% Crack3D_Fluid_Ele_Nodes_All = str2num(fgetl(file_cpno));
	Crack3D_Fluid_Ele_Nodes_All = load([Full_Pathname,'.cpno_',num2str(POST_Substep)]);
	
	cc_count  = 0;
	for i=1:num_Crack(POST_Substep) 
		Crack3D_Fluid_Ele_Nodes(i,1:Crack3D_Fluid_Ele_Num{i},1:3) = Crack3D_Fluid_Ele_Nodes_All(cc_count+1:cc_count+Crack3D_Fluid_Ele_Num{i},1:3);
		cc_count  = cc_count + Crack3D_Fluid_Ele_Num{i};
	end
	
	
	
	
	%读取裂缝面计算点坐标文件
	if exist([Full_Pathname,'.ccpx_',num2str(POST_Substep)], 'file') ==2
		file_Crack3D_CalP_X = fopen([Full_Pathname,'.ccpx_',num2str(POST_Substep)]);
		file_Crack3D_CalP_Y = fopen([Full_Pathname,'.ccpy_',num2str(POST_Substep)]);
		file_Crack3D_CalP_Z = fopen([Full_Pathname,'.ccpz_',num2str(POST_Substep)]);
		Crack3D_CalP_X = cell(num_Crack(POST_Substep));
		Crack3D_CalP_Y = cell(num_Crack(POST_Substep));
		Crack3D_CalP_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Crack3D_CalP_X{i} = str2num(fgetl(file_Crack3D_CalP_X));
			Crack3D_CalP_Y{i} = str2num(fgetl(file_Crack3D_CalP_Y));
			Crack3D_CalP_Z{i} = str2num(fgetl(file_Crack3D_CalP_Z));
		end
		fclose(file_Crack3D_CalP_X);
		fclose(file_Crack3D_CalP_Y);
		fclose(file_Crack3D_CalP_Z);	
    else
	    Crack3D_CalP_X =[];
		Crack3D_CalP_Y =[];
		Crack3D_CalP_Z =[];
	end
	
	%读取裂缝面计算点法向向量文件
	if exist([Full_Pathname,'.fnnx_',num2str(POST_Substep)], 'file') ==2
		file_Cracks_FluNode_Vector_3D_X = fopen([Full_Pathname,'.fnnx_',num2str(POST_Substep)]);
		file_Cracks_FluNode_Vector_3D_Y = fopen([Full_Pathname,'.fnny_',num2str(POST_Substep)]);
		file_Cracks_FluNode_Vector_3D_Z = fopen([Full_Pathname,'.fnnz_',num2str(POST_Substep)]);
		Cracks_FluNode_Vector_3D_X = cell(num_Crack(POST_Substep));
		Cracks_FluNode_Vector_3D_Y = cell(num_Crack(POST_Substep));
		Cracks_FluNode_Vector_3D_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Cracks_FluNode_Vector_3D_X{i} = str2num(fgetl(file_Cracks_FluNode_Vector_3D_X));
			Cracks_FluNode_Vector_3D_Y{i} = str2num(fgetl(file_Cracks_FluNode_Vector_3D_Y));
			Cracks_FluNode_Vector_3D_Z{i} = str2num(fgetl(file_Cracks_FluNode_Vector_3D_Z));
		end
		fclose(file_Cracks_FluNode_Vector_3D_X);
		fclose(file_Cracks_FluNode_Vector_3D_Y);
		fclose(file_Cracks_FluNode_Vector_3D_Z);	
	else
	    Cracks_FluNode_Vector_3D_X=[];
		Cracks_FluNode_Vector_3D_Y=[];
		Cracks_FluNode_Vector_3D_Z=[];
	end

	%读取裂缝面计算点上偏置点的位移向量，2021-02-11.
	if exist([Full_Pathname,'.fnux_',num2str(POST_Substep)], 'file') ==2
		file_Cracks_CalP_UpDis_3D_X = fopen([Full_Pathname,'.fnux_',num2str(POST_Substep)]);
		file_Cracks_CalP_UpDis_3D_Y = fopen([Full_Pathname,'.fnuy_',num2str(POST_Substep)]);
		file_Cracks_CalP_UpDis_3D_Z = fopen([Full_Pathname,'.fnuz_',num2str(POST_Substep)]);
		Cracks_CalP_UpDis_3D_X = cell(num_Crack(POST_Substep));
		Cracks_CalP_UpDis_3D_Y = cell(num_Crack(POST_Substep));
		Cracks_CalP_UpDis_3D_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Cracks_CalP_UpDis_3D_X{i} = str2num(fgetl(file_Cracks_CalP_UpDis_3D_X));
			Cracks_CalP_UpDis_3D_Y{i} = str2num(fgetl(file_Cracks_CalP_UpDis_3D_Y));
			Cracks_CalP_UpDis_3D_Z{i} = str2num(fgetl(file_Cracks_CalP_UpDis_3D_Z));
		end
		fclose(file_Cracks_CalP_UpDis_3D_X);
		fclose(file_Cracks_CalP_UpDis_3D_Y);
		fclose(file_Cracks_CalP_UpDis_3D_Z);	
	else
	    Cracks_CalP_UpDis_3D_X=[];
		Cracks_CalP_UpDis_3D_Y=[];
		Cracks_CalP_UpDis_3D_Z=[];
	end

	%读取裂缝面计算点下偏置点的位移向量，2021-02-11.
	if exist([Full_Pathname,'.fnlx_',num2str(POST_Substep)], 'file') ==2
		file_Cracks_CalP_LowDis_3D_X = fopen([Full_Pathname,'.fnlx_',num2str(POST_Substep)]);
		file_Cracks_CalP_LowDis_3D_Y = fopen([Full_Pathname,'.fnly_',num2str(POST_Substep)]);
		file_Cracks_CalP_LowDis_3D_Z = fopen([Full_Pathname,'.fnlz_',num2str(POST_Substep)]);
		Cracks_CalP_LowDis_3D_X = cell(num_Crack(POST_Substep));
		Cracks_CalP_LowDis_3D_Y = cell(num_Crack(POST_Substep));
		Cracks_CalP_LowDis_3D_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Cracks_CalP_LowDis_3D_X{i} = str2num(fgetl(file_Cracks_CalP_LowDis_3D_X));
			Cracks_CalP_LowDis_3D_Y{i} = str2num(fgetl(file_Cracks_CalP_LowDis_3D_Y));
			Cracks_CalP_LowDis_3D_Z{i} = str2num(fgetl(file_Cracks_CalP_LowDis_3D_Z));
		end
		fclose(file_Cracks_CalP_LowDis_3D_X);
		fclose(file_Cracks_CalP_LowDis_3D_Y);
		fclose(file_Cracks_CalP_LowDis_3D_Z);	
	else
	    Cracks_CalP_LowDis_3D_X=[];
		Cracks_CalP_LowDis_3D_Y=[];
		Cracks_CalP_LowDis_3D_Z=[];
	end

	%读取流体单元平均法向量文件(其实是Gauss点位置的)，2020-01-22
	if exist([Full_Pathname,'.fenx_',num2str(POST_Substep)], 'file') ==2
		file_FluidEle_GaussNor_3D_X = fopen([Full_Pathname,'.fenx_',num2str(POST_Substep)]);
		file_FluidEle_GaussNor_3D_Y = fopen([Full_Pathname,'.feny_',num2str(POST_Substep)]);
		file_FluidEle_GaussNor_3D_Z = fopen([Full_Pathname,'.fenz_',num2str(POST_Substep)]);
		FluidEle_GaussNor_3D_X = cell(num_Crack(POST_Substep));
		FluidEle_GaussNor_3D_Y = cell(num_Crack(POST_Substep));
		FluidEle_GaussNor_3D_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			FluidEle_GaussNor_3D_X{i} = str2num(fgetl(file_FluidEle_GaussNor_3D_X));
			FluidEle_GaussNor_3D_Y{i} = str2num(fgetl(file_FluidEle_GaussNor_3D_Y));
			FluidEle_GaussNor_3D_Z{i} = str2num(fgetl(file_FluidEle_GaussNor_3D_Z));
		end
		fclose(file_FluidEle_GaussNor_3D_X);
		fclose(file_FluidEle_GaussNor_3D_Y);
		fclose(file_FluidEle_GaussNor_3D_Z);	
	else
	    FluidEle_GaussNor_3D_X=[];
		FluidEle_GaussNor_3D_Y=[];
		FluidEle_GaussNor_3D_Z=[];
	end
	
	%读取流体单元Gauss点的局部坐标系,2020-01-22
	if exist([Full_Pathname,'.fexx_',num2str(POST_Substep)], 'file') ==2
		file_FluidEle_LCS_VectorX_X = fopen([Full_Pathname,'.fexx_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorX_Y = fopen([Full_Pathname,'.fexy_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorX_Z = fopen([Full_Pathname,'.fexz_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorY_X = fopen([Full_Pathname,'.feyx_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorY_Y = fopen([Full_Pathname,'.feyy_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorY_Z = fopen([Full_Pathname,'.feyz_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorZ_X = fopen([Full_Pathname,'.fezx_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorZ_Y = fopen([Full_Pathname,'.fezy_',num2str(POST_Substep)]);
		file_FluidEle_LCS_VectorZ_Z = fopen([Full_Pathname,'.fezz_',num2str(POST_Substep)]);		
		FluidEle_LCS_VectorX_X = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorX_Y = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorX_Z = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorY_X = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorY_Y = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorY_Z = cell(num_Crack(POST_Substep));	
		FluidEle_LCS_VectorZ_X = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorZ_Y = cell(num_Crack(POST_Substep));
		FluidEle_LCS_VectorZ_Z = cell(num_Crack(POST_Substep));		
		for i=1:num_Crack(POST_Substep) 
			FluidEle_LCS_VectorX_X{i} = str2num(fgetl(file_FluidEle_LCS_VectorX_X));
			FluidEle_LCS_VectorX_Y{i} = str2num(fgetl(file_FluidEle_LCS_VectorX_Y));
			FluidEle_LCS_VectorX_Z{i} = str2num(fgetl(file_FluidEle_LCS_VectorX_Z));
			FluidEle_LCS_VectorY_X{i} = str2num(fgetl(file_FluidEle_LCS_VectorY_X));
			FluidEle_LCS_VectorY_Y{i} = str2num(fgetl(file_FluidEle_LCS_VectorY_Y));
			FluidEle_LCS_VectorY_Z{i} = str2num(fgetl(file_FluidEle_LCS_VectorY_Z));
			FluidEle_LCS_VectorZ_X{i} = str2num(fgetl(file_FluidEle_LCS_VectorZ_X));
			FluidEle_LCS_VectorZ_Y{i} = str2num(fgetl(file_FluidEle_LCS_VectorZ_Y));
			FluidEle_LCS_VectorZ_Z{i} = str2num(fgetl(file_FluidEle_LCS_VectorZ_Z));			
		end
		fclose(file_FluidEle_LCS_VectorX_X);
		fclose(file_FluidEle_LCS_VectorX_Y);
		fclose(file_FluidEle_LCS_VectorX_Z);	
		fclose(file_FluidEle_LCS_VectorY_X);
		fclose(file_FluidEle_LCS_VectorY_Y);
		fclose(file_FluidEle_LCS_VectorY_Z);	
		fclose(file_FluidEle_LCS_VectorZ_X);
		fclose(file_FluidEle_LCS_VectorZ_Y);
		fclose(file_FluidEle_LCS_VectorZ_Z);			
	else
	    FluidEle_LCS_VectorX_X=[];
		FluidEle_LCS_VectorX_Y=[];
		FluidEle_LCS_VectorX_Z=[];
	    FluidEle_LCS_VectorY_X=[];
		FluidEle_LCS_VectorY_Y=[];
		FluidEle_LCS_VectorY_Z=[];
	    FluidEle_LCS_VectorZ_X=[];
		FluidEle_LCS_VectorZ_Y=[];
		FluidEle_LCS_VectorZ_Z=[];		
	end
	%读取流体单元Gauss点的局部坐标系,2020-01-22
	if exist([Full_Pathname,'.cfrx_',num2str(POST_Substep)], 'file') ==2
		file_FluidEle_Contact_Force_X = fopen([Full_Pathname,'.cfrx_',num2str(POST_Substep)]);
		file_FluidEle_Contact_Force_Y = fopen([Full_Pathname,'.cfry_',num2str(POST_Substep)]);
		file_FluidEle_Contact_Force_Z = fopen([Full_Pathname,'.cfrz_',num2str(POST_Substep)]);		
		FluidEle_Contact_Force_X = cell(num_Crack(POST_Substep));
		FluidEle_Contact_Force_Y = cell(num_Crack(POST_Substep));
		FluidEle_Contact_Force_Z = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			FluidEle_Contact_Force_X{i} = str2num(fgetl(file_FluidEle_Contact_Force_X));
			FluidEle_Contact_Force_Y{i} = str2num(fgetl(file_FluidEle_Contact_Force_Y));
			FluidEle_Contact_Force_Z{i} = str2num(fgetl(file_FluidEle_Contact_Force_Z));		
		end
		fclose(file_FluidEle_Contact_Force_X);
		fclose(file_FluidEle_Contact_Force_Y);
		fclose(file_FluidEle_Contact_Force_Z);			
	else
	    FluidEle_Contact_Force_X=[];
		FluidEle_Contact_Force_Y=[];
		FluidEle_Contact_Force_Z=[];		
	end	
	%读取裂缝面计算点开度文件
	if exist([Full_Pathname,'.cape_',num2str(POST_Substep)], 'file') ==2
		file_Crack3D_CalP_Aperture = fopen([Full_Pathname,'.cape_',num2str(POST_Substep)]);
		Crack3D_CalP_Aperture = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Crack3D_CalP_Aperture{i} = str2num(fgetl(file_Crack3D_CalP_Aperture));
		end
		fclose(file_Crack3D_CalP_Aperture);
		%获取每个裂缝最大开度值
		for i=1:num_Crack(POST_Substep) 
		    Max_Aperture_of_each_Crack(i) = max(Crack3D_CalP_Aperture{i});
			Min_Aperture_of_each_Crack(i) = min(Crack3D_CalP_Aperture{i});
		    % Max_Aperture_of_each_Crack(i) = max(abs(Crack3D_CalP_Aperture{i}));
			% Min_Aperture_of_each_Crack(i) = min(abs(Crack3D_CalP_Aperture{i}));						
		end
    else
    	Crack3D_CalP_Aperture=[];
	end

	%读取裂缝面计算点开度文件
	if exist([Full_Pathname,'.capf_',num2str(POST_Substep)], 'file') ==2
		file_Crack3D_FluEl_Aperture = fopen([Full_Pathname,'.capf_',num2str(POST_Substep)]);
		Crack3D_FluEl_Aperture = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Crack3D_FluEl_Aperture{i} = str2num(fgetl(file_Crack3D_FluEl_Aperture));
		end
		fclose(file_Crack3D_FluEl_Aperture);
    else
    	Crack3D_FluEl_Aperture=[];
	end
	
	%读取裂缝面离散节点开度文件
	if exist([Full_Pathname,'.cmap_',num2str(POST_Substep)], 'file') ==2
		file_Crack3D_Node_Aperture = fopen([Full_Pathname,'.cmap_',num2str(POST_Substep)]);
		Crack3D_Node_Aperture = cell(num_Crack(POST_Substep));
		for i=1:num_Crack(POST_Substep) 
			Crack3D_Node_Aperture{i} = str2num(fgetl(file_Crack3D_Node_Aperture));
		end
		fclose(file_Crack3D_Node_Aperture);
    else
    	Crack3D_Node_Aperture=[];
	end	
else
    Crack_X = [];   Crack_Y = [];  Crack_Z = [];
	Enriched_Node_Type = [];
	Elem_Type = [];x_cr_tip_nodes=[];y_cr_tip_nodes=[];
	POS = []; Coors_Element_Crack= [];Coors_Vertex= [];
    Coors_Junction= []; Coors_Tip= [];
	Crack_Tip_Type= [];Node_Jun_elem=[];Node_Cross_elem=[];
	Crack_node_X=[];Crack_node_Y=[];Crack_node_Z=[];
	Crack_Ele_1 =[];Crack_Ele_2 =[];Crack_Ele_3 =[];
	Crack_node_local_X=[];Crack_node_local_Y=[];Crack_node_local_Z=[];
	Crack_Node_in_Ele =[];
	Crack3D_Vector_S1_X =[];Crack3D_Vector_S1_Y =[];Crack3D_Vector_S1_Z =[];
	Crack3D_CalP_X =[];Crack3D_CalP_Y=[];Crack3D_CalP_Z=[];
end

% Read enriched nodes of cracks if cracks exist.
if num_Crack(POST_Substep) ~= 0
	disp('    > Reading enriched nodes of cracks....') 
	Post_Enriched_Nodes = load([Full_Pathname,'.ennd_',num2str(POST_Substep)]);
else
	Post_Enriched_Nodes =[];
end

% Read crushed elements if Key_Crush=1
if Key_Crush ==1
    crue_namefile = [Full_Pathname,'.crue_',num2str(POST_Substep)];
	% Check if "crue" file of the POST_Substep substep exist or not.
    if exist(crue_namefile,'file') ==2
        Crushed_element = load(crue_namefile);
	else
	    Crushed_element =[];
	end
else
    Crushed_element =[];
end

% Read nodal average stress file.
if Key_PLOT(3,1)~=0
    disp('    > Reading nodal average stress file....') 
	%读取节点应力
	if Key_Data_Format==1 
        Stress_Matrix = load([Full_Pathname,'.strn_',num2str(POST_Substep)]);
	elseif Key_Data_Format==2  %Binary
		% c_file = fopen([Full_Pathname,'.strn_',num2str(POST_Substep)],'rb');
		% [cc_Stress_Matrix,cc_count]   = fread(c_file,inf,'double');
		% fclose(c_file);
		%%%%转换成Matlab中的数据格式
		% for ccc_i=1:cc_count/4
			% Stress_Matrix(ccc_i,1) = ccc_i;
			% Stress_Matrix(ccc_i,2) = cc_Stress_Matrix(ccc_i*4-3);
			% Stress_Matrix(ccc_i,3) = cc_Stress_Matrix(ccc_i*4-2);
			% Stress_Matrix(ccc_i,4) = cc_Stress_Matrix(ccc_i*4-1);
			% Stress_Matrix(ccc_i,5) = cc_Stress_Matrix(ccc_i*4);
		% end
	end
end

% Get the total force matrix(fx ,fy ,fsum).
FORCE_Matrix = zeros(Num_Node,3);
for i=1:Num_Foc_x
    c_node = Foc_x(i,1);
	FORCE_Matrix(c_node,1) = Foc_x(i,2);
end
for i=1:Num_Foc_y
    c_node = Foc_y(i,1);
	FORCE_Matrix(c_node,2) = Foc_y(i,2);
end
for i=1:Num_Foc_z
    c_node = Foc_z(i,1);
	FORCE_Matrix(c_node,3) = Foc_z(i,2);
end
for i=1:Num_Node
    FORCE_Matrix(i,4) = sqrt(FORCE_Matrix(i,1)^2+FORCE_Matrix(i,2)^2++FORCE_Matrix(i,3)^2);
end

% Plot mesh.
if Key_PLOT(1,1)==1
    Plot_Mesh3D(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

% Plot mesh slice.
if Key_PLOT(1,1)==2 || Key_PLOT(1,1)==3 || Key_PLOT(1,1)==4
    Plot_Mesh3D_Slice(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

% Plot deformation.
if Key_PLOT(2,1)==1
    Plot_Deformation3D(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

% Plot nodal slice stress.
if Key_PLOT(3,1)==1 || Key_PLOT(3,1)==2 || Key_PLOT(3,1)==3
    Plot_Node_Stress_3D_Slice(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

% Plot nodal stress contour,2020-03-11.
if Key_PLOT(3,1)==4
    Plot_Node_Stress_3D(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS,Elem_Type)
end

% Plot nodal slice dispalcement.
if Key_PLOT(4,1)==1 || Key_PLOT(4,1)==2 || Key_PLOT(4,1)==3
    Plot_Node_Disp_3D_Slice(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

% Plot nodal displacement contour,2020-03-11.
if Key_PLOT(4,1)==4
    Plot_Node_Disp_3D(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS,Elem_Type)
end

% Plot crack contour,2020-03-12.
if Key_PLOT(5,1)~=0
    Plot_Crack_Contour_3D(POST_Substep,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
end

disp('    Plot completed.')
disp(' ')

clear DISP
clear Force_X_Matrix
clear Force_Y_Matrix
clear Boundary_X_Matrix
clear Boundary_Y_Matrix
clear Node_Jun_elem
clear POS
clear Crack_X
clear Crack_Y
clear Coors_Element_Crack
clear Coors_Vertex
clear Coors_Junction
clear Coors_Tip
clear Elem_Type