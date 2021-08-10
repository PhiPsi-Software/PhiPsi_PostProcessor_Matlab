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

function Animate_Deformation(Real_num_iteration)
% This function generates the animation of the mesh deformation.

global Node_Coor Elem_Node Outline
global Num_Elem Version Num_Node
global Key_PLOT Full_Pathname Key_Dynamic Foc_x Foc_y
global Size_Font num_Crack Color_Crack
global Color_outline_Udefor Color_Backgro_Defor_1 
global aveg_area_ele Time_Delay delt_time_NewMark Width_Crack
global Output_Freq num_of_Material
global Color_Backgro_Defor_1 Color_Backgro_Defor_2 Color_Backgro_Defor_3 Color_Backgro_Defor_4
global Color_Backgro_Defor_5 Color_Backgro_Defor_6 Color_Backgro_Defor_7
global Color_Backgro_Defor_8 Color_Backgro_Defor_9 Color_Backgro_Defor_10
global Elem_Material Key_Crush Color_Crushed_ele
global Num_Foc_x Num_Foc_y Itera_Num Itera_HF_Time
global Na_Crack_X Na_Crack_Y num_Na_Crack Key_HF_Analysis
global frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl POS_Hl Elem_Type_Hl
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl POS_Incl Elem_Type_Incl
global Color_Inclusion
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Key_Time_String
global Ele_Cross_Point_RABCD
global Key_Data_Format
global Key_Animation Key_Ani_Ave Arc_Crack_Coor
global Title_Font Key_Figure_Control_Widget

disp('    Generating deformation animations......')

%*************************************************	
% Get the max and min value of displacements
%*************************************************	
disp(['    > Calculating the coordinates range of the deformed mesh......']) 
Min_X_Coor_New = zeros(Real_num_iteration,1);
Max_X_Coor_New = zeros(Real_num_iteration,1);
Min_Y_Coor_New = zeros(Real_num_iteration,1);
Max_Y_Coor_New = zeros(Real_num_iteration,1);
for i=1:Real_num_iteration
    if mod(i,Output_Freq)==0
		% DISP = load([Full_Pathname,'.disn_',num2str(Itera_Num(i))]);
		if Key_Data_Format==1 
			DISP   = load([Full_Pathname,'.disn_',num2str(Itera_Num(i))]);
		elseif Key_Data_Format==2  %Binary
			c_file = fopen([Full_Pathname,'.disn_',num2str(Itera_Num(i))],'rb');
			[cc_DISP,cc_count]   = fread(c_file,inf,'double');
			fclose(c_file);
			%转换成Matlab中的数据格式
			% for ccc_i=1:cc_count/2
				% DISP(ccc_i,1) = ccc_i;
				% DISP(ccc_i,2) = cc_DISP(ccc_i*2-1);
				% DISP(ccc_i,3) = cc_DISP(ccc_i*2);
			% end
			ttt  = 1:cc_count/2;
			ttt1 = ttt*2-1;
			ttt2 = ttt*2;
			DISP(ttt,1) = ttt;
			DISP(ttt,2) = cc_DISP(ttt1);
			DISP(ttt,3) = cc_DISP(ttt2);						
		end	
		
		
		scale = Key_PLOT(2,6);
		% Get the new coordinates of all nodes
		New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
		New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);
		clear DISP
		% Get the maximum and minimum value of the new coordinates of all nodes
		Min_X_Coor_New(i) = min(min(New_Node_Coor(:,1)));
		Max_X_Coor_New(i) = max(max(New_Node_Coor(:,1)));
		Min_Y_Coor_New(i) = min(min(New_Node_Coor(:,2)));
		Max_Y_Coor_New(i) = max(max(New_Node_Coor(:,2)));
	end
end
Last_Min_X = min(Min_X_Coor_New);
Last_Max_X = max(Max_X_Coor_New);
Last_Min_Y = min(Min_Y_Coor_New);
Last_Max_Y = max(Max_Y_Coor_New);

%************************************
% Read Boundary file if necessary
%************************************
if Key_PLOT(2,7) == 2 || Key_PLOT(2,7) == 3
	Boundary_X_Matrix = load([Full_Pathname,'.boux']);
	Boundary_Y_Matrix = load([Full_Pathname,'.bouy']);
end

%********************************************
% Get the total force matrix(fx ,fy ,fsum).
%********************************************
FORCE_Matrix = zeros(Num_Node,3);
for i=1:Num_Foc_x
    c_node = Foc_x(i,1);
	FORCE_Matrix(c_node,1) = Foc_x(i,2);
end
for i=1:Num_Foc_y
    c_node = Foc_y(i,1);
	FORCE_Matrix(c_node,2) = Foc_y(i,2);
end
for i=1:Num_Node
    FORCE_Matrix(i,3) = sqrt(FORCE_Matrix(i,1)^2+FORCE_Matrix(i,2)^2);
end

%*****************************************************************	
% Loop through each frame to get the max pressure and quantity 
% of all frame for hydraulic fracturing simulation, i.e.	
% Max_Pressure_for_aimation and Max_Quantity_for_aimation.
%*****************************************************************	
Max_Pressure_for_aimation =0;
Max_Quantity_for_aimation =0;
if exist([Full_Pathname,'.cdpr_',num2str(Itera_Num(1))], 'file') ==2  
	for i=1:Real_num_iteration
		if num_Crack(i) ~= 0
			file_cdpr = fopen([Full_Pathname,'.cdpr_',num2str(Itera_Num(i))]);			
			division_pressure = cell(num_Crack(i));
			max_Pressure=zeros(num_Crack(i));
			for j=1:num_Crack(i) 
				division_pressure{j} = str2num(fgetl(file_cdpr));
				max_Pressure(j) = max(division_pressure{j});
			end
			fclose(file_cdpr);			
			max_max_Pressure(i) = max(max_Pressure);	
		end	
	end
	Max_Pressure_for_aimation = max(max_max_Pressure);
	disp(['    > The maximum pressure of flow of all frames is ',num2str(Max_Pressure_for_aimation),' MPa.']) 

	for i=1:Real_num_iteration
		if num_Crack(i) ~= 0
			file_cdqu = fopen([Full_Pathname,'.cdqu_',num2str(Itera_Num(i))]);
			division_quantity  = cell(num_Crack(i));
			max_Quantity=zeros(num_Crack(i));
			for j=1:num_Crack(i) 
				division_quantity{j} = str2num(fgetl(file_cdqu));
				max_Quantity(j) = max(division_quantity{j});
			end
			fclose(file_cdqu);
			max_max_Quantity(i) = max(max_Quantity);	
		end	
	end
	Max_Quantity_for_aimation = max(max_max_Quantity);
	disp(['    > The maximum quantity of flow of all frames is ',num2str(Max_Quantity_for_aimation),' m^2/s.']) 
end

%*************************************************			
% Loop through each frame to plot and store.
%*************************************************	
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
i_output=0;
for i=1:Real_num_iteration
% for i=5:65
	if mod(i,Output_Freq)==0
	    i_output = i_output + 1;
		disp(['    > Plotting and saving the deformation figure of frame ',num2str(i),'......']) 
		% DISP = load([Full_Pathname,'.disn_',num2str(Itera_Num(i))]);
		if Key_Data_Format==1 
			DISP   = load([Full_Pathname,'.disn_',num2str(Itera_Num(i))]);
		elseif Key_Data_Format==2  %Binary
			c_file = fopen([Full_Pathname,'.disn_',num2str(Itera_Num(i))],'rb');
			[cc_DISP,cc_count]   = fread(c_file,inf,'double');
			fclose(c_file);
			%转换成Matlab中的数据格式
			% for ccc_i=1:cc_count/2
				% DISP(ccc_i,1) = ccc_i;
				% DISP(ccc_i,2) = cc_DISP(ccc_i*2-1);
				% DISP(ccc_i,3) = cc_DISP(ccc_i*2);
			% end
			ttt  = 1:cc_count/2;
			ttt1 = ttt*2-1;
			ttt2 = ttt*2;
			DISP(ttt,1) = ttt;
			DISP(ttt,2) = cc_DISP(ttt1);
			DISP(ttt,3) = cc_DISP(ttt2);						
		end	
		
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read crushed elements if Key_Crush=1
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_Crush ==1
			crue_namefile = [Full_Pathname,'.crue_',num2str(Itera_Num(i))];
			% Check if "crue" file of the i substep exist or not.
			if exist(crue_namefile,'file') ==2
				Crushed_element = load(crue_namefile);
			else
				Crushed_element =[];
			end
		else
			Crushed_element =[];
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read kiel file.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		disp('    > Reading kiel (killed elements) file....') 
		num_Killed_Ele=0;
		if exist([Full_Pathname,'.kiel_',num2str(Itera_Num(i))], 'file') ==2  
			Killed_Elements = load([Full_Pathname,'.kiel_',num2str(Itera_Num(i))]);
			num_Killed_Ele  = size(Killed_Elements,1);  %被杀死单元的数目
		end

		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read coordinates and all other information of cracks if cracks exist.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% if num_Crack(Real_num_iteration) ~=0
		if num_Crack(i) ~= 0
			file_Crack_X = fopen([Full_Pathname,'.crax_',num2str(Itera_Num(i))]);
			file_Crack_Y = fopen([Full_Pathname,'.cray_',num2str(Itera_Num(i))]);
			Crack_X = cell(num_Crack(i));
			Crack_Y = cell(num_Crack(i));
			for j=1:num_Crack(i) 
				Crack_X{j} = str2num(fgetl(file_Crack_X));
				Crack_Y{j} = str2num(fgetl(file_Crack_Y));
			end
			fclose(file_Crack_X);
			fclose(file_Crack_Y);
			if Key_Data_Format==1 
				Enriched_Node_Type = load([Full_Pathname,'.ennd_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.ennd_',num2str(Itera_Num(i))],'rb');
				[cc_Enriched_Node_Type,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Enriched_Node_Type = (reshape(cc_Enriched_Node_Type,num_Crack(Itera_Num(i)),Num_Node))';
			end				
			
			if Key_Data_Format==1 
				POS = load([Full_Pathname,'.posi_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.posi_',num2str(Itera_Num(i))],'rb');
				[cc_POS,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				POS = (reshape(cc_POS,num_Crack(Itera_Num(i)),Num_Node))';
			end		
			
			
			if Key_Data_Format==1 
				Elem_Type = load([Full_Pathname,'.elty_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.elty_',num2str(Itera_Num(i))],'rb');
				[cc_Elem_Type,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Elem_Type = (reshape(cc_Elem_Type,num_Crack(Itera_Num(i)),Num_Elem))';
			end
			
			if Key_Data_Format==1 
				tem_Coors_Element_Crack = load([Full_Pathname,'.celc_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.celc_',num2str(Itera_Num(i))],'rb');
				[cc_tem_Coors_Element_Crack,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				for ccc_i=1:cc_count/4
					tem_Coors_Element_Crack(ccc_i,1) = cc_tem_Coors_Element_Crack(ccc_i*4-3);
					tem_Coors_Element_Crack(ccc_i,2) = cc_tem_Coors_Element_Crack(ccc_i*4-2);
					tem_Coors_Element_Crack(ccc_i,3) = cc_tem_Coors_Element_Crack(ccc_i*4-1);
					tem_Coors_Element_Crack(ccc_i,4) = cc_tem_Coors_Element_Crack(ccc_i*4);
				end
			end
			
			
			num_ele = size(tem_Coors_Element_Crack,1)/num_Crack(i);
			for i_Cr = 1:num_Crack(i)
				for i_Ele = 1:num_ele
					Coors_Element_Crack(i_Ele,i_Cr,1:4) = tem_Coors_Element_Crack((i_Cr-1)*num_ele + i_Ele,1:4);
				end
			end
			
			Node_Cross_elem = [];
			
			if Key_Data_Format==1 
				Node_Jun_elem = load([Full_Pathname,'.njel_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.njel_',num2str(Itera_Num(i))],'rb');
				[cc_Node_Jun_elem,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Node_Jun_elem = (reshape(cc_Node_Jun_elem,num_Crack(Itera_Num(i)),Num_Node))';
			end	
			if Key_Data_Format==1 
				Node_Jun_Hole = load([Full_Pathname,'.njhl_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.njhl_',num2str(Itera_Num(i))],'rb');
				[cc_Node_Jun_Hole,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Node_Jun_Hole = (reshape(cc_Node_Jun_Hole,num_Crack(Itera_Num(i)),Num_Node))';
			end	
			
			if Key_Data_Format==1 
				Coors_Vertex        = load([Full_Pathname,'.celv_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.celv_',num2str(Itera_Num(i))],'rb');
				[cc_Coors_Vertex,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Coors_Vertex = (reshape(cc_Coors_Vertex,2,num_ele))';
			end

			if Key_Data_Format==1 
				tem_Coors_Junction      = load([Full_Pathname,'.celj_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.celj_',num2str(Itera_Num(i))],'rb');
				[cc_tem_Coors_Junction,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				for ccc_i=1:cc_count/4
					tem_Coors_Junction(ccc_i,1) = cc_tem_Coors_Junction(ccc_i*4-3);
					tem_Coors_Junction(ccc_i,2) = cc_tem_Coors_Junction(ccc_i*4-2);
					tem_Coors_Junction(ccc_i,3) = cc_tem_Coors_Junction(ccc_i*4-1);
					tem_Coors_Junction(ccc_i,4) = cc_tem_Coors_Junction(ccc_i*4);
				end
			end
			
			
			num_ele = size(tem_Coors_Element_Crack,1)/num_Crack(i);
			for i_Cr = 1:num_Crack(i)
				for i_Ele = 1:num_ele
					Coors_Junction(i_Ele,i_Cr,1:4) = tem_Coors_Junction((i_Cr-1)*num_ele + i_Ele,1:4);
				end
			end
			
			disp('    > Reading celt file....');
			if Key_Data_Format==1 
				Coors_Tip           = load([Full_Pathname,'.celt_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.celt_',num2str(Itera_Num(i))],'rb');
				[cc_Coors_Tip,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Coors_Tip = (reshape(cc_Coors_Tip,2,num_ele))';
			end				
			
			Crack_Tip_Type      = load([Full_Pathname,'.ctty_',num2str(Itera_Num(i))]);
			% end
		else
		    % 1111
			Crack_X = [];   Crack_Y = [];
			Enriched_Node_Type = [];
			Crack_Tip_Type = [];
			Elem_Type = [];
			POS = []; Coors_Element_Crack= [];Coors_Vertex= [];Node_Jun_elem=[];Node_Cross_elem=[];
			Coors_Junction= []; Coors_Tip= []; Elem_Type= [];Node_Jun_Hole=[];
		end
		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read coordinates of arc (Arc_Crack_Coor(,,1:11)) crack if arc crack exist, 2017-07-26.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if num_Crack(Real_num_iteration) ~= 0
		    if num_Crack(i) ~= 0
				if exist([Full_Pathname,'.arcc_',num2str(Itera_Num(i))], 'file') ==2  
					disp('    > Reading coordinates of arc cracks....') 
					Yes_Arc_Crack = 1;
					file_Arc = fopen([Full_Pathname,'.arcc_',num2str(Itera_Num(i))]);
					for iii_crack=1:num_Crack(Itera_Num(i))
						nPt = size(Crack_X{iii_crack},2);
						for jjj=1:11
							Arc_Crack_Coor(iii_crack,1:nPt-1,jjj)= str2num(fgetl(file_Arc));
						end
					end
					fclose(file_Arc);
				else 
					Yes_Arc_Crack  = 0;
					Arc_Crack_Coor(1:1000,1:1000,1:11) = 0.0;
				end
			end
		else 
			% Arc_Crack_Coor(1:100,1:500,1:11) = 0.0;
			Yes_Arc_Crack  = 0;
			Arc_Crack_Coor(1:1000,1:1000,1:11) = 0.0;
		end
		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read coordinates and other info of holes if holes exist.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if num_Hole ~= 0
			disp('    > Reading coordinates of hole....') 
			file_Hole = fopen([Full_Pathname,'.hlcr']);
			for iHole=1:num_Hole
				Hole_Coor(iHole,1:3)= str2DOUBLE(fgetl(file_Hole));
			end
			fclose(file_Hole);
			disp('    > Reading ennh file of hole....') 
			Enriched_Node_Type_Hl = load([Full_Pathname,'.ennh_',num2str(Itera_Num(i))]);
			disp('    > Reading posh file of hole....') 
			POS_Hl = load([Full_Pathname,'.posh_',num2str(Itera_Num(i))]);
			
			disp('    > Reading elth file of hole....');
			Elem_Type_Hl = load([Full_Pathname,'.elth_',num2str(Itera_Num(i))]);
		end		
		
		% 读取圆形夹杂的坐标
		if num_Circ_Inclusion ~= 0
			disp('    > Reading coordinates of circle inclusions....') 
			file_Circ_Inclusion = fopen([Full_Pathname,'.jzcr']);
			for iiii=1:num_Circ_Inclusion
				Circ_Inclu_Coor(iiii,1:3)= str2DOUBLE(fgetl(file_Circ_Inclusion));
			end
			fclose(file_Circ_Inclusion);
			disp('    > Reading ennh file of circle inclusions....') 
			Enriched_Node_Type_Incl = load([Full_Pathname,'.ennj_',num2str(Itera_Num(i))]);
			disp('    > Reading posh file of circle inclusions....') 
			POS_Incl = load([Full_Pathname,'.posj_',num2str(Itera_Num(i))]);
			disp('    > Reading elth file of circle inclusions....');
			Elem_Type_Incl = load([Full_Pathname,'.eltj_',num2str(Itera_Num(i))]);
		end

		% 读取多边形夹杂的坐标,added on 2016-10-04.
		if num_Poly_Inclusion ~= 0
			disp('    > Reading coordinates of polygon inclusions....') 
			file_Poly_Incl_Coor_x = fopen([Full_Pathname,'.jzpx']);
			file_Poly_Incl_Coor_y = fopen([Full_Pathname,'.jzpy']);
			Poly_Incl_Coor_x = cell(num_Poly_Inclusion);
			Poly_Incl_Coor_y = cell(num_Poly_Inclusion);
			for iiii=1:num_Poly_Inclusion
				Poly_Incl_Coor_x{iiii} = str2DOUBLE(fgetl(file_Poly_Incl_Coor_x));
				Poly_Incl_Coor_y{iiii} = str2DOUBLE(fgetl(file_Poly_Incl_Coor_y));
			end
			fclose(file_Poly_Incl_Coor_x);
			fclose(file_Poly_Incl_Coor_y);
			
			disp('    > Reading ennh file of polygon inclusions....') 
			Enriched_Node_Type_Incl = load([Full_Pathname,'.ennj_',num2str(Itera_Num(i))]);
			disp('    > Reading posh file of polygon inclusions....') 
			POS_Incl = load([Full_Pathname,'.posj_',num2str(Itera_Num(i))]);
			disp('    > Reading elth file of polygon inclusions....');
			Elem_Type_Incl = load([Full_Pathname,'.eltj_',num2str(Itera_Num(i))]);
		end
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Read enriched nodes of cracks if cracks exist.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if num_Crack(i) ~= 0
			disp(['      ----- Reading enriched nodes of cracks......'])		
			if Key_Data_Format==1 
				Post_Enriched_Nodes = load([Full_Pathname,'.ennd_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.ennd_',num2str(Itera_Num(i))],'rb');
				[cc_Post_Enriched_Nodes,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Post_Enriched_Nodes = (reshape(cc_Post_Enriched_Nodes,num_Crack(Itera_Num(i)),Num_Node))';
			end	
				
		else
			Post_Enriched_Nodes =[];
		end
		scale = Key_PLOT(2,6);
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Get the new coordinates of all nodes
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		disp(['      ----- Get the new coordinates of all node......'])		
		New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
		New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);
		xi = zeros(4,Num_Elem); yi = xi;
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];   % Nodes for current element
			if Elem_Material(iElem)==1
				xi_1(:,iElem) = New_Node_Coor(NN',1);       % Initial x-coordinates of nodes
				yi_1(:,iElem) = New_Node_Coor(NN',2);       % Initial y-coordinates of nodes		
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
		end
		% New figure.
		Tools_New_Figure
		hold on;
		title('Deformation','FontName',Title_Font,'FontSize',Size_Font)
		axis off; 
		axis equal;
		delta = sqrt(aveg_area_ele);
		axis([Last_Min_X-delta  Last_Max_X+delta  Last_Min_Y-delta  Last_Max_Y+delta]); 
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	  
		% Plot deformed mesh.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		disp(['      ----- Plot deformed mesh......'])		
		patch(xi_1,yi_1,Color_Backgro_Defor_1) 
		patch(xi_2,yi_2,Color_Backgro_Defor_2)    
		patch(xi_3,yi_3,Color_Backgro_Defor_3)   
		patch(xi_4,yi_4,Color_Backgro_Defor_4)  
		patch(xi_5,yi_5,Color_Backgro_Defor_5)  
		patch(xi_6,yi_6,Color_Backgro_Defor_6)  
		patch(xi_7,yi_7,Color_Backgro_Defor_7)  
		patch(xi_8,yi_8,Color_Backgro_Defor_8)  
		patch(xi_9,yi_9,Color_Backgro_Defor_9)  
		patch(xi_10,yi_10,Color_Backgro_Defor_10)  
		%杀死的单元用白色
		if num_Killed_Ele>0
		    %num_Killed_Ele
			for iElem = 1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
				if ismember(iElem,Killed_Elements)
					xi_killed(1:4) = New_Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
					yi_killed(1:4) = New_Node_Coor(NN',2);  
					patch(xi_killed,yi_killed,'white') 
					% 1111
				end
			end 
		end		
		
		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot crushed elements.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if isempty(Crushed_element)==0
			xi_crushed =[];    yi_crushed =[];
			size(Crushed_element);
			for iElem = 1:Num_Elem
				if Crushed_element(iElem)==1
					NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
						  Elem_Node(iElem,3) Elem_Node(iElem,4)];          % Nodes for current element
					xi_crushed(:,iElem) = New_Node_Coor(NN',1);            % Deformed x-coordinates of nodes
					yi_crushed(:,iElem) = New_Node_Coor(NN',2);            % Deformed y-coordinates of nodes 
				end		
			end
			patch(xi_crushed,yi_crushed,Color_Crushed_ele) 
		end		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot undeformed outline
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,8)==1
			line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
				 [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1.5,'Color',Color_outline_Udefor)     			 
		end   
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		%       单元接触状态
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,13)==1 
			Yes_Contact1=0;  %裂纹面闭合
			Yes_Contact2=0;  %裂纹面由支撑剂支撑
			for iElem = 1:Num_Elem
				xi_Elcs1(1:4,iElem) = [0;0;0;0];
				yi_Elcs1(1:4,iElem) = [0;0;0;0];
				xi_Elcs2(1:4,iElem) = [0;0;0;0];
				yi_Elcs2(1:4,iElem) = [0;0;0;0];
			end
			%如果接触状态文件存在则:
			if exist([Full_Pathname,'.elcs_',num2str(Itera_Num(i))], 'file') ==2  
				disp(['      ----- Plotting contact status of elements...'])
				ELCS= load([Full_Pathname,'.elcs_',num2str(Itera_Num(i))]);
				for iElem = 1:Num_Elem
				
					NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
						  Elem_Node(iElem,3) Elem_Node(iElem,4)];         % Nodes for current element
					if ELCS(iElem) == 1    %裂纹面闭合
					% size(New_Node_Coor(NN',1))
						xi_Elcs1(:,iElem) = New_Node_Coor(NN',1);                                     
						yi_Elcs1(:,iElem) = New_Node_Coor(NN',2); 	
						Yes_Contact1=1;
					end
					if ELCS(iElem) == 2    %裂纹面滑动
						xi_Elcs2(:,iElem) = New_Node_Coor(NN',1);                                     
						yi_Elcs2(:,iElem) = New_Node_Coor(NN',2); 	
						Yes_Contact2=1;
					end
				end
				if Yes_Contact1==1
					patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])          %用粉红色标记该单�?
				end
				if Yes_Contact2==1
					patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])    %用紫色标记该单元
				end
			end
			%如果粘聚裂缝粘聚单元状态文件存在则(2017-04-24):
			if exist([Full_Pathname,'.elco_',num2str(Itera_Num(i))], 'file') ==2  
				disp(['      ----- Plotting cohesive crack status of elements...'])
				ELCS= load([Full_Pathname,'.elco_',num2str(Itera_Num(i))]);
				for iElem = 1:Num_Elem
					NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
						  Elem_Node(iElem,3) Elem_Node(iElem,4)];         % Nodes for current element
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
					patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])               %用粉红色标记该单�?
				end
				if Yes_Contact2==1
					patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])         %用紫色标记该单元
				end
			end
		end
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot boundary conditions by darwing triangles
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,7) == 2 || Key_PLOT(2,7) == 3
			% W = Max_X_Coor_New - Min_X_Coor_New;
			% H = Max_Y_Coor_New - Min_Y_Coor_New;
			% Length of the edge of the boundary triangle
			% l_edge = max(W,H)/20.0; 
			l_edge = 0.6*sqrt(aveg_area_ele);
			% Loop through boundary_x nodes
			for j = 1:size(Boundary_X_Matrix,1)
				c_node   = Boundary_X_Matrix(j);
				x_c_node = 	New_Node_Coor(c_node,1);
				y_c_node =  New_Node_Coor(c_node,2);
				% The x location of the boundary triangle	
				tri_X = [x_c_node   x_c_node-sqrt(3)/2*l_edge  x_c_node-sqrt(3)/2*l_edge  x_c_node]; 
				% The y location of the boundary triangle
				tri_Y = [y_c_node   y_c_node+1/2*l_edge        y_c_node-1/2*l_edge        y_c_node]; 
				% Plot the boundary triangle
				plot(tri_X,tri_Y,'color','black')                                                                  
			end
			% Loop through boundary_y nodes
			for j = 1:size(Boundary_Y_Matrix,1)
				c_node   = Boundary_Y_Matrix(j);
				x_c_node = 	New_Node_Coor(c_node,1);
				y_c_node =  New_Node_Coor(c_node,2);
				% The x location of the boundary triangle
				tri_X = [x_c_node   x_c_node-1/2*l_edge         x_c_node+1/2*l_edge        x_c_node]; 
				% The y location of the boundary triangle
				tri_Y = [y_c_node   y_c_node-sqrt(3)/2*l_edge   y_c_node-sqrt(3)/2*l_edge  y_c_node]; 
				% Plot the boundary triangle
				plot(tri_X,tri_Y,'color','black')                                                
			end
		end
		%<<<<<<<<<<<<<<<<<<<
		% Plot forces.
		%<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,7) == 1 || Key_PLOT(2,7) == 3
			disp(['      ----- Plotting forces of nodes......'])
			Max_x_Force = max(abs(FORCE_Matrix(:,1)));
			Max_y_Force = max(abs(FORCE_Matrix(:,2)));
			Max_Force   = max(Max_x_Force,Max_y_Force);
			W = Max_X_Coor_New(i) - Min_X_Coor_New(i);
			H = Max_Y_Coor_New(i) - Min_Y_Coor_New(i);
			% length of force arrow
			% REMOVE:length_arrow = sqrt(max_area_ele);
			length_arrow = max(W,H)/15.0;    
			% Loop through each node.
			for i_Node = 1:Num_Node
				if FORCE_Matrix(i_Node,3) ~=0           % If the nodes has force load, then:
					c_force_x   = FORCE_Matrix(i_Node,1);
					c_force_y   = FORCE_Matrix(i_Node,2);
					delta_L_x = c_force_x*length_arrow/Max_Force;
					delta_L_y = c_force_y*length_arrow/Max_Force;
					
					StartPoint = [New_Node_Coor(i_Node,1)-delta_L_x   New_Node_Coor(i_Node,2)-delta_L_y];
					EndPoint   = [New_Node_Coor(i_Node,1)             New_Node_Coor(i_Node,2)          ];
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

		%<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot Gauss points.
		%<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,4) == 1
			disp(['      ----- Plotting Gauss points......'])
			% Read gauss 坐标
			if Key_Data_Format==1 
			    Gauss_Coor = load([Full_Pathname,'.gcor_',num2str(i)]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.gcor_',num2str(i)],'rb');
				[cc_Read_Gauss_Coor,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				for ccc_i=1:cc_count/2
					Gauss_Coor(ccc_i,1) = ccc_i;
					Gauss_Coor(ccc_i,2) = cc_Read_Gauss_Coor(ccc_i*2-1);
					Gauss_Coor(ccc_i,3) = cc_Read_Gauss_Coor(ccc_i*2);
				end
			end			
			
			% 读取Gauss点位移.
			if Key_Data_Format==1 
			    Gauss_Disp = load([Full_Pathname,'.disg_',num2str(i)]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.disg_',num2str(i)],'rb');
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
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot cracks line if necessary
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,5) == 1
			if num_Crack(i)~=0
				for i_Crack = 1:num_Crack(i)
					nPt = size(Crack_X{i_Crack},2);
					for i_Seg = 1:nPt-1
					    x = [Crack_X{i_Crack}(i_Seg) Crack_X{i_Crack}(i_Seg+1)];
					    y = [Crack_Y{i_Crack}(i_Seg) Crack_Y{i_Crack}(i_Seg+1)];
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
		                %----------------------------
						%如果是弧形裂缝,2017-07-26
						%----------------------------
						%Arc_Crack_Coor:x,y,r,Radian_Start,Radian_End,Radian,Point_Start_x,Point_Start_y,Point_End_x,Point_End_y
						if abs(sum(Arc_Crack_Coor(i_Crack,i_Seg,1:11)))>=1.0e-10 
							c_R=Arc_Crack_Coor(i_Crack,i_Seg,4);
							c_Direcction  =Arc_Crack_Coor(i_Crack,i_Seg,3);
							c_Radian_Start=Arc_Crack_Coor(i_Crack,i_Seg,5)*pi/180;
							c_Radian_End  =Arc_Crack_Coor(i_Crack,i_Seg,6)*pi/180;
							%**********************
							%   如果是逆时针圆弧
							%**********************
							if c_Direcction >0.5
								%若结束角度大于起始角,则直接绘制圆弧
								if c_Radian_End>=c_Radian_Start 
									c_alpha=c_Radian_Start:pi/100:c_Radian_End;
									c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
								%若结束角度小于起始角,则分成两部分绘制圆弧
								else
									%第1部分:Radian_Start到360
									c_alpha1=c_Radian_Start:pi/100:2*pi;
									c_x=c_R*cos(c_alpha1)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y=c_R*sin(c_alpha1)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
									c_x_2=c_R*cos(c_alpha_2)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y_2=c_R*sin(c_alpha_2)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
								%若结束角度大于起始角,则分成两部分绘制圆弧
								if c_Radian_End>=c_Radian_Start 
									%第1部分:Radian_End到360
									c_alpha1=c_Radian_End:pi/100:2*pi;
									c_x=c_R*cos(c_alpha1)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y=c_R*sin(c_alpha1)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
									c_alpha_2=0.0:pi/100:c_Radian_Start;
									c_x_2=c_R*cos(c_alpha_2)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y_2=c_R*sin(c_alpha_2)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
								%若结束角度小于起始角,则直接绘制圆弧
								else
									c_alpha=c_Radian_End:pi/100:c_Radian_Start;
									c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i_Crack,i_Seg,1);
									c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i_Crack,i_Seg,2);
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
						elseif abs(sum(Arc_Crack_Coor(i_Crack,i_Seg,1:11)))<1.0e-10 					
							plot(last_x,last_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   
						end
					end
				end	
			end
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot holes.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		disp(['      ----- Plotting hole...'])
		if num_Hole ~=0
			for iHole = 1:num_Hole
				Coor_x  = Hole_Coor(iHole,1);
				Coor_y  = Hole_Coor(iHole,2);
				c_R  = Hole_Coor(iHole,3);
				num_fineness = 100;
				for j_P = 1:num_fineness+1
					alpha = 2*pi/num_fineness*(j_P-1);
					c_x(j_P) = Coor_x + c_R*cos(alpha);
					c_y(j_P) = Coor_y + c_R*sin(alpha);
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(j_P),c_y(j_P));
                    %Sometimes the element does not exist, for example, only part of the hole locates inside the model	
			        if c_Elem_Num~=0
						[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(j_P),c_y(j_P));
					    [c_dis_x(j_P),c_dis_y(j_P)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
					end
				end
				x_new = c_x + c_dis_x*scale;
				y_new = c_y + c_dis_y*scale;
				% plot(x_new,y_new,'-')
				patch(x_new,y_new,'white','edgecolor','black','LineWidth',0.1)	
			end	
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制圆形夹杂
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		disp(['      ----- Plotting circle inclusions...'])
		if num_Circ_Inclusion ~=0
			for i_Inclusion = 1:num_Circ_Inclusion
				Coor_x  = Circ_Inclu_Coor(i_Inclusion,1);
				Coor_y  = Circ_Inclu_Coor(i_Inclusion,2);
				c_R  = Circ_Inclu_Coor(i_Inclusion,3);
				num_fineness = 100;
				for j_P = 1:num_fineness+1
					alpha = 2*pi/num_fineness*(j_P-1);
					c_x(j_P) = Coor_x + c_R*cos(alpha);
					c_y(j_P) = Coor_y + c_R*sin(alpha);
					[Kesi,Yita] = Cal_KesiYita_by_Coors(c_x(j_P),c_y(j_P));
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(c_x(j_P),c_y(j_P));
					[c_dis_x(j_P),c_dis_y(j_P)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				end
				x_new = c_x + c_dis_x*scale;
				y_new = c_y + c_dis_y*scale;
				% plot(x_new,y_new,'-')
				patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'	
			end	
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制圆形夹杂,2016-10-19
		%<<<<<<<<<<<<<<<<<<<<<<<<<
		disp(['      ----- Plotting circle inclusion...'])
		if num_Circ_Inclusion ~=0
			for i_Inclusion = 1:num_Circ_Inclusion
				Coor_x  = Circ_Inclu_Coor(i_Inclusion,1);
				Coor_y  = Circ_Inclu_Coor(i_Inclusion,2);
				c_R  = Circ_Inclu_Coor(i_Inclusion,3);
				num_fineness = 100;
				for cc_j = 1:num_fineness+1
					alpha = 2*pi/num_fineness*(cc_j-1);
					cc_x(cc_j) = Coor_x + c_R*cos(alpha);
					cc_y(cc_j) = Coor_y + c_R*sin(alpha);
					[Kesi,Yita] = Cal_KesiYita_by_Coors(cc_x(cc_j),cc_y(cc_j));
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(cc_x(cc_j),cc_y(cc_j));
					[c_dis_x(cc_j),c_dis_y(cc_j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				end
				x_new = cc_x + c_dis_x*scale;
				y_new = cc_y + c_dis_y*scale;
				% 绘制
				plot(x_new,y_new,'-','color','black')
				% 填充
				% patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
			end	
		end
			
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制多边形夹杂,2016-10-19
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if num_Poly_Inclusion ~=0
			disp(['      ----- Plotting poly inclusion...'])
			for iii = 1:num_Poly_Inclusion
				nEdge = size(Poly_Incl_Coor_x{iii},2);
				Num_Diversion = 5;    %多边形的每条边拆分成5份
				cc_k = 0;
				for iEdge = 1:nEdge
					%获取边线的起点和终点
					if iEdge==nEdge
						Line_Edge(1,1) = Poly_Incl_Coor_x{iii}(iEdge);   %边线的起点
						Line_Edge(1,2) = Poly_Incl_Coor_y{iii}(iEdge);
						Line_Edge(2,1) = Poly_Incl_Coor_x{iii}(1);       %边线的终点
						Line_Edge(2,2) = Poly_Incl_Coor_y{iii}(1);
					else
						Line_Edge(1,1) = Poly_Incl_Coor_x{iii}(iEdge);   %边线的起点
						Line_Edge(1,2) = Poly_Incl_Coor_y{iii}(iEdge);
						Line_Edge(2,1) = Poly_Incl_Coor_x{iii}(iEdge+1); %边线的终点
						Line_Edge(2,2) = Poly_Incl_Coor_y{iii}(iEdge+1);
					end
					% 等分点.
					a_x = Line_Edge(1,1);
					a_y = Line_Edge(1,2);
					b_x = Line_Edge(2,1);
					b_y = Line_Edge(2,2);
					%计算边线起点的位置
					cc_k =cc_k+1;
					cc_x(cc_k) = a_x;
					cc_y(cc_k) = a_y;
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(cc_x(cc_k),cc_y(cc_k));
			        if c_Elem_Num~=0
					    [Kesi,Yita] = Cal_KesiYita_by_Coors(cc_x(cc_k),cc_y(cc_k));
					    [cc_dis_x(cc_k),cc_dis_y(cc_k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
					end
					%计算等分点的位移
					for cc_j =  1:Num_Diversion-1
						cc_k=cc_k+1;
						cc_x(cc_k) = (cc_j*b_x+(Num_Diversion-cc_j)*a_x)/Num_Diversion;
						cc_y(cc_k) = (cc_j*b_y+(Num_Diversion-cc_j)*a_y)/Num_Diversion;
						
						[c_Elem_Num] = Cal_Ele_Num_by_Coors(cc_x(cc_k),cc_y(cc_k));
						if c_Elem_Num~=0
					    	[Kesi,Yita] = Cal_KesiYita_by_Coors(cc_x(cc_k),cc_y(cc_k));
						    [cc_dis_x(cc_k),cc_dis_y(cc_k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
						end
					end												      
				end
				x_new = cc_x + cc_dis_x*scale;
				y_new = cc_y + cc_dis_y*scale;
				................
				% option1:填充
				%................
				patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
				................
				% option2:绘制边线
				%................
				% x_new(cc_k+1) = x_new(1);
				% y_new(cc_k+1) = y_new(1);
				% plot(x_new,y_new,'-','color','black')
			end	
		end	
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制天然裂缝.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Calculate and plot shaped cracks if necessary.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,5) == 2 & num_Crack(i)~=0
			disp(['      ----- Plotting shaped cracks......'])
			[Shaped_Crack_Points] = Cal_Shaped_Cracks(Crack_X,Crack_Y,i,Itera_Num(i),num_Crack,Crack_Tip_Type,POS,...
									Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
									Coors_Vertex,Coors_Junction,Coors_Tip,DISP,scale);
			Plot_Shaped_Cracks(Shaped_Crack_Points);
		end
		%<<<<<<<<<<<<<<<<<<<<<<
		%绘制支撑剂的圆球
		%<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(2,10)==1 
			%如果支撑剂直径和坐标文件存在
			if exist([Full_Pathname,'.epcr_',num2str(Itera_Num(i))], 'file') ==2   
                disp(['      ----- Plotting proppant......'])				
				Plotting shaped cracks
				%读取文件		
				c_D_Coors=load([Full_Pathname,'.epcr_',num2str(Itera_Num(i))]);
				num_Proppant = size(c_D_Coors,1);
				for tt_i=1:num_Proppant
					%绘制实心圆球
					alpha=0:pi/20:2*pi;
					c_Elem_Num = c_D_Coors(tt_i,1);
					R = c_D_Coors(tt_i,2)/2*Key_PLOT(2,6);
					old_Coor_x = c_D_Coors(tt_i,3);
					old_Coor_y = c_D_Coors(tt_i,4);
					omega = c_D_Coors(tt_i,5);  %对应裂纹片段的倾角角
					l_elem= sqrt(aveg_area_ele);
					offset_delta = 0.001*l_elem;
					%原支撑剂在点的上裂纹偏置微小量之后的坐标
					old_Coor_x_Up  = old_Coor_x - offset_delta*sin(omega);
					old_Coor_y_Up  = old_Coor_y + offset_delta*cos(omega);
					%原支撑剂在点的上裂纹偏置微小量之后下偏置
					old_Coor_x_Low = old_Coor_x + offset_delta*sin(omega);
					old_Coor_y_Low = old_Coor_y - offset_delta*cos(omega);
					%计算支撑剂中心坐标变形后的坐标
					[Kesi_Up,Yita_Up] = Cal_KesiYita_by_Coors(old_Coor_x_Up,old_Coor_y_Up);
					[dis_x_Up,dis_y_Up] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi_Up,Yita_Up...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
					[Kesi_Low,Yita_Low] = Cal_KesiYita_by_Coors(old_Coor_x_Low,old_Coor_y_Low);
					[dis_x_Low,dis_y_Low] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi_Low,Yita_Low...
																	 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																	  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
					%真实的支撑剂在点的位置											  
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
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot text on the left or the bottom of the figure.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		plot_string_MatMES = ['PhiPsi  ',Version];
		plot_string_Frame  = ['Frame ',num2str(i),' / ',num2str(Real_num_iteration)];
		plot_string_Time   = ['Time: ',num2str(i*delt_time_NewMark*1000),' ms'];
		% if Key_HF_Analysis==1
			% plot_string_Time_HF   = ['Time: ',num2str(Itera_HF_Time(i)/60.0,'%0.4f'),' mins'];
		% end
		if  exist([Full_Pathname,'.hftm'], 'file') ==2 
			if Key_Time_String==1
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output),'%0.4f'),' s'];
			elseif Key_Time_String==2
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output)/60.0,'%0.4f'),' mins'];
			elseif Key_Time_String==3
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output)/60.0/60.0,'%0.4f'),' hours'];
			elseif Key_Time_String==4
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output)/60.0/60.0/24,'%0.4f'),' days'];			
			elseif Key_Time_String==5
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output)/60.0/60.0/24/30.41666,'%0.4f'),' months'];					
			elseif Key_Time_String==6
				plot_string_Time   = ['Time: ',num2str(Itera_HF_Time(i_output)/60.0/60.0/24/30.41666/12.0,'%0.4f'),' years'];	
			end	
		end
		%plot_string_Time_HF   = ['Time: ms'];
		plot_string_Scale  = ['Scale factor: ',num2str(scale)];
		range_W = abs(Last_Max_X-Last_Min_X);
		range_H = abs(Last_Max_Y-Last_Min_Y);
		if range_H >= 0.75*range_W      % Left
			loc_x = -range_H/2+ Last_Min_X;
			loc_y =  Last_Max_Y-range_H*0.05;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_y =  Last_Max_Y-range_H*0.15;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_y =  Last_Max_Y-range_H*0.25;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_y =  Last_Max_Y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			% if Key_HF_Analysis==1
				% loc_y =  Last_Max_Y-range_H*0.35;
				% text(loc_x,loc_y,plot_string_Time_HF, 'color','black');
			% end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_y =  Last_Max_Y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
		else                            % Bottom
			loc_y =  Last_Min_Y-range_H*0.05;
			loc_x =  Last_Min_X;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_x =  Last_Min_X + range_W*0.25;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_x =  Last_Min_X + range_W*0.45;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_x =  Last_Min_X + range_W*0.60;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			% if Key_HF_Analysis==1
				% loc_x =  Last_Min_X + range_W*0.60;
				% text(loc_x,loc_y,plot_string_Time_HF, 'color','black');
			% end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_x =  Last_Min_X + range_W*0.60;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Save the current figure
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		deformation(i_output) = getframe(gcf);       
		im=frame2im(deformation(i_output));         
		[II,map]=rgb2ind(im,256);
		kkkk=i_output-0;
		str1='_deformation';
		str2=Full_Pathname;
		FileName =[str2,str1,'.gif'];
		if kkkk==1;
			imwrite(II,map,FileName,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
		else
			imwrite(II,map,FileName,'gif','WriteMode','append','DelayTime',Time_Delay);
		end
		
		close
    end
end