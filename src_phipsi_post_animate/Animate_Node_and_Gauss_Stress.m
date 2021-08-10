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

function Animate_Node_and_Gauss_Stress(Real_num_iteration)
% This function generates the animation of nodal stress contours.

global Node_Coor Elem_Node Outline Inline
global Num_Node Num_Elem Version
global Key_PLOT Full_Pathname Key_Dynamic
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global Color_outline_Udefor Color_Backgro_Defor_1 
global aveg_area_ele Time_Delay delt_time_NewMark Num_Contourlevel
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor
global Key_Ani_Ave num_Crack Width_Crack Color_Crack 
global Key_Problem Material_Para Key_Animation Num_Gauss_Points
global Key_Integral_Sol
global Num_Accuracy_Contour Key_Contour_Metd
global Output_Freq
global Color_Contourlevel Color_Mesh Key_Flipped_Gray Itera_Num Itera_HF_Time
global Na_Crack_X Na_Crack_Y num_Na_Crack Key_HF_Analysis
global frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global Num_Foc_x Num_Foc_y  Foc_x Foc_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl POS_Hl Elem_Type_Hl
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl POS_Incl Elem_Type_Incl
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Key_Time_String
global Ele_Cross_Point_RABCD
global Key_Data_Format
global Title_Font Key_Figure_Control_Widget

disp('    Generating animations of nodal stress contours......')

% Get the max and min value of displacements
disp(['    > Calculating the coordinates range of the deformed body......']) 
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
			%向量化编程
			ttt  = 1:cc_count/2;
			ttt1 = ttt*2-1;
			ttt2 = ttt*2;
			DISP(ttt,1) = ttt;
			DISP(ttt,2) = cc_DISP(ttt1);
			DISP(ttt,3) = cc_DISP(ttt2);	
		end	
		scale = Key_PLOT(3,6);

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

% Read Boundary file if necessary
if Key_PLOT(3,7) == 2 || Key_PLOT(3,7) == 3
	Boundary_X_Matrix = load([Full_Pathname,'.boux']);
	Boundary_Y_Matrix = load([Full_Pathname,'.bouy']);
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
for i=1:Num_Node
    FORCE_Matrix(i,3) = sqrt(FORCE_Matrix(i,1)^2+FORCE_Matrix(i,2)^2);
end

% Get the max and min value of stress
if Key_Animation(2)==1     % Node.
	if Key_Ani_Ave(2)==1
	    i_output=0;
		for i=1:Real_num_iteration
		    if mod(i,Output_Freq)==0
			    i_output=i_output+1;
				%读取节点应力
				if Key_Data_Format==1 
					Stress_Matrix = load([Full_Pathname,'.strn_',num2str(Itera_Num(i))]);
				elseif Key_Data_Format==2  %Binary
					c_file = fopen([Full_Pathname,'.strn_',num2str(Itera_Num(i))],'rb');
					[cc_Stress_Matrix,cc_count]   = fread(c_file,inf,'double');
					fclose(c_file);
					%转换成Matlab中的数据格式
					% for ccc_i=1:cc_count/4
						% Stress_Matrix(ccc_i,1) = ccc_i;
						% Stress_Matrix(ccc_i,2) = cc_Stress_Matrix(ccc_i*4-3);
						% Stress_Matrix(ccc_i,3) = cc_Stress_Matrix(ccc_i*4-2);
						% Stress_Matrix(ccc_i,4) = cc_Stress_Matrix(ccc_i*4-1);
						% Stress_Matrix(ccc_i,5) = cc_Stress_Matrix(ccc_i*4);
					% end
					%向量化编程
					ttt  = 1:cc_count/4;
					ttt1 = ttt*4-3;
					ttt2 = ttt*4-2;
					ttt3 = ttt*4-1;
					ttt4 = ttt*4;
					Stress_Matrix(ttt,1) = ttt;
					Stress_Matrix(ttt,2) = cc_Stress_Matrix(ttt1);
					Stress_Matrix(ttt,3) = cc_Stress_Matrix(ttt2);
					Stress_Matrix(ttt,4) = cc_Stress_Matrix(ttt3);
					Stress_Matrix(ttt,5) = cc_Stress_Matrix(ttt4);					
				end
				
				
				if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2 || Key_PLOT(3,3) ==3  % Principle stress.
				    [Stress_PS1,Stress_PS2,Node_theta] = Cal_Principal_Stress(Stress_Matrix(:,2),...
					                                                Stress_Matrix(:,3),Stress_Matrix(:,4),1);
				end
				if i_output==1
					MaxSxx =max(Stress_Matrix(:,2));MinSxx=min(Stress_Matrix(:,2));
					MaxSyy =max(Stress_Matrix(:,3));MinSyy=min(Stress_Matrix(:,3));
					MaxSxy =max(Stress_Matrix(:,4));MinSxy=min(Stress_Matrix(:,4));
					MaxSvm =max(Stress_Matrix(:,5));MinSvm=min(Stress_Matrix(:,5));
					
					if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2|| Key_PLOT(3,3) ==3   % Principle stress.
						MaxS_1 =max(Stress_PS1);MinS_1=min(Stress_PS1);
						MaxS_2 =max(Stress_PS2);MinS_2=min(Stress_PS2);				
					end
				elseif i_output > 1
					if max(Stress_Matrix(:,2)) > MaxSxx, MaxSxx=max(Stress_Matrix(:,2)); end
					if min(Stress_Matrix(:,2)) < MinSxx, MinSxx=min(Stress_Matrix(:,2)); end
					if max(Stress_Matrix(:,3)) > MaxSyy, MaxSyy=max(Stress_Matrix(:,3)); end
					if min(Stress_Matrix(:,3)) < MinSyy, MinSyy=min(Stress_Matrix(:,3)); end
					if max(Stress_Matrix(:,4)) > MaxSxy, MaxSxy=max(Stress_Matrix(:,4)); end
					if min(Stress_Matrix(:,4)) < MinSxy, MinSxy=min(Stress_Matrix(:,4)); end
					if max(Stress_Matrix(:,5)) > MaxSvm, MaxSvm=max(Stress_Matrix(:,5)); end
					if min(Stress_Matrix(:,5)) < MinSvm, MinSvm=min(Stress_Matrix(:,5)); end
					if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2  || Key_PLOT(3,3) ==3 % Principle stress.
						if max(Stress_PS1) > MaxS_1, MaxS_1=max(Stress_PS1); end
						if min(Stress_PS1) < MinS_1, MinS_1=min(Stress_PS1); end
						if max(Stress_PS2) > MaxS_2, MaxS_2=max(Stress_PS2); end
						if min(Stress_PS2) < MinS_2, MinS_2=min(Stress_PS2); end					
					end
				end
				clear Stress_Matrix
			end
		end
	end
elseif Key_Animation(2)==2 % Gauss points.
    if Key_Ani_Ave(2)==1
	    i_output=0;
	    for i=1:Real_num_iteration
		    i_output = i_output+1;
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
					%向量化编程
					ttt  = 1:cc_count/2;
					ttt1 = ttt*2-1;
					ttt2 = ttt*2;
					DISP(ttt,1) = ttt;
					DISP(ttt,2) = cc_DISP(ttt1);
					DISP(ttt,3) = cc_DISP(ttt2);								
				end	
				
				% Read gauss stress.
				if Key_Data_Format==1 
			  	    tem_Gauss_Stress = load([Full_Pathname,'.strg_',num2str(Itera_Num(i))]);
				elseif Key_Data_Format==2  %Binary
					c_file = fopen([Full_Pathname,'.strg_',num2str(Itera_Num(i))],'rb');
					[cc_Gauss_Stress,cc_count]   = fread(c_file,inf,'double');
					fclose(c_file);
					%转换成Matlab中的数据格式
					% for ccc_i=1:cc_count/4
						% tem_Gauss_Stress(ccc_i,1) = ccc_i;
						% tem_Gauss_Stress(ccc_i,2) = cc_Gauss_Stress(ccc_i*4-3);
						% tem_Gauss_Stress(ccc_i,3) = cc_Gauss_Stress(ccc_i*4-2);
						% tem_Gauss_Stress(ccc_i,4) = cc_Gauss_Stress(ccc_i*4-1);
						% tem_Gauss_Stress(ccc_i,5) = cc_Gauss_Stress(ccc_i*4);
					% end
					%向量化编程
					ttt  = 1:cc_count/4;
					ttt1 = ttt*4-3;
					ttt2 = ttt*4-2;
					ttt3 = ttt*4-1;
					ttt4 = ttt*4;
					tem_Gauss_Stress(ttt,1) = ttt;
					tem_Gauss_Stress(ttt,2) = cc_Gauss_Stress(ttt1);
					tem_Gauss_Stress(ttt,3) = cc_Gauss_Stress(ttt2);
					tem_Gauss_Stress(ttt,4) = cc_Gauss_Stress(ttt3);
					tem_Gauss_Stress(ttt,5) = cc_Gauss_Stress(ttt4);		
				end
				
				
				
				% Gauss点数目
				num_Total_Gauss = size(tem_Gauss_Stress,1);
				scale = Key_PLOT(2,6);
				% Read coordinates and all other information of cracks if cracks exist.
				% 111,num_Crack(i) 
				if num_Crack(i) ~= 0
				    disp('    > Reading coordinates and other information of cracks....') 
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
					
			    	% tem_Ele_Cross_Point_RABCD = load([Full_Pathname,'.crab_',num2str(Itera_Num(i))]); %2017-05-04
				    % tem_num_ele = size(tem_Ele_Cross_Point_RABCD,1)/2;
				    % Ele_Cross_Point_RABCD(1:tem_num_ele,1:10,1) = tem_Ele_Cross_Point_RABCD(1:2:tem_num_ele*2,1:10);
				    % Ele_Cross_Point_RABCD(1:tem_num_ele,1:10,2) = tem_Ele_Cross_Point_RABCD(2:2:tem_num_ele*2,1:10);
					% Node_Cross_elem = load([Full_Pathname,'.ncel_',num2str(Itera_Num(i))]);%Cross增强节点对应的Cross单元号,added on 2017-05-04
					Node_Cross_elem =[];
					
	                % Node_Jun_elem = load([Full_Pathname,'.njel_',num2str(Itera_Num(i))]);%Junction增强节点对应的Junction单元号,added on 2016-07-11
					% Node_Jun_Hole = load([Full_Pathname,'.njhl_',num2str(Itera_Num(i))]);%Junction增强节点(Crack and hole)对应的Hole号
					if Key_Data_Format==1 
						Node_Jun_elem = load([Full_Pathname,'.njel_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元号,added on 2016-07-11
					elseif Key_Data_Format==2  %Binary
						c_file = fopen([Full_Pathname,'.njel_',num2str(Itera_Num(i))],'rb');
						[cc_Node_Jun_elem,cc_count]   = fread(c_file,inf,'int');
						fclose(c_file);
						%转换成Matlab中的数据格式
						Node_Jun_elem = (reshape(cc_Node_Jun_elem,num_Crack(Itera_Num(i)),Num_Node))';
					end	
					if Key_Data_Format==1 
						Node_Jun_Hole = load([Full_Pathname,'.njhl_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元号,added on 2016-07-11
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

					if Key_Data_Format==1 
						Coors_Tip           = load([Full_Pathname,'.celt_',num2str(Itera_Num(i))]);
					elseif Key_Data_Format==2  %Binary
						c_file = fopen([Full_Pathname,'.celt_',num2str(Itera_Num(i))],'rb');
						[cc_Coors_Tip,cc_count]   = fread(c_file,inf,'double');
						fclose(c_file);
						%转换成Matlab中的数据格式
						Coors_Tip = (reshape(cc_Coors_Tip,2,num_ele))';
					end							
				else
					Crack_X = [];   Crack_Y = [];
					Enriched_Node_Type = [];
					Elem_Type = [];
					POS = []; Coors_Element_Crack= [];Coors_Vertex= [];Node_Jun_elem=[];Node_Cross_elem=[];
					Coors_Junction= []; Coors_Tip= []; Elem_Type= [];Node_Jun_Hole=[];
				end
				%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				% Read coordinates and other info of holes if holes exist.
				%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				if num_Hole ~= 0
					disp('    > Reading coordinates of hole....') 
					file_Hole = fopen([Full_Pathname,'.hlcr']);
					for iHole=1:num_Hole
						Hole_Coor(iHole,1:3)= str2num(fgetl(file_Hole));
					end
					fclose(file_Hole);
					disp('    > Reading ennh file of hole....') 
					Enriched_Node_Type_Hl = load([Full_Pathname,'.ennh_',num2str(Itera_Num(i))]);
					disp('    > Reading posh file of hole....') 
					POS_Hl = load([Full_Pathname,'.posh_',num2str(Itera_Num(i))]);
					
					disp('    > Reading elth file of hole....');
					Elem_Type_Hl = load([Full_Pathname,'.elth_',num2str(Itera_Num(i))]);
				end	
				% Get the stress of all Gauss points.
				Stress_xx(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,2);
				Stress_yy(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,3);
				Stress_xy(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,4);
				Stress_vm(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,5);
						
			    % Calculate principal stress if Key_PLOT(3,3) ==1 or 2 or 3.
				if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2|| Key_PLOT(3,3) ==3 
					[Stress_PS1,Stress_PS2,Gauss_theta] = Cal_Principal_Stress(Stress_xx,Stress_yy,Stress_xy,1);
				end
			
				clear tem_Gauss_Stress
				if i_output==1
					MaxSxx =max(Stress_xx);MinSxx=min(Stress_xx);
					MaxSyy =max(Stress_yy);MinSyy=min(Stress_yy);
					MaxSxy =max(Stress_xy);MinSxy=min(Stress_xy);
					MaxSvm =max(Stress_vm);MinSvm=min(Stress_vm);
					if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2 || Key_PLOT(3,3) ==3  % Principle stress.
						MaxS_1 =max(Stress_PS1);MinS_1=min(Stress_PS1);
						MaxS_2 =max(Stress_PS2);MinS_2=min(Stress_PS2);				
					end
				elseif i_output > 1
					if max(Stress_xx) > MaxSxx, MaxSxx=max(Stress_xx); end
					if min(Stress_xx) < MinSxx, MinSxx=min(Stress_xx); end
					if max(Stress_yy) > MaxSyy, MaxSyy=max(Stress_yy); end
					if min(Stress_yy) < MinSyy, MinSyy=min(Stress_yy); end
					if max(Stress_xy) > MaxSxy, MaxSxy=max(Stress_xy); end
					if min(Stress_xy) < MinSxy, MinSxy=min(Stress_xy); end
					if max(Stress_vm) > MaxSvm, MaxSvm=max(Stress_vm); end
					if min(Stress_vm) < MinSvm, MinSvm=min(Stress_vm); end
					if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2|| Key_PLOT(3,3) ==3   % Principle stress.
						if max(Stress_PS1) > MaxS_1, MaxS_1=max(Stress_PS1); end
						if min(Stress_PS1) < MinS_1, MinS_1=min(Stress_PS1); end
						if max(Stress_PS2) > MaxS_2, MaxS_2=max(Stress_PS2); end
						if min(Stress_PS2) < MinS_2, MinS_2=min(Stress_PS2); end					
					end
				end
				clear Stress_xx;clear Stress_yy;clear Stress_xy;clear Stress_vm;
            end
		end
	end
end

% Loop through each frame to get the max pressure and quantity of all frame for hydraulic fracturing simulation, i.e.	
% Max_Pressure_for_aimation and Max_Quantity_for_aimation.
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

% Loop through each frame to plot and store
i_output =0;
for i=1:Real_num_iteration
    if mod(i,Output_Freq)==0
	    i_output = i_output +1;
		disp(['    > Plotting and saving stress contours of frame ',num2str(i),'......']) 
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
			%向量化编程
			ttt  = 1:cc_count/2;
			ttt1 = ttt*2-1;
			ttt2 = ttt*2;
			DISP(ttt,1) = ttt;
			DISP(ttt,2) = cc_DISP(ttt1);
			DISP(ttt,3) = cc_DISP(ttt2);						
		end		
		
		% Read kiel file.
		num_Killed_Ele=0;
		if exist([Full_Pathname,'.kiel_',num2str(Itera_Num(i))], 'file') ==2  
			Killed_Elements = load([Full_Pathname,'.kiel_',num2str(Itera_Num(i))]);
			num_Killed_Ele  = size(Killed_Elements,1);  %被杀死单元的数目
		end
		
		%读取节点应力
		if Key_Data_Format==1 
			Stress_Matrix = load([Full_Pathname,'.strn_',num2str(Itera_Num(i))]);
		elseif Key_Data_Format==2  %Binary
			c_file = fopen([Full_Pathname,'.strn_',num2str(Itera_Num(i))],'rb');
			[cc_Stress_Matrix,cc_count]   = fread(c_file,inf,'double');
			fclose(c_file);
			%转换成Matlab中的数据格式
			% for ccc_i=1:cc_count/4
				% Stress_Matrix(ccc_i,1) = ccc_i;
				% Stress_Matrix(ccc_i,2) = cc_Stress_Matrix(ccc_i*4-3);
				% Stress_Matrix(ccc_i,3) = cc_Stress_Matrix(ccc_i*4-2);
				% Stress_Matrix(ccc_i,4) = cc_Stress_Matrix(ccc_i*4-1);
				% Stress_Matrix(ccc_i,5) = cc_Stress_Matrix(ccc_i*4);
			% end
			%向量化编程
			ttt  = 1:cc_count/4;
			ttt1 = ttt*4-3;
			ttt2 = ttt*4-2;
			ttt3 = ttt*4-1;
			ttt4 = ttt*4;
			Stress_Matrix(ttt,1) = ttt;
			Stress_Matrix(ttt,2) = cc_Stress_Matrix(ttt1);
			Stress_Matrix(ttt,3) = cc_Stress_Matrix(ttt2);
			Stress_Matrix(ttt,4) = cc_Stress_Matrix(ttt3);
			Stress_Matrix(ttt,5) = cc_Stress_Matrix(ttt4);			
		end		
		
		%如果是绘制Gauss点应力
		if Key_Animation(2)==2
			% 读取Gauss点位移.
			if Key_Data_Format==1 
			    tem_Gauss_dis = load([Full_Pathname,'.disg_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.disg_',num2str(Itera_Num(i))],'rb');
				[cc_tem_Gauss_dis,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				% for ccc_i=1:cc_count/2
					% tem_Gauss_dis(ccc_i,1) = ccc_i;
					% tem_Gauss_dis(ccc_i,2) = cc_tem_Gauss_dis(ccc_i*2-1);
					% tem_Gauss_dis(ccc_i,3) = cc_tem_Gauss_dis(ccc_i*2);
				% end
				%向量化编程
				ttt  = 1:cc_count/2;
				ttt1 = ttt*2-1;
				ttt2 = ttt*2;
				tem_Gauss_dis(ttt,1) = ttt;
				tem_Gauss_dis(ttt,2) = cc_tem_Gauss_dis(ttt1);
				tem_Gauss_dis(ttt,3) = cc_tem_Gauss_dis(ttt2);					
			end
			% Read gauss stress.
			if Key_Data_Format==1 
				tem_Gauss_Stress = load([Full_Pathname,'.strg_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.strg_',num2str(Itera_Num(i))],'rb');
				[cc_Gauss_Stress,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				% for ccc_i=1:cc_count/4
					% tem_Gauss_Stress(ccc_i,1) = ccc_i;
					% tem_Gauss_Stress(ccc_i,2) = cc_Gauss_Stress(ccc_i*4-3);
					% tem_Gauss_Stress(ccc_i,3) = cc_Gauss_Stress(ccc_i*4-2);
					% tem_Gauss_Stress(ccc_i,4) = cc_Gauss_Stress(ccc_i*4-1);
					% tem_Gauss_Stress(ccc_i,5) = cc_Gauss_Stress(ccc_i*4);
				% end
				%向量化编程
				ttt  = 1:cc_count/4;
				ttt1 = ttt*4-3;
				ttt2 = ttt*4-2;
				ttt3 = ttt*4-1;
				ttt4 = ttt*4;
				tem_Gauss_Stress(ttt,1) = ttt;
				tem_Gauss_Stress(ttt,2) = cc_Gauss_Stress(ttt1);
				tem_Gauss_Stress(ttt,3) = cc_Gauss_Stress(ttt2);
				tem_Gauss_Stress(ttt,4) = cc_Gauss_Stress(ttt3);
				tem_Gauss_Stress(ttt,5) = cc_Gauss_Stress(ttt4);					
			end			
			
		    % Gauss点数目
		    num_Total_Gauss = size(tem_Gauss_Stress,1);
			% 读取Gauss点的坐标
			if Key_Data_Format==1 
                tem_Gauss_cor = load([Full_Pathname,'.gcor_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.gcor_',num2str(Itera_Num(i))],'rb');
				[cc_tem_Gauss_cor,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
				%转换成Matlab中的数据格式
				% for ccc_i=1:cc_count/2
					% tem_Gauss_cor(ccc_i,1) = ccc_i;
					% tem_Gauss_cor(ccc_i,2) = cc_tem_Gauss_cor(ccc_i*2-1);
					% tem_Gauss_cor(ccc_i,3) = cc_tem_Gauss_cor(ccc_i*2);
				% end
				%向量化编程
				ttt  = 1:cc_count/2;
				ttt1 = ttt*2-1;
				ttt2 = ttt*2;
				tem_Gauss_cor(ttt,1) = ttt;
				tem_Gauss_cor(ttt,2) = cc_tem_Gauss_cor(ttt1);
				tem_Gauss_cor(ttt,3) = cc_tem_Gauss_cor(ttt2);						
			end
		end
	
		scale = Key_PLOT(2,6);

		% Read coordinates and all other information of cracks if cracks exist.
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
			
			% tem_Ele_Cross_Point_RABCD = load([Full_Pathname,'.crab_',num2str(Itera_Num(i))]); %2017-05-04
			% tem_num_ele = size(tem_Ele_Cross_Point_RABCD,1)/2;
			% Ele_Cross_Point_RABCD(1:tem_num_ele,1:10,1) = tem_Ele_Cross_Point_RABCD(1:2:tem_num_ele*2,1:10);
			% Ele_Cross_Point_RABCD(1:tem_num_ele,1:10,2) = tem_Ele_Cross_Point_RABCD(2:2:tem_num_ele*2,1:10);			
			% Node_Cross_elem     = load([Full_Pathname,'.ncel_',num2str(Itera_Num(i))]);
			Node_Cross_elem = [];
			
			% Node_Jun_elem       = load([Full_Pathname,'.njel_',num2str(Itera_Num(i))]);%Junction增强节点对应的Junction单元号,added on 2016-07-11
			% Node_Jun_Hole       = load([Full_Pathname,'.njhl_',num2str(Itera_Num(i))]);
			if Key_Data_Format==1 
				Node_Jun_elem = load([Full_Pathname,'.njel_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元号,added on 2016-07-11
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.njel_',num2str(Itera_Num(i))],'rb');
				[cc_Node_Jun_elem,cc_count]   = fread(c_file,inf,'int');
				fclose(c_file);
				%转换成Matlab中的数据格式
				Node_Jun_elem = (reshape(cc_Node_Jun_elem,num_Crack(Itera_Num(i)),Num_Node))';
			end	
			if Key_Data_Format==1 
				Node_Jun_Hole = load([Full_Pathname,'.njhl_',num2str(Itera_Num(i))]); %Junction增强节点对应的Junction单元号,added on 2016-07-11
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
				tem_Coors_Junction  = load([Full_Pathname,'.celj_',num2str(Itera_Num(i))]);
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
		else
			Crack_X = [];   Crack_Y = [];
			Enriched_Node_Type = [];
			Crack_Tip_Type     = [];
			Elem_Type = [];
			POS = []; Coors_Element_Crack= [];Coors_Vertex= [];Node_Jun_elem=[];Node_Cross_elem=[];
			Coors_Junction= []; Coors_Tip= []; Elem_Type= []; Node_Jun_Hole =[];
		end
		
		% Read coordinates of arc (Arc_Crack_Coor(,,1:11)) crack if arc crack exist, 2017-07-26.
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
		else
			Yes_Arc_Crack  = 0;
			Arc_Crack_Coor(1:1000,1:1000,1:11) = 0.0;	
		end
		
		% Read coordinates and other info of holes if holes exist.
		if num_Hole ~= 0
			disp('    > Reading coordinates of hole....') 
			file_Hole = fopen([Full_Pathname,'.hlcr']);
			for iHole=1:num_Hole
				Hole_Coor(iHole,1:3)= str2num(fgetl(file_Hole));
			end
			fclose(file_Hole);
			disp('    > Reading ennh file of hole....') 
			Enriched_Node_Type_Hl = load([Full_Pathname,'.ennh_',num2str(Itera_Num(i))]);
			disp('    > Reading posh file of hole....') 
			POS_Hl = load([Full_Pathname,'.posh_',num2str(Itera_Num(i))]);
					
			disp('    > Reading elth file of hole....');
			Elem_Type_Hl = load([Full_Pathname,'.elth_',num2str(Itera_Num(i))]);
		end	
		% 读取圆形夹杂的坐标.
		if num_Circ_Inclusion ~= 0
			disp('    > Reading coordinates of circle inclusions....') 
			file_Circ_Inclusion = fopen([Full_Pathname,'.jzcr']);
			for iiii=1:num_Circ_Inclusion
				Circ_Inclu_Coor(iiii,1:3)= str2num(fgetl(file_Circ_Inclusion));
			end
			fclose(file_Circ_Inclusion);
			disp('    > Reading ennh file of circle inclusions....') 
			Enriched_Node_Type_Incl = load([Full_Pathname,'.ennj_',num2str(Itera_Num(i))]);
			disp('    > Reading posh file of circle inclusions....') 
			POS_Incl = load([Full_Pathname,'.posj_',num2str(Itera_Num(i))]);
			disp('    > Reading elth file of circle inclusions....');
			Elem_Type_Incl = load([Full_Pathname,'.eltj_',num2str(Itera_Num(i))]);
		end

		% 读取多边形夹杂的坐标,2016-10-04.
		if num_Poly_Inclusion ~= 0
			disp('    > Reading coordinates of polygon inclusions....') 
			file_Poly_Incl_Coor_x = fopen([Full_Pathname,'.jzpx']);
			file_Poly_Incl_Coor_y = fopen([Full_Pathname,'.jzpy']);
			Poly_Incl_Coor_x = cell(num_Poly_Inclusion);
			Poly_Incl_Coor_y = cell(num_Poly_Inclusion);
			for iiii=1:num_Poly_Inclusion
				Poly_Incl_Coor_x{iiii} = str2num(fgetl(file_Poly_Incl_Coor_x));
				Poly_Incl_Coor_y{iiii} = str2num(fgetl(file_Poly_Incl_Coor_y));
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
		% Prepare for Gauss stress plot.
		if Key_Animation(2)==2
			Stress_xx(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,2);
			Stress_yy(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,3);
			Stress_xy(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,4);
			Stress_vm(1:num_Total_Gauss) = tem_Gauss_Stress(1:num_Total_Gauss,5);
				
			dis_x(1:num_Total_Gauss) = tem_Gauss_dis(1:num_Total_Gauss,2);
			dis_y(1:num_Total_Gauss) = tem_Gauss_dis(1:num_Total_Gauss,3);		
			
			New_Gauss_Coor(1:num_Total_Gauss,1) = tem_Gauss_cor(1:num_Total_Gauss,2) + scale*dis_x(1:num_Total_Gauss)';
			New_Gauss_Coor(1:num_Total_Gauss,2) = tem_Gauss_cor(1:num_Total_Gauss,3) + scale*dis_y(1:num_Total_Gauss)';	
			
			clear tem_Gauss_dis
			clear tem_Gauss_Stress
			clear tem_Gauss_cor
			
			% Calculate principal stress if Key_PLOT(3,3) ==1 or 2 or 3.
			if Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2|| Key_PLOT(3,3) ==3 
				[Stress_PS1,Stress_PS2,Gauss_theta] = Cal_Principal_Stress(Stress_xx,Stress_yy,Stress_xy,1);
			end
		end
		
		% Get the new coordinates of all nodes
		New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
		New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);
		
		% Calculate shaped cracks if necessary.
		if Key_PLOT(3,5) == 2 & num_Crack(i)~=0
			disp(['      > Calculating shaped crack points......'])
			[Shaped_Crack_Points] = Cal_Shaped_Cracks(Crack_X,Crack_Y,i,Itera_Num(i),num_Crack,Crack_Tip_Type,POS,...
										   Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
										   Coors_Vertex,Coors_Junction,Coors_Tip,DISP,scale);
		end
		
		% Get resample coors.
		delta = sqrt(aveg_area_ele)/Num_Accuracy_Contour;
		if Key_PLOT(3,8)==0
			gx    = Min_X_Coor:delta:Max_X_Coor; 
			gy    = Min_Y_Coor:delta:Max_Y_Coor;
		elseif Key_PLOT(3,8)==1
			gx    = Last_Min_X:delta:Last_Max_X; 
			gy    = Last_Min_Y:delta:Last_Max_Y; 
		end
		
		% Plot the stress contours.
		Start_Stress_Type =1;      % Plot Sxx, Sxy, Syy and Mises.
		End_Stress_Type   =4;    
		if Key_PLOT(3,2) ==1       % Only plot Mises stress.
			Start_Stress_Type =4;
			End_Stress_Type   =4;
		end
		if Key_PLOT(3,2) ==2       % Only plot stress-x.
			Start_Stress_Type =1;
			End_Stress_Type   =1;
		end
		if Key_PLOT(3,2) ==3       % Only plot stress-y.
			Start_Stress_Type =2;
			End_Stress_Type   =2;
		end
		if Key_PLOT(3,2) ==4       % Only plot stress-xy.
			Start_Stress_Type =3;
			End_Stress_Type   =3;
		end
		if (Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2) && Key_PLOT(3,2)==1       % Plot principle stress and misses.
			Start_Stress_Type =4;
			End_Stress_Type   =6;
		elseif (Key_PLOT(3,3) ==1 || Key_PLOT(3,3) ==2) && Key_PLOT(3,2)==0   % Plot principle stress.
			Start_Stress_Type =5;
			End_Stress_Type   =6;
		end
		if Key_PLOT(3,3) ==3    % Only plot max principle stress.
			Start_Stress_Type =5;
			End_Stress_Type   =5;
		end
		for j=Start_Stress_Type:End_Stress_Type
			% Resample nodal stress.
			if Key_Animation(2)==1
				if Key_PLOT(3,8)==0
					% Resample with a grid, THIS IS VERY SLOW if Node_Coor is large
					if j==1
						disp('      Resample Sxx....')
						 [Sxx,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_Matrix(:,2),gx,gy);					
					elseif j==2
						disp('      Resample Syy....')				
						[Syy,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_Matrix(:,3),gx,gy);						
					elseif j==3	
						disp('      Resample Sxy....')				
						 [Sxy,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_Matrix(:,4),gx,gy);					
					elseif j==4
						disp('      Resample Svm....')				
						[Svm,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_Matrix(:,5),gx,gy);	
					elseif j==5					 
						disp('      Resample major principal stress....') 	
						[Stress_1,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_PS1(:),gx,gy);		
					elseif j==6					 
						disp('      Resample minor principal stress....') 	
						[Stress_2,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Stress_PS2(:),gx,gy);						
					end
				elseif Key_PLOT(3,8)==1
					% Resample with a grid, THIS IS VERY SLOW if Node_Coor is large
					if j==1
						disp('      Resample Sxx....')
						[Sxx,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_Matrix(:,2),gx,gy);						 
					elseif j==2
						disp('      Resample Syy....')				
						[Syy,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_Matrix(:,3),gx,gy);					 
					elseif j==3	
						disp('      Resample Sxy....')				
						[Sxy,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_Matrix(:,4),gx,gy);					 
					elseif j==4	
						disp('      Resample Svm....')				
						[Svm,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_Matrix(:,5),gx,gy);
					elseif j==5					 
						disp('      Resample major principal stress....') 	
						[Stress_1,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_PS1(:),gx,gy);					
					elseif j==6					 
						disp('      Resample minor principal stress....') 	
						[Stress_2,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Stress_PS2(:),gx,gy);						
					end
				end
			% Resample Gauss stress.
			elseif Key_Animation(2)==2
				if Key_PLOT(3,8)==0
					% Resample with a grid, THIS IS VERY SLOW if Gauss_Coor is large
					if j==1
						disp('      Resample Sxx....') 
						[Sxx,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_xx(:),gx,gy);					 
					elseif j==2					 
						disp('      Resample Syy....') 
						[Syy,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_yy(:),gx,gy);					 
					elseif j==3					 
						disp('      Resample Sxy....') 
						[Sxy,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_xy(:),gx,gy);					 
					elseif j==4					 
						disp('      Resample Svm....') 
						[Svm,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_vm(:),gx,gy);
					elseif j==5					 
						disp('      Resample major principal stress....') 	
						[Stress_1,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_PS1(:),gx,gy);		
					elseif j==6					 
						disp('      Resample minor principal stress....') 	
						[Stress_2,X,Y] = Tools_gridfit(Gauss_Coor(:,1),Gauss_Coor(:,2),Stress_PS2(:),gx,gy);								
					end
				elseif Key_PLOT(3,8)==1
					% Resample with a grid, THIS IS VERY SLOW if Gauss_Coor is large
					if j==1
						disp('      Resample Sxx....') 			
						[Sxx,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_xx(:),gx,gy);					 
					elseif j==2
						disp('      Resample Syy....') 			
						[Syy,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_yy(:),gx,gy);					 
					elseif j==3
						disp('      Resample Sxy....') 			
						[Sxy,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_xy(:),gx,gy);					 
					elseif j==4
						disp('      Resample Svm....')			 
						[Svm,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_vm(:),gx,gy);	
					elseif j==5					 
						disp('      Resample major principal stress....') 	
						[Stress_1,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_PS1(:),gx,gy);					
					elseif j==6					 
						disp('      Resample minor principal stress....') 	
						[Stress_2,X,Y] = Tools_gridfit(New_Gauss_Coor(:,1),New_Gauss_Coor(:,2),Stress_PS2(:),gx,gy);								
					end
				end
			end
			% New figure.
			Tools_New_Figure
			
			hold on;
			if j == 1
				disp('      Contouring Sxx....') 
				contourf(X,Y,Sxx,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_x_x','FontName',Title_Font,'FontSize',Size_Font);
				clear Sxx
			elseif j == 2
				disp('      Contouring Syy....') 
				contourf(X,Y,Syy,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_y_y','FontName',Title_Font,'FontSize',Size_Font);
				clear Syy
			elseif j == 3
				disp('      Contouring Sxy....') 
				contourf(X,Y,Sxy,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_x_y','FontName',Title_Font,'FontSize',Size_Font);
				clear Sxy
			elseif j == 4
				disp('      Contouring Svm....') 
				contourf(X,Y,Svm,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_v_m','FontName',Title_Font,'FontSize',Size_Font);
				clear Svm
			elseif j == 5
				disp('      Contouring major principal stress....') 
				contourf(X,Y,Stress_1,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_1','FontName',Title_Font,'FontSize',Size_Font);
				clear Stress_1
			elseif j == 6
				disp('      Contouring minor principal stress....') 
				contourf(X,Y,Stress_2,Num_Contourlevel,'LineStyle','none')
				title('Stress plot, \sigma_2','FontName',Title_Font,'FontSize',Size_Font);
				clear Stress_2
			end
			% Set colormap.
			if Key_Flipped_Gray==0
				colormap(Color_Contourlevel)
			elseif Key_Flipped_Gray==1
				colormap(flipud(gray))
			end
			% Control of the range of the colour bar
			if Key_Ani_Ave(2)==1
				if j==1
					caxis([MinSxx, MaxSxx]); 
				elseif j==2
					caxis([MinSxy, MaxSxy]); 
				elseif j==3
					caxis([MinSyy, MaxSyy]); 
				elseif j==4
					caxis([MinSvm, MaxSvm]); 
				elseif j==5
					caxis([MinS_1, MaxS_1]); 
				elseif j==6
					caxis([MinS_2, MaxS_2]); 
				end
			end
			% 如果各用各的云图绘制范围,且是Gauss点云图，则需要修正,因为有可能网格化之后(绘图方法2)的数据范围发生了变化
			% if Key_Ani_Ave(2)==0 && Key_Animation(2)==2
				% if j==1
					% caxis([min(Stress_xx), max(Stress_xx)]); 
				% elseif j==2
					% caxis([min(Stress_yy), max(Stress_yy)]); 
				% elseif j==3
					% caxis([min(Stress_xy), max(Stress_xy)]); 
				% elseif j==4
					% caxis([min(Stress_vm), max(Stress_vm)]); 
				% elseif j==5
					% caxis([min(Stress_PS1), max(Stress_PS1)]); 
				% elseif j==6
					% caxis([min(Stress_PS2), max(Stress_PS2)]); 
				% end
			% end
			% Plot the mesh
			if Key_PLOT(3,9) ==1
				if Key_PLOT(3,8)==0
					for iElem =1:Num_Elem
						NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
							  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)]; % Nodes for current element
						xi = Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
						yi = Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
						plot(xi,yi,'LineWidth',1,'Color',Color_Mesh)
					end
				elseif Key_PLOT(3,8)==1
					for iElem =1:Num_Elem
						NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
							  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)];     % Nodes for current element
						xi = New_Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
						yi = New_Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
						plot(xi,yi,'LineWidth',1,'Color',Color_Mesh)
					end
				end
			end
			%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			% Plot boundary conditions by darwing triangles
			%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if Key_PLOT(3,7) == 2 || Key_PLOT(3,7) == 3
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
			if Key_PLOT(3,7) == 1 || Key_PLOT(3,7) == 3
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


			%<<<<<<<<<<<<<<<<<<<
			% 绘制天然裂缝.
			%<<<<<<<<<<<<<<<<<<<
			if Key_PLOT(3,12) == 1
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
							
							% plot(last_x,last_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   
							plot(last_x,last_y,'w','LineWidth',1,'Color','white') 
							% plot(x,y,'--','Color',Color_Crack,'LineWidth',Width_Crack)  
						end
					end	
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
							Line_Edge(1,1) = Poly_Incl_Coor_x{iii}(iEdge); %边线的起点
							Line_Edge(1,2) = Poly_Incl_Coor_y{iii}(iEdge);
							Line_Edge(2,1) = Poly_Incl_Coor_x{iii}(1);     %边线的终点
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
						%计算边线起点的位移
						cc_k =cc_k+1;
						cc_x(cc_k) = a_x;
						cc_y(cc_k) = a_y;
						[Kesi,Yita] = Cal_KesiYita_by_Coors(cc_x(cc_k),cc_y(cc_k));
						[c_Elem_Num] = Cal_Ele_Num_by_Coors(cc_x(cc_k),cc_y(cc_k));
						[cc_dis_x(cc_k),cc_dis_y(cc_k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																		 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																		  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
						%计算等分点的位移
						for cc_j =  1:Num_Diversion-1
							cc_k=cc_k+1;
							cc_x(cc_k) = (cc_j*b_x+(Num_Diversion-cc_j)*a_x)/Num_Diversion;
							cc_y(cc_k) = (cc_j*b_y+(Num_Diversion-cc_j)*a_y)/Num_Diversion;
							[Kesi,Yita] = Cal_KesiYita_by_Coors(cc_x(cc_k),cc_y(cc_k));
							[c_Elem_Num] = Cal_Ele_Num_by_Coors(cc_x(cc_k),cc_y(cc_k));
							[cc_dis_x(cc_k),cc_dis_y(cc_k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,Itera_Num(i),DISP,Kesi,Yita...
																		 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																		  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
						end												      
					end
					x_new = cc_x + cc_dis_x*scale;
					y_new = cc_y + cc_dis_y*scale;
					................
					% option1:填充
					%................
					% patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
					................
					% option2:绘制边线
					%................
					x_new(cc_k+1) = x_new(1);
					y_new(cc_k+1) = y_new(1);
					plot(x_new,y_new,'-','color','black')
				end	
			end	
			% The range of the plot.
			delta = sqrt(aveg_area_ele);
			if Key_PLOT(3,8)==0
				min_x = Min_X_Coor-delta; max_x = Max_X_Coor+delta; 
				min_y = Min_Y_Coor-delta; max_y = Max_Y_Coor+delta;
			elseif Key_PLOT(3,8)==1
				min_x = Last_Min_X-delta; max_x = Last_Max_X+delta;
				min_y = Last_Min_Y-delta; max_y = Last_Max_Y+delta;
			end
			axis([min_x max_x min_y max_y]);	
			
			% Fill the outside of the domain by white color.
			if Key_PLOT(3,8)==0
				Tools_Fillout(Node_Coor(Outline(:,1),1),Node_Coor(Outline(:,1),2),[min_x max_x min_y max_y],'w'); 
			elseif Key_PLOT(3,8)==1
				Tools_Fillout(New_Node_Coor(Outline(:,1),1),New_Node_Coor(Outline(:,1),2),[min_x max_x min_y max_y],'w'); 
			end
			
			% Fill the inside of the domain by white color, holes, for example.
			if isempty(Inline)==0
				if Key_PLOT(3,8)==0
					fill(Node_Coor(Inline(:,1),1),Node_Coor(Inline(:,1),2),'w');
				elseif Key_PLOT(3,8)==1 
					fill(New_Node_Coor(Inline(:,1),1),New_Node_Coor(Inline(:,1),2),'w')	
				end
			end
			
			% Plot the outline of the mesh.
			if Key_PLOT(3,8)==0
				line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
					 [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
			elseif Key_PLOT(3,8)==1
				line([New_Node_Coor(Outline(:,1),1) New_Node_Coor(Outline(:,2),1)]', ...
					 [New_Node_Coor(Outline(:,1),2) New_Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
			end
			
			axis off; 
			axis equal;

			%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
			colorbar('FontName',Title_Font,'FontSize',Size_Font);
			
			% Plot text on the left or the bottom of the figure
			plot_string_MatMES = ['PhiPsi  ',Version];
			plot_string_Frame  = ['Frame ',num2str(i_output),' / ',num2str(Real_num_iteration)];
			% plot_string_Time   = ['Time: ',num2str(i*delt_time_NewMark*1000),' ms'];
			% Plot text on the left or the bottom of the figure.
			plot_string_MatMES = ['PhiPsi ',Version];
			plot_string_Frame  = ['Frame ',num2str(i_output),' / ',num2str(Real_num_iteration)];
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
			plot_string_Scale  = ['Scale factor: ',num2str(scale)];
		
			% if Key_HF_Analysis==1
				% plot_string_Time_HF   = ['Time: ',num2str(Itera_HF_Time(i)/60.0,'%0.4f'),' mins'];
			% end
			plot_string_Scale  = ['Scale factor: ',num2str(scale)];
			if Key_PLOT(3,8)==0
				range_W = abs(Max_X_Coor-Min_X_Coor);
				range_H = abs(Max_Y_Coor-Min_Y_Coor);
				if range_H >= 0.75*range_W         % Left
					loc_x = -range_H/2+ Min_X_Coor;
					loc_y =  Max_Y_Coor-range_H*0.05;
					text(loc_x,loc_y,plot_string_MatMES,'color','black');
					loc_y =  Max_Y_Coor-range_H*0.15;
					text(loc_x,loc_y,plot_string_Scale,'color','black');
					loc_y =  Max_Y_Coor-range_H*0.25;
					text(loc_x,loc_y,plot_string_Frame,'color','black');
					if Key_Dynamic==1
						loc_y =  Max_Y_Coor-range_H*0.35;
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
				else                               % Bottom
					loc_y =  Min_Y_Coor - range_H*0.2;
					loc_x =  Min_X_Coor;
					text(loc_x,loc_y,plot_string_MatMES,'color','black');
					loc_x =  Min_X_Coor + range_W*0.25;
					text(loc_x,loc_y,plot_string_Scale,'color','black');
					loc_x =  Min_X_Coor + range_W*0.45;
					text(loc_x,loc_y,plot_string_Frame,'color','black');
					if Key_Dynamic==1
						loc_x =  Min_X_Coor + range_W*0.60;
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
			elseif Key_PLOT(3,8)==1
				range_W = abs(Last_Max_X-Last_Min_X);
				range_H = abs(Last_Max_Y-Last_Min_Y);
				if range_H >= 0.75*range_W             % Left
					loc_x = -range_H/2+ Last_Min_X;
					loc_y =  Last_Max_Y-range_H*0.05;
					text(loc_x,loc_y,plot_string_MatMES,'color','black');
					loc_y =  Last_Max_Y-range_H*0.15;
					text(loc_x,loc_y,plot_string_Scale,'color','black');
					loc_y =  Last_Max_Y-range_H*0.25;
					text(loc_x,loc_y,plot_string_Frame,'color','black');
					if Key_Dynamic==1
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
				else                                   % Bottom
					loc_y =  Last_Min_Y - range_H*0.05; 
					loc_x =  Last_Min_X;
					text(loc_x,loc_y,plot_string_MatMES,'color','black');
					loc_x =  Last_Min_X + range_W*0.25;
					text(loc_x,loc_y,plot_string_Scale,'color','black');
					loc_x =  Last_Min_X + range_W*0.45;
					text(loc_x,loc_y,plot_string_Frame,'color','black');
					if Key_Dynamic==1
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
			end

			
			% Plot cracks line if necessary
			if Key_PLOT(3,5) == 1
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
									%若结束角度大于起始角度,则直接绘制圆弧
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
									%若结束角度小于起始角度,则分成两部分绘制圆弧
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
									%若结束角度大于起始角度,则分成两部分绘制圆弧
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
									%若结束角度小于起始角度,则直接绘制圆弧
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
			
			% Plot shaped cracks if necessary.
			if Key_PLOT(3,5) == 2 & num_Crack(i)~=0
			disp(['      > Plotting shaped cracks......'])
				Plot_Shaped_Cracks(Shaped_Crack_Points);
			end
		
		
			% Save the current figure
			if j==1
				stress_xx(i_output) = getframe(gcf);       
				im=frame2im(stress_xx(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-xx';
				str2=Full_Pathname;
				FileName1 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName1,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName1,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			elseif j==2
				stress_yy(i_output) = getframe(gcf);       
				im=frame2im(stress_yy(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-yy';
				str2=Full_Pathname;
				FileName2 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName2,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName2,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			elseif j==3
				stress_xy(i_output) = getframe(gcf);       
				im=frame2im(stress_xy(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-xy';
				str2=Full_Pathname;
				FileName3 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName3,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName3,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			elseif j==4
				stress_vm(i_output) = getframe(gcf);       
				im=frame2im(stress_vm(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-vm';
				str2=Full_Pathname;
				FileName4 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName4,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName4,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			elseif j==5
				stress_S_1(i_output) = getframe(gcf);       
				im=frame2im(stress_S_1(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-1';
				str2=Full_Pathname;
				FileName4 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName4,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName4,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			elseif j==6
				stress_S_2(i_output) = getframe(gcf);       
				im=frame2im(stress_S_2(i_output));         
				[I,map]=rgb2ind(im,256);
				k=i_output-0;
				str1='_stress-2';
				str2=Full_Pathname;
				FileName4 =[str2,str1,'.gif'];
				if k==1;
					imwrite(I,map,FileName4,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
				else
					imwrite(I,map,FileName4,'gif','WriteMode','append','DelayTime',Time_Delay);
				end
			end
			close
		end
	end
end

clear Stress_Matrix



