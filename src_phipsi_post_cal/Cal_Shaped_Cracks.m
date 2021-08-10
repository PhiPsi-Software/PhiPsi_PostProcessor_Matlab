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

function [Shaped_Crack_Points]=Cal_Shaped_Cracks(Crack_X,Crack_Y,ifra,isub,num_Crack,Crack_Tip_Type,POS,...
                            Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
					        Coors_Vertex,Coors_Junction,Coors_Tip,DISP,scale)
% This function calculates the coodinates of points along the surface line of the shaped crack which are then used to plot the shaped cracks.

global aveg_area_ele
global Num_Ele_Bore Node_Coor Elem_Node Num_Elem
global Edge_crack_inter Full_Pathname
global Key_Data_Format

%计算增强单元的平均大小
Total_Enr_Area =0.0;
total_enr_Ele  =0;

if Key_Data_Format==1 
	Elem_Type = load([Full_Pathname,'.elty_',num2str(isub)]);
elseif Key_Data_Format==2  %Binary
	c_file = fopen([Full_Pathname,'.elty_',num2str(isub)],'rb');
	[cc_Elem_Type,cc_count]   = fread(c_file,inf,'int');
	fclose(c_file);
	%转换成Matlab中的数据格式
	Elem_Type = (reshape(cc_Elem_Type,num_Crack(isub),Num_Elem))';
end		

for i=1:Num_Elem
    N1  = Elem_Node(i,1);                                                  % Node 1 for current element
    N2  = Elem_Node(i,2);                                                  % Node 2 for current element
    N3  = Elem_Node(i,3);                                                  % Node 3 for current element
    N4  = Elem_Node(i,4);                                                  % Node 4 for current element
	NN = [N1 N2 N3 N4];                                                    % Four nodes of the current element  
	X_NODES = Node_Coor(NN,1);                                             % Coordinates of the four nodes of the current element
	Y_NODES = Node_Coor(NN,2);
	%如果是增强单元
	if sum(Elem_Type(i,:)) >= 1
	    total_enr_Ele = total_enr_Ele +1;
		Total_Enr_Area = Total_Enr_Area + polyarea(X_NODES,Y_NODES);
	end
end

aveg_area_ele_Enr = Total_Enr_Area/total_enr_Ele;

c_Num_Divis_Elment = 6;    % 11,6
offset_delta_Factor = 0.001;           %delta_L=offset_delta_Factor*sqrt(aveg_area_ele);
if isempty(Crack_X)==0
	for i = 1:num_Crack(ifra)
		%----------------------------------------------
		% Case 1: The crack has no junction tip, then:
		%----------------------------------------------
		if sum(Crack_Tip_Type(i,:)) <= 0  % Notice that if Crack_Tip_Type(iCrack,iTip)=1, then the tip in J tip.
			nPt = size(Crack_X{i},2);
			first_Tip = [];
			second_Tip= [];
			Offsetted_UP = [];
			Offsetted_DOWN = [];
			Shaped_Points =[];
			Last_Shaped_Points=[];
			Offs_Edge_Tip_Up=[];
			Offs_Edge_Tip_Down=[];
			Flag_Edge_Tip =0;
			tt=[];
			for iPt = 2:nPt
				x = [Crack_X{i}(iPt-1) Crack_X{i}(iPt)];
				y = [Crack_Y{i}(iPt-1) Crack_Y{i}(iPt)];
				% Length of the current segment.
				l_seg = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
				% Length of element. 
				l_elem= sqrt(aveg_area_ele_Enr);
				% Length of each division
				l_division = l_elem/c_Num_Divis_Elment;        % 5 can be changed!!!
				% Number of division.
				Num_Division = round(l_seg/l_division);
				if (Num_Division<=1)
				    Num_Division = 2;
				end
				% The value of the offset.
				offset_delta = offset_delta_Factor*l_elem;
				
				if iPt ==2 
				    if Crack_Tip_Type(i,1) ~=-2  % Not edge crack tip.
					    first_Tip = [x(1) y(1)];
					else                         % Edge crack tip.
					    Flag_Edge_Tip = 1;       % The first tip is the edge crack.
					    l_AB = [x(1) y(1); x(2) y(2)];
					    [Offs_Edge_Tip_Up,Offs_Edge_Tip_Down] = Cal_Offseted_Single_Point([x(1) y(1)],l_AB,offset_delta);
					end		    
				end
				if iPt ==nPt 
				    if Crack_Tip_Type(i,2) ~=-2  % Not edge crack tip.
					    second_Tip= [x(2) y(2)];
                    else                         % Edge crack tip.
					    Flag_Edge_Tip = 2;       % The second tip is the edge crack.
					    l_AB = [x(1) y(1); x(2) y(2)];
					    [Offs_Edge_Tip_Up,Offs_Edge_Tip_Down] = Cal_Offseted_Single_Point([x(2) y(2)],l_AB,offset_delta);				
					end
				end
				

				% Get the equal diversion points of line AB, then offset them by a small delta.
				[Div_Points,Offsetted_D_P_Up,Offsetted_D_P_Down] = Cal_Equal_Division_Points(Num_Division, ...
				                                                        [x(1) y(1); x(2) y(2)],offset_delta,0);								 
				Offsetted_UP   = [Offsetted_UP;   Offsetted_D_P_Up];
				Offsetted_DOWN = [Offsetted_DOWN; Offsetted_D_P_Down];
								
				tt = [tt;Div_Points];
			end
			
			
            Offsetted_UP  =unique(Offsetted_UP,'rows','stable');
			Offsetted_DOWN=unique(Offsetted_DOWN,'rows','stable');
			
			
			if Flag_Edge_Tip==0       % No edge crack.
			    Shaped_Points = [first_Tip;Offsetted_UP;second_Tip;flipud(Offsetted_DOWN)];
			elseif Flag_Edge_Tip ==1  % The first tip is the edge crack.
			    Shaped_Points = [Offs_Edge_Tip_Up;first_Tip;Offsetted_UP;second_Tip;flipud(Offsetted_DOWN);Offs_Edge_Tip_Down];
			elseif Flag_Edge_Tip ==2  % The second tip is the edge crack.
			    Shaped_Points = [first_Tip;Offsetted_UP;Offs_Edge_Tip_Up;second_Tip;Offs_Edge_Tip_Down;flipud(Offsetted_DOWN)];
			end
			
			% Remove duplicate rows of Shaped_Points.
            Shaped_Points=unique(Shaped_Points,'rows','stable');												  
			
			
			for j =1:size(Shaped_Points,1)
				% Get the local coordinates of the points of the crack. 
				[Kesi,Yita] = Cal_KesiYita_by_Coors(Shaped_Points(j,1),Shaped_Points(j,2));
				% Get the element number which contains the points of the crack. 
				[c_Elem_Num] = Cal_Ele_Num_by_Coors(Shaped_Points(j,1),Shaped_Points(j,2));
				% Calculate the displacement of the points of the crack.		
				[dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,ifra,DISP,Kesi,Yita...
													    ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
													    Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
                				
				Last_Shaped_Points(j,1) = Shaped_Points(j,1)+dis_x(j)*scale;
				Last_Shaped_Points(j,2) = Shaped_Points(j,2)+dis_y(j)*scale;
			end
			Shaped_Crack_Points{i}=Last_Shaped_Points;
		%----------------------------------------------
		% Case 2: The crack has junction tip, then:
		%---------------------------------------------- 
		else
		    % --------------------------------------
			% If the first tip is the junction tip.
			% --------------------------------------
			if Crack_Tip_Type(i,1)==1
				nPt = size(Crack_X{i},2);
				first_Tip_UP = [];first_Tip_DOWN=[];
				second_Tip2 = [];
				Offsetted_UP2 = [];
				Offsetted_DOWN2 = [];
				Shaped_Points2 =[];
				Last_Shaped_Points2=[];				
				for iPt = 2:nPt
					x = [Crack_X{i}(iPt-1) Crack_X{i}(iPt)];
					y = [Crack_Y{i}(iPt-1) Crack_Y{i}(iPt)];
					% Length of the current segment.
					l_seg = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
					% Length of element. 
					l_elem= sqrt(aveg_area_ele_Enr);
					% Length of each division.
					l_division = l_elem/c_Num_Divis_Elment;        % 5 can be changed!!!
					% Number of division of each segment.
					Num_Division = round(l_seg/l_division);
					if (Num_Division<=1)
						Num_Division = 2;
					end
					% The value of the offset.
				    offset_delta = offset_delta_Factor*l_elem;
					if iPt ==2
					    l_AB = [x(1) y(1); x(2) y(2)];
						[first_Tip_UP,first_Tip_DOWN] = ...
						                  Cal_Offseted_Single_Point([x(1) y(1)],l_AB,1.0*offset_delta);
										  
					end
					if iPt ==nPt
						second_Tip2= [x(2) y(2)];
					end
					% Length of the current crack segment, i.e. the value of the offset.
					% offset_delta = 1e-3 * sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
					
					% Get the equal diversion points of line AB, then offset them by a small delta.
					[Div_Points,Offsetted_D_P_Up,Offsetted_D_P_Down] = Cal_Equal_Division_Points(Num_Division, ...
																[x(1) y(1); x(2) y(2)],offset_delta,0);
					Offsetted_UP2   = [Offsetted_UP2;   Offsetted_D_P_Up];
					Offsetted_DOWN2 = [Offsetted_DOWN2; Offsetted_D_P_Down];
				end
				
				% Get the first point of Offsetted_UP2, name as Fir_Poit_Offsetted_UP2.
				Fir_Poit_Offsetted_UP2 = Offsetted_UP2(1,:);
				% Calculate signed distance of Fir_Poit_Offsetted_UP2 to l_AB.
				[Signed_Distance_1] = Cal_Signed_Distance(l_AB,Fir_Poit_Offsetted_UP2);
			    % Calculate signed distance of first_Tip_UP to l_AB.
				[Signed_Distance_2] = Cal_Signed_Distance(l_AB,first_Tip_UP);
				% If first_Tip_UP and Fir_Poit_Offsetted_UP2 are at the same side of l_AB, then:
				if Signed_Distance_1*Signed_Distance_2 >0.0
				    Shaped_Points2 = [first_Tip_DOWN;first_Tip_UP;Offsetted_UP2;second_Tip2;flipud(Offsetted_DOWN2)];	
				else
                    Shaped_Points2 = [first_Tip_UP;first_Tip_DOWN;Offsetted_UP2;second_Tip2;flipud(Offsetted_DOWN2)];	
				end
				% Remove duplicate rows of Shaped_Points2.
				Shaped_Points2=unique(Shaped_Points2,'rows','stable');
				for j =1:size(Shaped_Points2,1)
					% Get the local coordinates of the points of the crack. 
					[Kesi,Yita] = Cal_KesiYita_by_Coors(Shaped_Points2(j,1),Shaped_Points2(j,2));
					% Get the element number which contains the points of the crack. 
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(Shaped_Points2(j,1),Shaped_Points2(j,2));
					% Calculate the displacement of the points of the crack.		
					[dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,ifra,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
					Last_Shaped_Points2(j,1) = Shaped_Points2(j,1)+dis_x(j)*scale;
					Last_Shaped_Points2(j,2) = Shaped_Points2(j,2)+dis_y(j)*scale;
				end
				
				Shaped_Crack_Points{i}=Last_Shaped_Points2;		
				
			end
			% --------------------------------------
			% If the second tip is the junction tip.
			% --------------------------------------
			if Crack_Tip_Type(i,2)==1
				nPt = size(Crack_X{i},2);
				first_Tip3 = [];
				second_Tip_UP = [];second_Tip_DOWN = [];
				Offsetted_UP3 = [];
				Offsetted_DOWN3 = [];
				Shaped_Points3 =[];
				Last_Shaped_Points3=[];	
				for iPt = 2:nPt
					x = [Crack_X{i}(iPt-1) Crack_X{i}(iPt)];
					y = [Crack_Y{i}(iPt-1) Crack_Y{i}(iPt)];
					% Length of the current segment.
					l_seg = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
					% Length of element. 
					l_elem= sqrt(aveg_area_ele_Enr);
					% Length of each division.
					l_division = l_elem/c_Num_Divis_Elment;     % 5 can be changed!!!
					% Number of division of each segment.
					Num_Division = round(l_seg/l_division);
					if (Num_Division<=1)
						Num_Division = 2;
					end
					% The value of the offset.
				    offset_delta = offset_delta_Factor*l_elem;
					
					if iPt ==2
						first_Tip3 = [x(1) y(1)];
					end
					if iPt ==nPt
					    l_AB = [x(1) y(1); x(2) y(2)];
						[second_Tip_UP,second_Tip_DOWN] = ...
						                  Cal_Offseted_Single_Point([x(2) y(2)],l_AB,offset_delta);
					end

					% Get the equal diversion points of line AB, then offset them by a small delta.
					[Div_Points,Offsetted_D_P_Up,Offsetted_D_P_Down] = Cal_Equal_Division_Points(Num_Division, ...
																      [x(1) y(1); x(2) y(2)],offset_delta,0);
																			 
					Offsetted_UP3   = [Offsetted_UP3;   Offsetted_D_P_Up];
					Offsetted_DOWN3 = [Offsetted_DOWN3; Offsetted_D_P_Down];
				end
				% Get the last point of Offsetted_UP3, name as Last_Poit_Offsetted_UP3.
				Last_Poit_Offsetted_UP3 = Offsetted_UP3(size(Offsetted_UP3,1),:);
				% Calculate signed distance of Last_Poit_Offsetted_UP3 to l_AB.
				[Signed_Distance_1] = Cal_Signed_Distance(l_AB,Last_Poit_Offsetted_UP3);
			    % Calculate signed distance of second_Tip_UP to l_AB.
				[Signed_Distance_2] = Cal_Signed_Distance(l_AB,second_Tip_UP);
				% If first_Tip_UP and Fir_Poit_Offsetted_UP2 are at the same side of l_AB, then:
				if Signed_Distance_1*Signed_Distance_2 >0
				    Shaped_Points3 = [first_Tip3;Offsetted_UP3;second_Tip_UP;second_Tip_DOWN;flipud(Offsetted_DOWN3)];
				else
                    Shaped_Points3 = [first_Tip3;Offsetted_UP3;second_Tip_DOWN;second_Tip_UP;flipud(Offsetted_DOWN3)];		
				end
				% Remove duplicate rows of Shaped_Points3.
				Shaped_Points3=unique(Shaped_Points3,'rows','stable');
				for j =1:size(Shaped_Points3,1)
					% Get the local coordinates of the points of the crack. 
					[Kesi,Yita] = Cal_KesiYita_by_Coors(Shaped_Points3(j,1),Shaped_Points3(j,2));
					% Get the element number which contains the points of the crack. 
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(Shaped_Points3(j,1),Shaped_Points3(j,2));
					% Calculate the displacement of the points of the crack.		
					[dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,ifra,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
					Last_Shaped_Points3(j,1) = Shaped_Points3(j,1)+dis_x(j)*scale;
					Last_Shaped_Points3(j,2) = Shaped_Points3(j,2)+dis_y(j)*scale;
				end
				Shaped_Crack_Points{i}=Last_Shaped_Points3;
			end				
		end
	end	
end
