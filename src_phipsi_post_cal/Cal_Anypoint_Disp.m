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

function [Disp_x,Disp_y] = Cal_Anypoint_Disp(iElem,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita,...
					                         Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
					                         Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y)
% This function calculates displacement of a given point.											 

global Node_Coor Elem_Node Elem_Material
global Num_Elem num_Crack
global Material_Para Key_Problem
global G_NN G_X_NODES G_Y_NODES G_X_Min G_X_Max G_Y_Min G_Y_Max
global Key_TipEnrich
global num_Hole Hole_Coor Enriched_Node_Type_Hl POS_Hl Elem_Type_Hl
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl POS_Incl Elem_Type_Incl
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Key_Heaviside_Value 
global Key_Hole_Value 
global Cross_Point_RABCD
global num_Cross Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Ellipse_Hole_Coor num_Ellipse_Hole
%global Node_Cross_elem

% Calculate the displacement of the points of the crack. 
n1  = Elem_Node(iElem,1);                                                  
n2  = Elem_Node(iElem,2);                                                  
n3  = Elem_Node(iElem,3);                                                  
n4  = Elem_Node(iElem,4);

% Material number of the current element.
mat_num    = Elem_Material(iElem);
% Material type of the current element.
c_mat_type = 1;

% Four nodes of the current element
NODES_iElem = G_NN(:,iElem)';
                                           
U = [DISP(n1,2) DISP(n1,3) DISP(n2,2) DISP(n2,3)...
	 DISP(n3,2) DISP(n3,3) DISP(n4,2) DISP(n4,3)];

% Coordinates of the four nodes of the current element
X_NODES = G_X_NODES(:,iElem);
Y_NODES = G_Y_NODES(:,iElem);

% Calculates N, dNdkesi, J and the determinant of Jacobian matrix.
[N,dNdkesi,J,detJ]  = Cal_N_dNdkesi_J_detJ(Kesi,Yita,X_NODES,Y_NODES);

% The inverse matrix of J.
% Inverse_J = inv(J);

% Calculate dNdx.
% dNdx = dNdkesi * Inverse_J;
% size(dNdx)
% Disp_x = U(1)*N(1,1) + U(3)*N(1,3) + U(5)*N(1,5) + U(7)*N(1,7);
% Disp_y = U(2)*N(1,1) + U(4)*N(1,3) + U(6)*N(1,5) + U(8)*N(1,7);  
Disp_x = U(1:2:7)*N(1,1:2:7)';
Disp_y = U(2:2:8)*N(1,1:2:7)';

tem_N(1:4) = N(1,1:2:7);

% Global coordinates of the gauss point.
% Global_coor_Gauss(1) = N(1,1)*X_NODES(1)+N(1,3)*X_NODES(2)+N(1,5)*X_NODES(3)+N(1,7)*X_NODES(4);
% Global_coor_Gauss(2) = N(1,1)*Y_NODES(1)+N(1,3)*Y_NODES(2)+N(1,5)*Y_NODES(3)+N(1,7)*Y_NODES(4);
Global_coor_Gauss(1) = N(1,1:2:7)*X_NODES(1:4);
Global_coor_Gauss(2) = N(1,1:2:7)*Y_NODES(1:4);
c_G_x = Global_coor_Gauss(1);
c_G_y = Global_coor_Gauss(2);

ttt=0;
% *****************************
% Loop through each crack.
% *****************************
for iCrack=1:num_Crack(isub)									
    % Conventional element.
	if any(Enriched_Node_Type(NODES_iElem,iCrack)) ==0
	    Disp_x = Disp_x;
	    Disp_y = Disp_y;
	% Enriched elements.
	elseif any(Enriched_Node_Type(NODES_iElem,iCrack)) ~=0
		% Loop through each node.
		for iNode = 1:4
			% Initialize flag of kinked crack of the tip segment.
			Flag_Check_Kink_Tip = 0;
			Yes_Gauss_in_BCD    = 0;
			%-----------------------------------------------------------
			%----------- Case 1: Heaviside enriched node ---------------
			%-----------------------------------------------------------
			if Enriched_Node_Type(NODES_iElem(iNode),iCrack) ==2        % Heaviside enriched node
			    % Get the number of the current enriched node.
				Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack);
				% if 2:fully cracked element without kink point or 3:fully cracked element with kink point, then:
				if Elem_Type(iElem,iCrack) == 2 || Elem_Type(iElem,iCrack) == 3
					ref_elem = iElem;
					% Coordinates of the 2 intersections(A,B) of the crack segment and these 4 element(fully cracked) lines.
					Coor_AB = [Coors_Element_Crack(ref_elem,iCrack,1) Coors_Element_Crack(ref_elem,iCrack,2);
							   Coors_Element_Crack(ref_elem,iCrack,3) Coors_Element_Crack(ref_elem,iCrack,4)];
					% if 3:fully cracked element with kink point, then:
					if Elem_Type(iElem,iCrack) == 3
						% Calculates the distance of the kink point(K) to the line(AB) composed of the 2 intersections(A,B).
						distance_Vertex = Cal_Signed_Distance(Coor_AB,Coors_Vertex(ref_elem,:));	
						H_Vertex = sign(distance_Vertex);	
						% Calculates the distance of the gauss point to the line(AB) composed of the 2 intersections(A,B).
						distance_Gauss = Cal_Signed_Distance(Coor_AB,Global_coor_Gauss);	
						H_gp = sign(distance_Gauss);	
						if Key_Heaviside_Value ==0
						    if H_gp <= 0.0
							    H_gp = 0.0;
							end
						end
						% Gauss point is under the triangle KAB. 
						if H_Vertex*H_gp <= 0
							H_gp = H_gp;
						else
							% The triangle KAB.
							Tri_KAB = [Coor_AB(1,:); Coor_AB(2,:);Coors_Vertex(ref_elem,:)];
							% Judge if the gauss point is inside the triangle KAB.
							Yes_in_KAB = inpolygon(Global_coor_Gauss(1),Global_coor_Gauss(2), ...
												   Tri_KAB(:,1),Tri_KAB(:,2));
							% If the gauss point is inside the triangle KAB, then:
							if Yes_in_KAB == 1
								H_gp = -H_gp;
							else
								H_gp =  H_gp;
							end
						end
						% Calculates the distance of the node to the line(AB).
						distance_Node = Cal_Signed_Distance(Coor_AB,[X_NODES(iNode),Y_NODES(iNode)]);
						H_i = sign(distance_Node);
						if Key_Heaviside_Value ==0
						    if H_i <= 0.0
							    H_i = 0.0;
							end
						end
					% if 2:fully cracked element without kink point, then:
					else
						distance_Gauss = Cal_Signed_Distance(Coor_AB,Global_coor_Gauss);
						H_gp = sign(distance_Gauss);	
						if Key_Heaviside_Value ==0
						    if H_gp <= 0.0
							    H_gp = 0.0;
							end
						end
						distance_Node = Cal_Signed_Distance(Coor_AB,[X_NODES(iNode),Y_NODES(iNode)]);
						H_i = sign(distance_Node);
						if Key_Heaviside_Value ==0
						    if H_i <= 0.0
							    H_i = 0.0;
							end
						end
					end
					% Calculate displacements.
					Disp_x = Disp_x + (H_gp-H_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
					Disp_y = Disp_y + (H_gp-H_i)*tem_N(iNode)*DISP(Num_Enriched_Node,3);
					ttt=ttt+(H_gp-H_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
					
				else
				    Disp_x = Disp_x;
					Disp_y = Disp_y;
				end

			%-----------------------------------------------------------
			%----------- Case 2: Junction enriched node ----------------
			%-----------------------------------------------------------	
			elseif Enriched_Node_Type(NODES_iElem(iNode),iCrack) ==3    % Junction enriched node
				% Get the number of the current enriched node.
				Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack);
		        % Get the number of the Junction element.
				Jun_elem = Node_Jun_elem(NODES_iElem(iNode),iCrack);
				% Coordinates of the 2 intersections(A,B) of the main crack segment and these 4 element(fully cracked) lines.
				Coor_AB = [Coors_Element_Crack(Jun_elem,iCrack,1) Coors_Element_Crack(Jun_elem,iCrack,2);
						   Coors_Element_Crack(Jun_elem,iCrack,3) Coors_Element_Crack(Jun_elem,iCrack,4)];
				% Coordinates of the 2 intersections(C,D) of the minor crack segment and these 4 element(fully cracked) lines.
				Coor_CD = [Coors_Junction(Jun_elem,iCrack,1) Coors_Junction(Jun_elem,iCrack,2);
						   Coors_Junction(Jun_elem,iCrack,3) Coors_Junction(Jun_elem,iCrack,4)];
				x0 = [Coors_Junction(Jun_elem,iCrack,1) Coors_Junction(Jun_elem,iCrack,2)];
				% Calculates the distance of nodes and the gauss point to the main crack segment.
				distance_Node_M  = Cal_Signed_Distance(Coor_AB,[X_NODES(iNode),Y_NODES(iNode)]);
				distance_Gauss_M = Cal_Signed_Distance(Coor_AB,Global_coor_Gauss);
				% Calculates the distance of nodes and the gauss point to the minor crack segment.
				distance_Node_m  = Cal_Signed_Distance(Coor_CD,[X_NODES(iNode),Y_NODES(iNode)]);
				distance_Gauss_m = Cal_Signed_Distance(Coor_CD,Global_coor_Gauss);
				% Calculates the distance of x0 to the main crack segment.
				distance_x0  = Cal_Signed_Distance(Coor_AB,x0);
				H_gp_M = sign(distance_Gauss_M);
				H_i_M  = sign(distance_Node_M); 
				H_gp_m = sign(distance_Gauss_m);
				H_i_m  = sign(distance_Node_m); 
				H_x0   = sign(distance_x0);
				if H_gp_M*H_x0 > 0
					if H_gp_m > 0
						J_gp=  1;
					else
						J_gp= -1;		
					end
				else
					J_gp=  0;
				end
				if H_i_M*H_x0 > 0
					if H_i_m > 0
						J_i=  1;
					else
						J_i= -1;
					end
				else
					J_i=  0;
				end
				% Calculate displacements.
				Disp_x = Disp_x + (J_gp-J_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
				Disp_y = Disp_y + (J_gp-J_i)*tem_N(iNode)*DISP(Num_Enriched_Node,3);
			%-----------------------------------------------------------------------------
			%----------- Case 3: Junction enriched node of crack and hole ----------------
			%----------- Added on 2017-05-19                              ----------------
			%-----------------------------------------------------------------------------	
			elseif Enriched_Node_Type(NODES_iElem(iNode),iCrack) ==6    % Junction enriched node of crack and hole
				% Get the number of the current enriched node.
				Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack);
		        % Get the number of the Junction element.
				Jun_elem = Node_Jun_elem(NODES_iElem(iNode),iCrack);
				% Hole number.
                c_Hole   = Node_Jun_Hole(NODES_iElem(iNode),iCrack); 
				if num_Hole >=1
					x0_Hole =  Hole_Coor(c_Hole,1);
					y0_Hole =  Hole_Coor(c_Hole,2);
					R0_Hole =  Hole_Coor(c_Hole,3);
				end
				if num_Ellipse_Hole >=1
					x0_Hole =  Ellipse_Hole_Coor(c_Hole,1);
					y0_Hole =  Ellipse_Hole_Coor(c_Hole,2);
					R0_Hole =  Ellipse_Hole_Coor(c_Hole,3);
				end
				% Coordinates of the 2 intersections(C,D) of the minor crack segment and these 4 element(fully cracked) lines.
				Coor_CD = [Coors_Junction(Jun_elem,iCrack,1) Coors_Junction(Jun_elem,iCrack,2);
						   Coors_Junction(Jun_elem,iCrack,3) Coors_Junction(Jun_elem,iCrack,4)];
				x0 = [Coors_Junction(Jun_elem,iCrack,1) Coors_Junction(Jun_elem,iCrack,2)];
				
                % Calculates the distance of nodes and the gauss point to the hole.
                distance_Node_M  = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,[X_NODES(iNode),Y_NODES(iNode)]);
                distance_Gauss_M = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,Global_coor_Gauss);
				
				
				% Calculates the distance of nodes and the gauss point to the minor crack segment.
				distance_Node_m  = Cal_Signed_Distance(Coor_CD,[X_NODES(iNode),Y_NODES(iNode)]);
				distance_Gauss_m = Cal_Signed_Distance(Coor_CD,Global_coor_Gauss);
				
				
				% Calculates the distance of x0 to the main crack segment.
				distance_x0  = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,x0);
				
				
				H_gp_M = sign(distance_Gauss_M);
				H_i_M  = sign(distance_Node_M); 
				H_gp_m = sign(distance_Gauss_m);
				H_i_m  = sign(distance_Node_m); 
				H_x0   = sign(distance_x0);
				
				%...............
				%类似Heaviside
				%...............
				% if H_gp_m > 0
					% J_gp=  1;
				% else
					% J_gp= -1;		
				% end
				% if H_i_m > 0
					% J_i=  1;
				% else
					% J_i= -1;
				% end
				
				%.................
				%采用Junction增强
				%.................
				if H_gp_M*H_x0 > 0
					if H_gp_m > 0
						J_gp=  1;
					else
						J_gp= -1;		
					end
				else
					J_gp=  0;
				end
				if H_i_M*H_x0 > 0
					if H_i_m > 0
						J_i=  1;
					else
						J_i= -1;
					end
				else
					J_i=  0;
				end
				
				% Calculate displacements.
				Disp_x = Disp_x + (J_gp-J_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
				Disp_y = Disp_y + (J_gp-J_i)*tem_N(iNode)*DISP(Num_Enriched_Node,3);
	        %-----------------------------------------------------------
			%-------- Case 4: Cross enriched node (2017-05-04) ---------
			%-----------------------------------------------------------	
			elseif Enriched_Node_Type(NODES_iElem(iNode),iCrack) ==4    % Cross enriched node
				% give up this scheme		
				

			%-----------------------------------------------------------
			%----------- Case 5: Tip enriched node ---------------------
			%-----------------------------------------------------------
			elseif Enriched_Node_Type(NODES_iElem(iNode),iCrack) ==1    % Tip enriched node
			    ref_elem = 0;
				% Find the fully tip enriched element which the current tip enriched node attached to. 
				% Fully tip enriched element.
				if Elem_Type(iElem,iCrack) == 1
					ref_elem = iElem;
					ref_elem_1=ref_elem;
				% Partly tip enriched element.
				else
					[find_element,xx] = find(Elem_Node == NODES_iElem(iNode));
					[find_location,xx] = find(Elem_Type(find_element,:) ==1);
					if isempty(find_location) ~=1
						ref_elem = find_element(find_location);
						Flag_Check_Kink_Tip =1;
					else
						[find_element_1,xx1] = find(Enriched_Node_Type(NODES_iElem(find_element,:),iCrack) == 1);
						node_p = [];
						for i =1:size(find_element_1,1);
							node_p = [node_p; NODES_iElem(find_element_1,xx1(i))];
						end
						for i =1:size(node_p,1)
							[find_element,xx] = find(Elem_Node == node_p(i));
							[find_location,xx] = find(Elem_Type(find_element,:) ==1);
							ref_elem = find_element(find_location);
						end
					end
					ref_elem_2=ref_elem;
				end
				% Coordinates of the 2 intersections(A,B) of the main crack segment and these 4 element(fully cracked) lines
				% and the coordinates of the tip.
				% Coors_Element_Crack(47,:)
				% Coors_Element_Crack(93,:)
				Coor_AB = [Coors_Element_Crack(ref_elem,iCrack,1) Coors_Element_Crack(ref_elem,iCrack,2);
						   Coors_Element_Crack(ref_elem,iCrack,3) Coors_Element_Crack(ref_elem,iCrack,4)];
				Segment     = Coor_AB(2,:) - Coor_AB(1,:);
				% The angle of the crack tip segment.
				omega   = atan2(Segment(2),Segment(1));
				% Coordinates of the tip.
				Coor_Tip = [Coor_AB(2,1) Coor_AB(2,2)];
				% The transformation matrix from the local tip coordinates to the global coordinates. 
				Trans_Matrix = [cos(omega),sin(omega);
							   -sin(omega),cos(omega)];
                %##############################################################################
				% The following part is very important.
				% For more information, see my doctoral dissertation on page 19.
				% num_points_crack = size(Crack(iCrack).coor,1);Crack_X{i} 
				num_points_crack = size(Crack_X{iCrack},2);
				%裂尖�?��单元�?个节�?
				c_X_NODES = G_X_NODES(:,ref_elem);
                c_Y_NODES = G_Y_NODES(:,ref_elem);
				if Flag_Check_Kink_Tip == 1 && num_points_crack > 2
					% if element is a fully cracked element with or without kink point.
					if Elem_Type(iElem,iCrack) == 2 || Elem_Type(iElem,iCrack) == 3
						coor_A = Coor_Tip;          % About point A and point B, see my doctoral dissertation.
						% 获得当前单元内的裂尖坐标.
						for iii=1:2
						    if iii==1
						        ttt_x = Crack_X{iCrack}(1);
								ttt_y = Crack_Y{iCrack}(1);
								% max(c_X_NODES(1:4))
								% min(c_X_NODES(1:4)) 
								% if ttt_x < max(c_X_NODES(1:4)) & ttt_x > min(c_X_NODES(1:4))
									% if ttt_y < max(c_Y_NODES(1:4)) & ttt_y > min(c_Y_NODES(1:4))
									    tip_x = ttt_x;
										tip_y = ttt_y;
								    % end
								% end
							elseif iii==2
						        ttt_x = Crack_X{iCrack}(num_points_crack);
								ttt_y = Crack_Y{iCrack}(num_points_crack);
								% if ttt_x < max(c_X_NODES(1:4)) & ttt_x > min(c_X_NODES(1:4))
									% if ttt_y < max(c_Y_NODES(1:4)) & ttt_y > min(c_Y_NODES(1:4))
									    tip_x = ttt_x;
										tip_y = ttt_y;
								    % end
								% end
							end
							
						end
						% tip_x = x_cr_tip_nodes(iCrack,NODES_iElem(iNode));
						% tip_y = y_cr_tip_nodes(iCrack,NODES_iElem(iNode));
						
						if Crack_X{iCrack}(1) == tip_x && Crack_Y{iCrack}(1) == tip_y	
							coor_B = [Crack_X{iCrack}(2),Crack_Y{iCrack}(2)];
						elseif Crack_X{iCrack}(num_points_crack) == tip_x && Crack_Y{iCrack}(num_points_crack)== tip_y	
							coor_B = [Crack_X{iCrack}(num_points_crack-1),Crack_Y{iCrack}(num_points_crack-1)];
						end
						
						len_AB = norm([coor_A(1)-coor_B(1),coor_A(2)-coor_B(2)]);
						
						% iCrack
						
						% coor_A(1), Crack_X{iCrack}(1)
						% coor_A(2), Crack_Y{iCrack}(1)
						% coor_A(1), Crack_X{iCrack}(num_points_crack)
						% coor_A(2), Crack_Y{iCrack}(num_points_crack)
						

						% if  coor_A(1) == Crack_X{iCrack}(1) && coor_A(2) == Crack_Y{iCrack}(1)
						if  (abs(coor_A(1) - Crack_X{iCrack}(1))<=1.0e-10) && (abs(coor_A(2) - Crack_Y{iCrack}(1))<=1.0e-10)
							coor_C = [Crack_X{iCrack}(3),Crack_Y{iCrack}(3)];
						elseif abs(coor_A(1) - Crack_X{iCrack}(num_points_crack))<=1.0e-10 & ...
							   abs(coor_A(2) - Crack_Y{iCrack}(num_points_crack))<=1.0e-10
							coor_C = [Crack_X{iCrack}(num_points_crack-2),Crack_Y{iCrack}(num_points_crack-2)];
						end
						
						
						% D is an added point on the extended line of AB from A to B.
						% D is added t form the triangle BCD.
						[~,coor_D] = Cal_Shorten_or_Extend_Line([coor_A;coor_B],2*len_AB,'B');
						Tri_BCD = [coor_B;coor_C;coor_D];
						Yes_Gauss_in_BCD = inpolygon(Global_coor_Gauss(1),Global_coor_Gauss(2),Tri_BCD(:,1),Tri_BCD(:,2));
						% If the current gauss point is in the triangle BCD, then:
						if Yes_Gauss_in_BCD==1
							[angle_AB_BC] = Cal_Angle_of_AB_and_BC([coor_A;coor_B],[coor_B;coor_C]);
							% ngle_AB_BC_360 = angle_AB_BC*180/pi
							[angle_AB_BGauss] = Cal_Angle_of_AB_and_BC([coor_A;coor_B],[coor_B;Global_coor_Gauss]);
							% angle_AB_BGauss_360 = angle_AB_BGauss*180/pi
							len_BGauss  = norm([coor_B(1)-Global_coor_Gauss(1),coor_B(2)-Global_coor_Gauss(2)]);
							if angle_AB_BGauss >= angle_AB_BC  % See equation (2-20) of my doctoral dissertation.  
								thete_Doc   = pi/2*(angle_AB_BGauss-angle_AB_BC)/(1.5*pi-angle_AB_BC);
							else
								thete_Doc   = pi/2*(angle_AB_BGauss-angle_AB_BC)/(angle_AB_BC-0.5*pi);
							end
							Global_coor_Gauss_Doc = [-len_AB-len_BGauss*cos(thete_Doc),-len_BGauss*sin(thete_Doc)];
							% Global_coor_Gauss     = [-len_AB-len_BGauss*cos(thete_Doc),-len_BGauss*sin(thete_Doc)];
							% iGauss
						end
					end			
				end
				%##############################################################################
				if Yes_Gauss_in_BCD==0
					xp_Gauss    = Trans_Matrix*(Global_coor_Gauss-Coor_Tip)';
					r_Gauss     = sqrt(xp_Gauss(1)^2 + xp_Gauss(2)^2);
					theta_Gauss = atan2(xp_Gauss(2),xp_Gauss(1));			
				elseif Yes_Gauss_in_BCD==1
					r_Gauss     = sqrt(Global_coor_Gauss_Doc(1)^2 + Global_coor_Gauss_Doc(2)^2);
					theta_Gauss = atan2(Global_coor_Gauss_Doc(2),Global_coor_Gauss_Doc(1));
				end
			    
				% Calculate the value(r,theta) of the gauss point in the local tip coordinates.
				% xp_Gauss = Trans_Matrix*(Global_coor_Gauss-Coor_Tip)';
				% r_Gauss  =    sqrt(xp_Gauss(1)^2 + xp_Gauss(2)^2);
				% theta_Gauss = atan2(xp_Gauss(2),xp_Gauss(1));
				% Check theta.
				if (theta_Gauss >pi  | theta_Gauss < -pi)
					disp('    Error :: Angle is wrong when calculates r and theta!') 
					Error_Message
				end
				
				% Calculate mu of orthotropic material.
				if c_mat_type ==1   % ISO material.
					mu =[];
				elseif c_mat_type ==2 |  c_mat_type ==3
					[mu,p1,p2,q1,q2] = Cal_Orth_mu(omega,mat_num);
				end
			
				% Calculate the tip enrichment functions F and their derivative, dFdx and dFdy.
				[F_Gauss,~,~] = Cal_F_dFdx_dFdy(r_Gauss,theta_Gauss,omega,mu,c_mat_type);
				
				% Calculate the value(r,theta) of node in the local tip coordinates.
				xp_Node = Trans_Matrix*([X_NODES(iNode),Y_NODES(iNode)]-Coor_Tip)';
				r_Node  =    sqrt(xp_Node(1)^2 + xp_Node(2)^2);
				theta_Node = atan2(xp_Node(2),xp_Node(1)); 
				[F_Node,~,~] = Cal_F_dFdx_dFdy(r_Node,theta_Node,omega,mu,c_mat_type);
				
				BI_enr = [];
				% for i_F =1:4
					% Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack)+i_F-1;
                    % F = F_Gauss(i_F)-F_Node(i_F);			
					% Disp_x = Disp_x + F*tem_N(iNode)*DISP(Num_Enriched_Node,2);
					% Disp_y = Disp_y + F*tem_N(iNode)*DISP(Num_Enriched_Node,3);
				% end
				
				
				% Calculate displacements for conventional crack.
				if Key_TipEnrich ==1
					for i_F =1:4
						% Get the number of the current enriched node.
						Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack)+i_F-1;
						F = F_Gauss(i_F)-F_Node(i_F);		
						Disp_x = Disp_x + F*tem_N(iNode)*DISP(Num_Enriched_Node,2);
						Disp_y = Disp_y + F*tem_N(iNode)*DISP(Num_Enriched_Node,3);
					end
				%裂尖增强方案2和方�?						
				elseif Key_TipEnrich ==2 | Key_TipEnrich ==3| Key_TipEnrich ==4
					% Get the number of the current enriched node.
					Num_Enriched_Node = POS(NODES_iElem(iNode),iCrack);
					F = F_Gauss(1)-F_Node(1);		
					Disp_x = Disp_x + F*tem_N(iNode)*DISP(Num_Enriched_Node,2);
					Disp_y = Disp_y + F*tem_N(iNode)*DISP(Num_Enriched_Node,3);
				end	
					

				
			end
		end
	end
end
% ttt
% ttt
% *****************************
% Loop through each hole.
% *****************************
if num_Hole ~=0
	for iHole=1:num_Hole									
		% Conventional element.
		if any(Enriched_Node_Type_Hl(NODES_iElem,iHole)) ==0
			Disp_x = Disp_x;
			Disp_y = Disp_y;
		% Enriched elements.
		elseif any(Enriched_Node_Type_Hl(NODES_iElem,iHole)) ~=0
            c_Hole_x = Hole_Coor(iHole,1);
            c_Hole_y = Hole_Coor(iHole,2);
            c_Hole_r = Hole_Coor(iHole,3);
			%圆上的点Hl_gp显然等于1
            Hl_gp = 1.0;
			% Loop through each node.
			for iNode = 1:4
				% ----------------------------------------------
				% ----------- Hole enriched node ---------------
				% ----------------------------------------------
				if Enriched_Node_Type_Hl(NODES_iElem(iNode),iHole) ==1  
                    if Elem_Type_Hl(iElem,iHole) == 1			
						% Get the number of the current enriched node.
						Num_Enriched_Node = POS_Hl(NODES_iElem(iNode),iHole);
						% iHole
						c_Dis_N = sqrt((X_NODES(iNode)-c_Hole_x)^2+(Y_NODES(iNode)-c_Hole_y)^2);
						if c_Dis_N<c_Hole_r
							if Key_Hole_Value ==0
							    Hl_i = 0.0;
							elseif Key_Hole_Value ==-1
							    Hl_i = -1.0;
							end
						else
							Hl_i = 1.0;
						end
						Disp_x = Disp_x + (Hl_gp-Hl_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
						Disp_y = Disp_y + (Hl_gp-Hl_i)*tem_N(iNode)*DISP(Num_Enriched_Node,3);
					end
				else
				    Disp_x = Disp_x;
					Disp_y = Disp_y;
				end
			end
		end
	end
end

% *****************************
% Loop through each cross.
% *****************************
if num_Cross ~=0
	for i_Cross=1:num_Cross									
		% Conventional element.
		if any(Enriched_Node_Type_Cross(NODES_iElem,i_Cross)) ==0
			%Disp_x = Disp_x;
			%Disp_y = Disp_y;
		% Enriched elements.
		elseif any(Enriched_Node_Type_Cross(NODES_iElem,i_Cross)) ~=0
			% Loop through each node.
			for iNode = 1:4
				% ----------------------------------------------
				% ----------- Hole enriched node ---------------
				% ----------------------------------------------
				if Enriched_Node_Type_Cross(NODES_iElem(iNode),i_Cross) ==1  		
					 % Get the number of the current enriched node.
					 Num_Enriched_Node = POS_Cross(NODES_iElem(iNode),i_Cross);
                     c_Cross_elem      = Node_Cross_elem(NODES_iElem(iNode),i_Cross);
                     % print *,Num_Enri_Node
		             if Elem_Type_Cross(c_Cross_elem,i_Cross) ==1
                         %Coordinates of the 2 intersections(A,B) of the main crack segment and these 4 element(fully cracked) lines.
                          Coor_AB(1,1:2)=Cross_Point_RABCD(i_Cross,2,1:2);
                          Coor_AB(2,1:2)=Cross_Point_RABCD(i_Cross,3,1:2);
                          %Coordinates of the 2 intersections(C,D) of the minor crack segment and these 4 element(fully cracked) lines.
                          Coor_CD(1,1:2)=Cross_Point_RABCD(i_Cross,4,1:2);
                          Coor_CD(2,1:2)=Cross_Point_RABCD(i_Cross,5,1:2);
				
                          %Calculates the distance of nodes and the gauss point to the main crack segment.
                          distance_Node_M  = Cal_Signed_Distance(Coor_AB,[X_NODES(iNode),Y_NODES(iNode)]);   
                          distance_Gauss_M = Cal_Signed_Distance(Coor_AB,Global_coor_Gauss);                  
                          %Calculates the distance of nodes and the gauss point to the minor crack segment.
                          distance_Node_sm = Cal_Signed_Distance(Coor_CD,[X_NODES(iNode),Y_NODES(iNode)]);
                          distance_Gauss_sm= Cal_Signed_Distance(Coor_CD,Global_coor_Gauss);
						  
						  
						  if Key_Heaviside_Value==-1
                              Cross_gp         = sign(distance_Gauss_M*distance_Gauss_sm);
                              Cross_i          = sign(distance_Node_M*distance_Node_sm);
						  elseif Key_Heaviside_Value==0
                              Cross_gp         = Cal_Sign_1_and_0(distance_Gauss_M*distance_Gauss_sm);
                              Cross_i          = Cal_Sign_1_and_0(distance_Node_M*distance_Node_sm);
						  end
						  
						  
					%测试(裂缝和Hole相遇,增加一个Junction增强节点) 
                      Coor_CD(1,1:2) = [10.00e-3,10.500e-3];
                      Coor_CD(2,1:2) = [11.00e-3,10.500e-3];
                     x0_Hole =  Hole_Coor(1,1);
                     y0_Hole =  Hole_Coor(1,2);
                     R0_Hole =  Hole_Coor(1,3);
				  
                      x0 = [11.000e-3,10.500e-3];
					  
					distance_Node_M  = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,[X_NODES(iNode),Y_NODES(iNode)]);
					distance_Gauss_M = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,Global_coor_Gauss);
					distance_Node_m  = Cal_Signed_Distance(Coor_CD,[X_NODES(iNode),Y_NODES(iNode)]);
					distance_Gauss_m = Cal_Signed_Distance(Coor_CD,Global_coor_Gauss);
					distance_x0  = Cal_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,x0);
					
					
					H_gp_M = sign(distance_Gauss_M);
					H_i_M  = sign(distance_Node_M); 
					H_gp_m = sign(distance_Gauss_m);
					H_i_m  = sign(distance_Node_m); 
					H_x0   = sign(distance_x0);
					if H_gp_M*H_x0 > 0
						if H_gp_m > 0
							Cross_gp=  1;
						else
							Cross_gp= -1;		
						end
					else
						Cross_gp=  0;
					end
					if H_i_M*H_x0 > 0
						if H_i_m > 0
							Cross_i=  1;
						else
							Cross_i= -1;
						end
					else
						Cross_i=  0;
					end
					
					  Cross_gp
					  Cross_i
                      % distance_Gauss =  Cal_Signed_Distance(Coor_AB,Global_coor_Gauss)
                      % Cross_gp       =  sign(distance_Gauss)

                      % distance_Node =  Cal_Signed_Distance(Coor_AB,[X_NODES(iNode),Y_NODES(iNode)])
                      % Cross_i = sign(distance_Node)
					  
	
								  
					      Disp_x = Disp_x + (Cross_gp-Cross_i)*tem_N(iNode)*DISP(Num_Enriched_Node,2);
					      Disp_y = Disp_y + (Cross_gp-Cross_i)*tem_N(iNode)*DISP(Num_Enriched_Node,3);
					end
				end
			end
		end
	end
end

% *************************************
% Loop through each circle inclusion.
% *************************************
%说明:因为材料界面处的增强函数值为0,�?��以下代码不起作用,�?��注释掉了
% if num_Circ_Inclusion ~=0
	% for iIncl=1:num_Circ_Inclusion									
		% if any(Enriched_Node_Type_Incl(NODES_iElem,iIncl)) ==0;
			% Disp_x = Disp_x;
			% Disp_y = Disp_y;
		% elseif any(Enriched_Node_Type_Incl(NODES_iElem,iIncl)) ~=0;
            % c_Incl_x = Circ_Inclu_Coor(iIncl,1);
            % c_Incl_y = Circ_Inclu_Coor(iIncl,2);
            % c_Incl_r = Circ_Inclu_Coor(iIncl,3);
            % Zeta_Node(1:4) = 0;
            % for i_N = 1,4; 
                % c_NODE  = Elem_Node(iElem,i_N); 
                % c_Dis_N = sqrt((X_NODES(i_N)-c_Incl_x)^2+(Y_NODES(i_N)-c_Incl_y)^2);
                % Zeta_Node(i_N) = c_Dis_N-c_Incl_r;
                % if abs(Zeta_Node(i_N))<1.0D-8 
                    % Zeta_Node(i_N) = 0;
                % end
            % end
            % Zeta_Node_abs = abs(Zeta_Node);
            % c_N(1) = tem_N(1);
            % c_N(2) = tem_N(2);
            % c_N(3) = tem_N(3);
            % c_N(4) = tem_N(4);
            % Zm = Zeta_Node(1)*tem_N(1)+ Zeta_Node(2)*tem_N(2)+Zeta_Node(3)*tem_N(3) + Zeta_Node(4)*tem_N(4);
            % Za = Zeta_Node_abs(1)*tem_N(1)+Zeta_Node_abs(2)*tem_N(2)+ Zeta_Node_abs(3)*tem_N(3)+Zeta_Node_abs(4)*tem_N(4);
			% for iNode = 1:4
				% if Enriched_Node_Type_Incl(NODES_iElem(iNode),iIncl) ==1  
                    % if Elem_Type_Incl(iElem,iIncl) == 1			
						% Num_Enriched_Node = POS_Incl(NODES_iElem(iNode),iIncl);
						% c_Dis_N = sqrt((X_NODES(iNode)-c_Incl_x)^2+(Y_NODES(iNode)-c_Incl_y)^2);

						% Incl_gp = Za-abs(Zm);   %实际�?圆上的Incl_gp=0
						% Disp_x = Disp_x + Incl_gp*tem_N(iNode)*DISP(Num_Enriched_Node,2);
						% Disp_y = Disp_y + Incl_gp*tem_N(iNode)*DISP(Num_Enriched_Node,3);
					% end
				% else
				    % Disp_x = Disp_x;
					% Disp_y = Disp_y;
				% end
			% end
		% end
	% end
% end