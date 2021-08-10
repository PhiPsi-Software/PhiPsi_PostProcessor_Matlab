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

function Plot_Node_Disp(DISP,isub,Crack_X,Crack_Y,POS,Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
					      Coors_Vertex,Coors_Junction,Coors_Tip,Crack_Tip_Type,Shaped_Crack_Points)
% This function plots the displacement contours of all nodes.

global Node_Coor Elem_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor
global Key_PLOT Num_Node Key_Flipped_Gray
global Size_Font num_Crack Outline Inline
global Num_Contourlevel Full_Pathname aveg_area_ele
global Color_Crack Width_Crack Num_Accuracy_Contour Key_Contour_Metd
global Color_Contourlevel Color_Mesh
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global num_Hole Hole_Coor
global num_Circ_Inclusion Circ_Inclu_Coor
global Color_Inclusion
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global num_Cross Cross_Coor Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Arc_Crack_Coor Yes_Arc_Crack
global Killed_Elements num_Killed_Ele Ellipse_Hole_Coor num_Ellipse_Hole
global Num_Colorbar_Level Colorbar_Font
global Title_Font Key_Figure_Control_Widget

disp('    > Plotting displacement contours....') 

scale = Key_PLOT(4,6);

% Get the new coordinates of all nodes
New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);

% Get the maximum and minimum value of the new coordinates of all nodes
Min_X_Coor_New = min(min(New_Node_Coor(:,1)));
Max_X_Coor_New = max(max(New_Node_Coor(:,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(:,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(:,2)));

% Get resample coors.
delta = sqrt(aveg_area_ele)/Num_Accuracy_Contour;
if Key_PLOT(4,8)==0
	gx    = Min_X_Coor:delta:Max_X_Coor; 
	gy    = Min_Y_Coor:delta:Max_Y_Coor;
elseif Key_PLOT(4,8)==1
	gx    = Min_X_Coor_New:delta:Max_X_Coor_New; 
	gy    = Min_Y_Coor_New:delta:Max_Y_Coor_New;
end

% 绘图设置.
Start_Dis_Type =1;      % Plot Sxx, Sxy, Syy and Mises.
End_Dis_Type   =2; 
if Key_PLOT(4,2) ==1       % Only plot位移分量平方和开根号.
    Start_Dis_Type =3;
	End_Dis_Type   =3;
end

% Plot the contours.
for i = Start_Dis_Type:End_Dis_Type
    % New figure.
    Tools_New_Figure
	hold on;
	if i == 1
	    if Key_PLOT(4,8)==0
		    disp('      Resample Ux....') 
			if 	Key_Contour_Metd==1
			    [X,Y,Ux] = griddata(Node_Coor(:,1),Node_Coor(:,2),DISP(1:Num_Node,2),...
	                    unique(Node_Coor(:,1)),unique(Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [Ux,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),DISP(1:Num_Node,2),gx,gy);
            end				
		elseif Key_PLOT(4,8)==1
			disp('      Resample Ux....') 
			if 	Key_Contour_Metd==1
			    [X,Y,Ux] = griddata(New_Node_Coor(:,1),New_Node_Coor(:,2),DISP(1:Num_Node,2),...
	                    unique(New_Node_Coor(:,1)),unique(New_Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [Ux,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),DISP(1:Num_Node,2),gx,gy);
            end				
		end
		disp('      Contouring Ux....') 
		contourf(X,Y,Ux,Num_Contourlevel,'LineStyle','none')
		clear Ux
	elseif i == 2
	    if Key_PLOT(4,8)==0
        	disp('      Resample Uy....') 
			if 	Key_Contour_Metd==1
			[X,Y,Uy] = griddata(Node_Coor(:,1),Node_Coor(:,2),DISP(1:Num_Node,3),...
	                    unique(Node_Coor(:,1)),unique(Node_Coor(:,2))');	
			elseif 	Key_Contour_Metd==2
			    [Uy,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),DISP(1:Num_Node,3),gx,gy);	
            end			
		elseif Key_PLOT(4,8)==1
        	disp('      Resample Uy....') 
			if 	Key_Contour_Metd==1
				[X,Y,Uy] = griddata(New_Node_Coor(:,1),New_Node_Coor(:,2),DISP(1:Num_Node,3),...
	                    unique(New_Node_Coor(:,1)),unique(New_Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [Uy,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),DISP(1:Num_Node,3),gx,gy);
			end
        end		
		disp('      Contouring Uy....') 
		contourf(X,Y,Uy,Num_Contourlevel,'LineStyle','none')
		clear Uy
	%位移分量平方和开根号
	elseif i == 3
        c_SQRT_Ux_Uy =sqrt(DISP(1:Num_Node,2).*DISP(1:Num_Node,2) + DISP(1:Num_Node,3).*DISP(1:Num_Node,3));
	    if Key_PLOT(4,8)==0
        	disp('      Resample sqrt(dis_x^2 + dis_y^2)....') 
			if 	Key_Contour_Metd==1
			[X,Y,Uy] = griddata(Node_Coor(:,1),Node_Coor(:,2),c_SQRT_Ux_Uy,...
	                    unique(Node_Coor(:,1)),unique(Node_Coor(:,2))');	
			elseif 	Key_Contour_Metd==2
			    [Uy,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),c_SQRT_Ux_Uy,gx,gy);	
            end			
		elseif Key_PLOT(4,8)==1
        	disp('      Resample sqrt(dis_x^2 + dis_y^2)....') 
			if 	Key_Contour_Metd==1
				[X,Y,Uy] = griddata(New_Node_Coor(:,1),New_Node_Coor(:,2),c_SQRT_Ux_Uy,...
	                    unique(New_Node_Coor(:,1)),unique(New_Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [Uy,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),c_SQRT_Ux_Uy,gx,gy);
			end
        end		
		disp('      Contouring sqrt(dis_x^2 + dis_y^2)....') 
		contourf(X,Y,Uy,Num_Contourlevel,'LineStyle','none')
		clear Uy
	end
	
	% Set colormap.
	if Key_Flipped_Gray==0
		colormap(Color_Contourlevel)
	elseif Key_Flipped_Gray==1
		colormap(flipud(gray))
	end
	
	% Plot the outline of the mesh
	% if Key_PLOT(4,8)==0
		% line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
			 % [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
	% elseif Key_PLOT(4,8)==1
		% line([New_Node_Coor(Outline(:,1),1) New_Node_Coor(Outline(:,2),1)]', ...
			 % [New_Node_Coor(Outline(:,1),2) New_Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
	% end
	
	% The range of the plot.
	if Key_PLOT(4,8)==0
	    min_x = Min_X_Coor; max_x = Max_X_Coor; min_y = Min_Y_Coor; max_y = Max_Y_Coor;
	elseif Key_PLOT(4,8)==1
	    min_x = Min_X_Coor_New; max_x = Max_X_Coor_New;
		min_y = Min_Y_Coor_New; max_y = Max_Y_Coor_New;
	end
	axis([min_x max_x min_y max_y]);
	
	range_x = max_x -min_x; 
    range_y = max_y -min_y; 
	% Fill the outside of the domain by white color.
	if Key_PLOT(4,8)==0
	    Tools_Fillout(Node_Coor(Outline(:,1),1),Node_Coor(Outline(:,1),2),[min_x-0.01*range_x max_x+0.01*range_x min_y-0.01*range_y max_y+0.01*range_y],'w'); 
	elseif Key_PLOT(4,8)==1
	    Tools_Fillout(New_Node_Coor(Outline(:,1),1),New_Node_Coor(Outline(:,1),2),[min_x-0.01*range_x max_x+0.01*range_x min_y-0.01*range_y max_y+0.01*range_y],'w'); 
	end
	
	% Fill the inside of the domain by white color, holes, for example.
	if isempty(Inline)==0
		if Key_PLOT(4,8)==0
			fill(Node_Coor(Inline(:,1),1),Node_Coor(Inline(:,1),2),'w');
		elseif Key_PLOT(4,8)==1 
			fill(New_Node_Coor(Inline(:,1),1),New_Node_Coor(Inline(:,1),2),'w')	
		end
	end	
	
	% Plot the mesh
    if Key_PLOT(4,9) ==1
		if Key_PLOT(4,8)==0
			for iElem =1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)]; % Nodes for current element
				xi = Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
				yi = Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
			end
		elseif Key_PLOT(4,8)==1
			for iElem =1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)]; % Nodes for current element
				xi = New_Node_Coor(NN',1);                                       % Deformed x-coordinates of nodes
				yi = New_Node_Coor(NN',2);                                       % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
			end
		end
    end	
	% Plot cracks line if necessary
	if Key_PLOT(4,5) == 1
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
								c_alpha=c_Radian_Start:pi/50:c_Radian_End;
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
								c_alpha=c_Radian_Start:pi/50:2*pi;
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
								%第2部分:0到Radian_End
								c_alpha_2=0.0:pi/50:c_Radian_End;
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
								c_alpha=c_Radian_End:pi/50:2*pi;
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
		for iii = 1:num_Hole
			Coor_x  = Hole_Coor(iii,1);
			Coor_y  = Hole_Coor(iii,2);
			c_R  = Hole_Coor(iii,3);
			num_fineness = 100;
			for j = 1:num_fineness+1
				alpha = 2*pi/num_fineness*(j-1);
				x(j) = Coor_x + c_R*cos(alpha);
				y(j) = Coor_y + c_R*sin(alpha);
				[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(j),y(j));
                %Sometimes the element does not exist, for example, only part of the hole locates inside the model	
			    if c_Elem_Num~=0
				    [Kesi,Yita] = Cal_KesiYita_by_Coors(x(j),y(j));
			    	[dis_x(j),dis_y(j)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y);
                end																  
			end
			x_new = x + dis_x*Key_PLOT(4,6);
			y_new = y + dis_y*Key_PLOT(4,6);
			% plot(x_new,y_new,'-')
			patch(x_new,y_new,'white','edgecolor','black','LineWidth',0.1)	
		end	
	end
	
	% Plot ellipse holes, deformation has not been considered (2020-08-10).
	if num_Ellipse_Hole ~=0
		disp(['      ----- Plotting ellipse hole...'])
		for i_H = 1:num_Ellipse_Hole
			Coor_x  = Ellipse_Hole_Coor(i_H,1);
			Coor_y  = Ellipse_Hole_Coor(i_H,2);
			a_Hole  = Ellipse_Hole_Coor(i_H,3);
			b_Hole  = Ellipse_Hole_Coor(i_H,4);
			c_Hole  = sqrt(a_Hole^2-b_Hole^2);
			theta_Hole  = Ellipse_Hole_Coor(i_H,5);
			% Plot the ellipse.
			hold on;
			syms x y;
			theta_Hole = theta_Hole*pi/180.0;
			f=(a_Hole^2-c_Hole^2*cos(theta_Hole)^2)*(x-Coor_x).^2+(a_Hole^2-c_Hole^2*sin(theta_Hole)^2)*(y-Coor_y).^2-c_Hole^2*sin(2*theta_Hole)*(x-Coor_x).*(y-Coor_y)-a_Hole^2*b_Hole^2;
			h=ezplot(f,[-1.5*a_Hole+Coor_x,1.5*a_Hole+Coor_x,-1.5*a_Hole+Coor_y,1.5*a_Hole+Coor_y]);
			h.Color = 'black';
            % Fill the ellipse.
			data=get(h,'contourMatrix');
			x=data(1,2:end);
			y=data(2,2:end);
			fill(x,y,'white')	
		end	
	end
	
	% 绘制圆形夹杂.
	if num_Circ_Inclusion ~=0
	    disp(['      ----- Plotting circle inclusion...'])
		for iii = 1:num_Circ_Inclusion
			Coor_x  = Circ_Inclu_Coor(iii,1);
			Coor_y  = Circ_Inclu_Coor(iii,2);
			c_R  = Circ_Inclu_Coor(iii,3);
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
			x_new = x + dis_x*Key_PLOT(4,6);
			y_new = y + dis_y*Key_PLOT(4,6);
			% plot(x_new,y_new,'-')
			plot(x_new,y_new,'-','color','black')
		end	
	end
	% 绘制多边形夹杂,2016-10-04
	if num_Poly_Inclusion ~=0
		disp(['      ----- Plotting poly inclusion...'])
		for iii = 1:num_Poly_Inclusion
			nEdge = size(Poly_Incl_Coor_x{iii},2);
			Num_Diversion = 5;    %多边形的每条边拆分成5份
			k = 0;
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
				k =k+1;
				t_x(k) = a_x;
				t_y(k) = a_y;
				[Kesi,Yita] = Cal_KesiYita_by_Coors(t_x(k),t_y(k));
				[c_Elem_Num] = Cal_Ele_Num_by_Coors(t_x(k),t_y(k));
				[t_dis_x(k),t_dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				%计算等分点的位移
				for jjj =  1:Num_Diversion-1
					k=k+1;
					t_x(k) = (jjj*b_x+(Num_Diversion-jjj)*a_x)/Num_Diversion;
					t_y(k) = (jjj*b_y+(Num_Diversion-jjj)*a_y)/Num_Diversion;
					[Kesi,Yita] = Cal_KesiYita_by_Coors(t_x(k),t_y(k));
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(t_x(k),t_y(k));
					[t_dis_x(k),t_dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				end												      
			end
			x_new = t_x + t_dis_x*Key_PLOT(4,6);
			y_new = t_y + t_dis_y*Key_PLOT(4,6);
			................
			% option1:填充
			%................
			patch(x_new,y_new,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
			................
			% option2:绘制边线
			%................
			% x_new(k+1) = x_new(1);
			% y_new(k+1) = y_new(1);
			% plot(x_new,y_new,'-','color','black')
		end	
	end	
	% Plot shaped cracks if necessary
	if Key_PLOT(4,5) == 2 & num_Crack(isub)~=0
	    disp(['      > Plotting shaped cracks......'])
	    Plot_Shaped_Cracks(Shaped_Crack_Points)
	end
	%绘制支撑剂的圆球
	if Key_PLOT(4,10)==1 
		%如果支撑剂直径和坐标文件存在：
		if exist([Full_Pathname,'.epcr_',num2str(isub)], 'file') ==2 
			disp(['      ----- Plotting proppant......'])  
			%读取文件		
			c_D_Coors=load([Full_Pathname,'.epcr_',num2str(isub)]);
			num_Proppant = size(c_D_Coors,1);
			for i=1:num_Proppant
				%绘制实心圆
				alpha=0:pi/20:2*pi;
				c_Elem_Num = c_D_Coors(i,1);
				R = c_D_Coors(i,2)/2*Key_PLOT(4,6);
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
																  Coors_Vertex,Coors_Junction,Coors_Tip);
				[Kesi_Low,Yita_Low] = Cal_KesiYita_by_Coors(old_Coor_x_Low,old_Coor_y_Low);
				[dis_x_Low,dis_y_Low] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi_Low,Yita_Low...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip);
				%真实的支撑剂所在点的位移											  
				c_dis_x = (dis_x_Up + dis_x_Low)/2;			
				c_dis_y = (dis_y_Up + dis_y_Low)/2;	
				c_Coor_x = old_Coor_x  + c_dis_x*Key_PLOT(4,6);
				c_Coor_y = old_Coor_y  + c_dis_y*Key_PLOT(4,6);
				
				c_x = c_Coor_x + R*cos(alpha);
				c_y = c_Coor_y + R*sin(alpha);
				plot(c_x,c_y,'-')
				axis equal
				fill(c_x,c_y,[128/255,138/255,135/255])
			end
		end
	end
	%绘制破裂区
	if Key_PLOT(4,15)==1 
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
	%杀死的单元用白色填充
	if num_Killed_Ele>0
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
			if ismember(iElem,Killed_Elements)
				xi_killed(:) = New_Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
				yi_killed(:) = New_Node_Coor(NN',2);  
				if Key_PLOT(4,9) == 1
				    patch(xi_killed,yi_killed,'white') 
				elseif Key_PLOT(4,9) == 0
				    patch(xi_killed,yi_killed,'white','EdgeColor','none') 
				end
			end
		end 
    end
	
    % Set Title.
    if i == 1
		title('Displacement x','FontName',Title_Font,'FontSize',Size_Font)
	elseif i == 2
		title('Displacement y','FontName',Title_Font,'FontSize',Size_Font)
	%位移分量平方和开根号
	elseif i == 3
		title('Squared displacement: sqrt(disx^2 + dispy^2)','FontName',Title_Font,'FontSize',Size_Font)
	end
	
	axis equal; 
	%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
	% colorbar('FontName',Title_Font,'FontSize',Size_Font);
	t1=caxis;
	t1=linspace(t1(1),t1(2),Num_Colorbar_Level);
	my_handle=colorbar('FontName',Colorbar_Font,'FontSize',Size_Font,'ytick',t1);	
    set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')


	% Active Figure control widget (2021-08-01)
	% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
	% Press q to exit.
	% Press r (or double-click) to reset to the initial.
	if Key_Figure_Control_Widget==1
		fcw(gca);
	end
	
	% Save pictures.
    if i == 1
	    Save_Picture(c_figure,Full_Pathname,'Dpxn')
	elseif i ==2
	    Save_Picture(c_figure,Full_Pathname,'Dpyn')
	end	
end


