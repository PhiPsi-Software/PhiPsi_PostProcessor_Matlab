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

function Plot_Node_Field_Value(DISP,Field_Value,isub,Crack_X,Crack_Y,POS,Enriched_Node_Type,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,...
					           Node_Cross_elem,Coors_Vertex,Coors_Junction,Coors_Tip,Crack_Tip_Type,Shaped_Crack_Points)
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
global Field_Flux_x Field_Flux_y
global max_area_ele
global num_Cross Cross_Coor Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Num_Colorbar_Level Colorbar_Font
global Title_Font Key_Figure_Control_Widget

disp('    > Plotting field problem....') 

scale = Key_PLOT(5,6);
% Get the new coordinates of all nodes
if isempty(DISP) ==0
    New_Node_Coor(:,1) = Node_Coor(:,1) + scale*DISP(1:Num_Node,2);
    New_Node_Coor(:,2) = Node_Coor(:,2) + scale*DISP(1:Num_Node,3);
else
    New_Node_Coor(:,1) = Node_Coor(:,1);
    New_Node_Coor(:,2) = Node_Coor(:,2);
end
% Get the maximum and minimum value of the new coordinates of all nodes
Min_X_Coor_New = min(min(New_Node_Coor(:,1)));
Max_X_Coor_New = max(max(New_Node_Coor(:,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(:,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(:,2)));

% Get resample coors.
delta = sqrt(aveg_area_ele)/Num_Accuracy_Contour;
if Key_PLOT(5,8)==0
	gx    = Min_X_Coor:delta:Max_X_Coor; 
	gy    = Min_Y_Coor:delta:Max_Y_Coor;
elseif Key_PLOT(5,8)==1
	gx    = Min_X_Coor_New:delta:Max_X_Coor_New; 
	gy    = Min_Y_Coor_New:delta:Max_Y_Coor_New;
end

% Plot the contours
for i = 1:1
    %...............
    % New figure.
	%...............
    Tools_New_Figure
	hold on;
	%...............
	%绘制标量场云图
	%...............
	hold on;
	if Key_PLOT(5,1)==1 | Key_PLOT(5,1)==3
	    if Key_PLOT(5,8)==0
		    disp('      Resample field value....') 
			if 	Key_Contour_Metd==1
			    [X,Y,field] = griddata(Node_Coor(:,1),Node_Coor(:,2),Field_Value(1:Num_Node),...
	                    unique(Node_Coor(:,1)),unique(Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [field,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Field_Value(1:Num_Node),gx,gy);
            end				
		elseif Key_PLOT(5,8)==1
			disp('      Resample field value....') 
			if 	Key_Contour_Metd==1
			    [X,Y,field] = griddata(New_Node_Coor(:,1),New_Node_Coor(:,2),Field_Value(1:Num_Node),...
	                    unique(New_Node_Coor(:,1)),unique(New_Node_Coor(:,2))');
			elseif 	Key_Contour_Metd==2
			    [field,X,Y] = Tools_gridfit(New_Node_Coor(:,1),New_Node_Coor(:,2),Field_Value(1:Num_Node),gx,gy);
            end				
		end
		disp('      Contouring field value....') 
		contourf(X,Y,field,Num_Contourlevel,'LineStyle','none')
		title('Field value','FontName',Title_Font,'FontSize',Size_Font)
		clear field
		
		%Set colormap.
		%colormap(gray)
		%colormap(hot)
		%colormap(cool)
		colormap(jet)
		%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
		% colorbar('FontName',Title_Font,'FontSize',Size_Font);
	    t1=caxis;
	    t1=linspace(t1(1),t1(2),Num_Colorbar_Level);
	    my_handle=colorbar('FontName',Colorbar_Font,'FontSize',Size_Font,'ytick',t1);		
		set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
	end
	%...................
	% Plot the mesh
	%...................
    if Key_PLOT(5,9) ==1
		if Key_PLOT(5,8)==0
			for iElem =1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)]; % Nodes for current element
				xi = Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
				yi = Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
			end
		elseif Key_PLOT(5,8)==1
			for iElem =1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
					  Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)];     % Nodes for current element
				xi = New_Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
				yi = New_Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
			end
		end
    end	
	%...................
	%绘制流量场矢量图
	%...................
	if Key_PLOT(5,1)==2 | Key_PLOT(5,1)==3
		Max_Flux_x = max(abs(Field_Flux_x));
		Max_Flux_y = max(abs(Field_Flux_y));
		Max_Flux   = max(Max_Flux_x,Max_Flux_y);
		
		W = Max_X_Coor - Min_X_Coor;
		H = Max_Y_Coor - Min_Y_Coor;
		
		% length of force arrow
		length_arrow = sqrt(max_area_ele)*0.8;
		%length_arrow = max(W,H)/15.0;          
		
		% Loop through each node.
		for i = 1:Num_Node
			c_flux_x   = Field_Flux_x(i);
			c_flux_y   = Field_Flux_y(i);

			delta_L_x = c_flux_x*length_arrow/Max_Flux;
			delta_L_y = c_flux_y*length_arrow/Max_Flux;
			
			StartPoint = [Node_Coor(i,1)-delta_L_x   Node_Coor(i,2)-delta_L_y     0];
			EndPoint   = [Node_Coor(i,1)             Node_Coor(i,2)               0];

			line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','black','LineWidth',1.5)
			
			% The length of the head of the arrow.
			length_arrow_head = length_arrow/3;
			
			% Plot the head of the arrow.
			theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
			theta_1 = pi/2 - theta - pi/3;
			delta_x = -length_arrow_head*cos(theta_1);
			delta_y =  length_arrow_head*sin(theta_1);
			line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','black','LineWidth',1.5);
			theta_2 = 3*pi/2 - theta + pi/3;
			delta_x = -length_arrow_head*cos(theta_2);
			delta_y =  length_arrow_head*sin(theta_2);
			line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','black','LineWidth',1.5);
		end
		set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
	end		
	
	% Plot the outline of the mesh
	if Key_PLOT(5,8)==0
		line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
			 [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
	elseif Key_PLOT(5,8)==1
		line([New_Node_Coor(Outline(:,1),1) New_Node_Coor(Outline(:,2),1)]', ...
			 [New_Node_Coor(Outline(:,1),2) New_Node_Coor(Outline(:,2),2)]','LineWidth',1,'Color','black')
	end
	
	% The range of the plot.
	if Key_PLOT(5,8)==0
	    min_x = Min_X_Coor; max_x = Max_X_Coor; min_y = Min_Y_Coor; max_y = Max_Y_Coor;
	elseif Key_PLOT(5,8)==1
	    min_x = Min_X_Coor_New; max_x = Max_X_Coor_New;
		min_y = Min_Y_Coor_New; max_y = Max_Y_Coor_New;
	end
	axis([min_x max_x min_y max_y]);	
	
	% Fill the outside of the domain by white color.
	% if Key_PLOT(5,8)==0
	    % Tools_Fillout(Node_Coor(Outline(:,1),1),Node_Coor(Outline(:,1),2),[min_x max_x min_y max_y],'w'); 
	% elseif Key_PLOT(5,8)==1
	    % Tools_Fillout(New_Node_Coor(Outline(:,1),1),New_Node_Coor(Outline(:,1),2),[min_x max_x min_y max_y],'w'); 
	% end
	
	% Fill the inside of the domain by white color, holes, for example.
	if isempty(Inline)==0
		if Key_PLOT(4,8)==0
			fill(Node_Coor(Inline(:,1),1),Node_Coor(Inline(:,1),2),'w');
		elseif Key_PLOT(4,8)==1 
			fill(New_Node_Coor(Inline(:,1),1),New_Node_Coor(Inline(:,1),2),'w')	
		end
	end	

	% Plot cracks line if necessary
	if Key_PLOT(5,5) == 1
		if num_Crack(isub)~=0
			for i = 1:num_Crack(isub)
				nPt = size(Crack_X{i},2);
				for iPt = 2:nPt
					x = [Crack_X{i}(iPt-1) Crack_X{i}(iPt)];
					y = [Crack_Y{i}(iPt-1) Crack_Y{i}(iPt)];
					if isempty(DISP) ==0
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
					else
						last_x = [ x(1)  x(2)];
						last_y = [ y(1)  y(2)];					
					end
					
					plot(last_x,last_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   
				end
			end	
		end
	end
	% Plot holes.
	disp(['      ----- Plotting hole...'])
	if num_Hole ~=0
		for iii = 1:num_Hole
			Coor_x  = Hole_Coor(iii,1);
			Coor_y  = Hole_Coor(iii,2);
			c_R  = Hole_Coor(iii,3);
			num_fineness = 100;
			if isempty(DISP) ==0
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
				x_new = x + dis_x*Key_PLOT(5,6);
				y_new = y + dis_y*Key_PLOT(5,6);
			else
				for j = 1:num_fineness+1
					alpha = 2*pi/num_fineness*(j-1);
					x(j) = Coor_x + c_R*cos(alpha);
					y(j) = Coor_y + c_R*sin(alpha);
				end
				x_new = x;  
				y_new = y;
			end
			% plot(x_new,y_new,'-')
			patch(x_new,y_new,'white','edgecolor','black','LineWidth',0.1)	
		end	
	end
	% 绘制圆形夹杂.
	disp(['      ----- Plotting circle inclusion...'])
	if num_Circ_Inclusion ~=0
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
			x_new = x + dis_x*Key_PLOT(5,6);
			y_new = y + dis_y*Key_PLOT(5,6);
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
				x(k) = a_x;
				y(k) = a_y;
				[Kesi,Yita] = Cal_KesiYita_by_Coors(x(k),y(k));
				[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(k),y(k));
				[dis_x(k),dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				%计算等分点的位移
				for j =  1:Num_Diversion-1
					k=k+1;
					x(k) = (j*b_x+(Num_Diversion-j)*a_x)/Num_Diversion;
					y(k) = (j*b_y+(Num_Diversion-j)*a_y)/Num_Diversion;
					[Kesi,Yita] = Cal_KesiYita_by_Coors(x(k),y(k));
					[c_Elem_Num] = Cal_Ele_Num_by_Coors(x(k),y(k));
					[dis_x(k),dis_y(k)] = Cal_Anypoint_Disp(c_Elem_Num,Enriched_Node_Type,POS,isub,DISP,Kesi,Yita...
																 ,Elem_Type,Coors_Element_Crack,Node_Jun_elem,Node_Jun_Hole,Node_Cross_elem,...
																  Coors_Vertex,Coors_Junction,Coors_Tip,Crack_X,Crack_Y); 
				end												      
			end
			x_new = x + dis_x*Key_PLOT(5,6);
			y_new = y + dis_y*Key_PLOT(5,6);
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
	if Key_PLOT(5,5) == 2 & num_Crack(isub)~=0
	    disp(['      > Plotting shaped cracks......'])
	    Plot_Shaped_Cracks(Shaped_Crack_Points)
	end
	%绘制支撑剂的圆球
	if Key_PLOT(5,10)==1 
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
				c_Coor_x = old_Coor_x  + c_dis_x*Key_PLOT(5,6);
				c_Coor_y = old_Coor_y  + c_dis_y*Key_PLOT(5,6);
				
				c_x = c_Coor_x + R*cos(alpha);
				c_y = c_Coor_y + R*sin(alpha);
				plot(c_x,c_y,'-')
				axis equal
				fill(c_x,c_y,[128/255,138/255,135/255])
			end
		end
	end
	%绘制破裂区
	if Key_PLOT(5,15)==1 
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
	axis equal; 

	% Active Figure control widget (2021-08-01)
	% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
	% Press q to exit.
	% Press r (or double-click) to reset to the initial.
	if Key_Figure_Control_Widget==1
		fcw(gca);
	end
	
	% Save pictures.
    if i == 1
	    Save_Picture(c_figure,Full_Pathname,'field_value')
	elseif i ==2
	end	
end


