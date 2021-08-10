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

function Plot_Mesh(isub,Crack_X,Crack_Y,Post_Enriched_Nodes,POS)
% This function plots the initial geometry.

global Node_Coor Elem_Node Key_POST_HF
global Num_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor
global Key_PLOT aveg_area_ele Outline
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global num_Crack num_of_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Elem_Material Num_Step_to_Plot
global Na_Crack_X Na_Crack_Y num_Na_Crack
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl
global Color_Inclusion
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl Elem_Type_Incl POS_Incl
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Yes_Field_Problem
global Field_Boundary_Value
global Field_Boundary_Qn
global num_Cross Enriched_Node_Type_Cross POS_Cross Elem_Type_Cross
global Arc_Crack_Coor Yes_Arc_Crack
global Key_Data_Format
global num_HC HC_Coor Enriched_Node_Type_HC Node_HC_elem POS_HC Elem_Type_HC
global Killed_Elements num_Killed_Ele Ellipse_Hole_Coor num_Ellipse_Hole
global Title_Font Key_Figure_Control_Widget

% Preparing.
xi = zeros(4,Num_Elem); yi = xi;
disp(['      ----- Plotting undeformed mesh......'])
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
for iElem = 1:Num_Elem
    NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
	      Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
	if Elem_Material(iElem)==1
		xi_1(:,iElem) = Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
		yi_1(:,iElem) = Node_Coor(NN',2);                                     % Initial y-coordinates of nodes		
	elseif Elem_Material(iElem)==2
		xi_2(:,iElem) = Node_Coor(NN',1);                                     
		yi_2(:,iElem) = Node_Coor(NN',2);   
	elseif Elem_Material(iElem)==3
		xi_3(:,iElem) = Node_Coor(NN',1);                                     
		yi_3(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==4
		xi_4(:,iElem) = Node_Coor(NN',1);                                     
		yi_4(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==5
		xi_5(:,iElem) = Node_Coor(NN',1);                                     
		yi_5(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==6
		xi_6(:,iElem) = Node_Coor(NN',1);                                     
		yi_6(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==7
		xi_7(:,iElem) = Node_Coor(NN',1);                                     
		yi_7(:,iElem) = Node_Coor(NN',2); 	
	elseif Elem_Material(iElem)==8
		xi_8(:,iElem) = Node_Coor(NN',1);                                     
		yi_8(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==9
		xi_9(:,iElem) = Node_Coor(NN',1);                                     
		yi_9(:,iElem) = Node_Coor(NN',2); 
	elseif Elem_Material(iElem)==10
		xi_10(:,iElem) = Node_Coor(NN',1);                                     
		yi_10(:,iElem) = Node_Coor(NN',2); 		
	end
end

% New figure.
Tools_New_Figure
hold on;
title('Finite Element Mesh','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
delta = sqrt(aveg_area_ele);
axis([Min_X_Coor-delta Max_X_Coor+delta Min_Y_Coor-delta Max_Y_Coor+delta]);

% Plot elements.  
disp(['      ----- Plotting elements...'])
if Key_PLOT(1,9) ==1
	patch(xi_1,yi_1,Color_Backgro_Mesh_1) 
	patch(xi_2,yi_2,Color_Backgro_Mesh_2)    
	patch(xi_3,yi_3,Color_Backgro_Mesh_3)   
	patch(xi_4,yi_4,Color_Backgro_Mesh_4)  
	patch(xi_5,yi_5,Color_Backgro_Mesh_5)  
	patch(xi_6,yi_6,Color_Backgro_Mesh_6)  
	patch(xi_7,yi_7,Color_Backgro_Mesh_7)     
	patch(xi_8,yi_8,Color_Backgro_Mesh_8) 
	patch(xi_9,yi_9,Color_Backgro_Mesh_9) 
	patch(xi_10,yi_10,Color_Backgro_Mesh_10)
	%杀死的单元用白色
	if num_Killed_Ele>0
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                                 % Nodes for current element
			if ismember(iElem,Killed_Elements)
				xi_killed(:) = Node_Coor(NN',1);                                     % Initial x-coordinates of nodes
				yi_killed(:) = Node_Coor(NN',2);  
				patch(xi_killed,yi_killed,'white') 
			end
		end 
	end
	
elseif Key_PLOT(1,9) ==0 %绘制模型的外边框
	line([Node_Coor(Outline(:,1),1) Node_Coor(Outline(:,2),1)]', ...
		 [Node_Coor(Outline(:,1),2) Node_Coor(Outline(:,2),2)]','LineWidth',1.0,'Color','black')
end

% 单元接触状态
if Key_PLOT(1,13)==1 
    Yes_Contact1=0;  %裂纹面闭合
	Yes_Contact2=0;  %裂纹面由支撑剂支撑
	%如果接触状态文件存在则：
	if exist([Full_Pathname,'.elcs_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting contact status of elements...'])
		ELCS= load([Full_Pathname,'.elcs_',num2str(Num_Step_to_Plot)]);
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
			if ELCS(iElem) == 1    %裂纹面闭合
				xi_Elcs1(:,iElem) = Node_Coor(NN',1);                                     
				yi_Elcs1(:,iElem) = Node_Coor(NN',2); 	
				Yes_Contact1=1;
			end
			if ELCS(iElem) == 2    %裂纹面由支撑剂支撑
				xi_Elcs2(:,iElem) = Node_Coor(NN',1);                                     
				yi_Elcs2(:,iElem) = Node_Coor(NN',2); 	
				Yes_Contact2=1;
			end
		end
		if Yes_Contact1==1
			patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])               %用粉红色标记该单元
		end
		if Yes_Contact2==1
			patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])         %用紫色标记该单元
		end
	end
	%如果粘聚裂缝粘聚单元状态文件存在则(2017-04-24)：
	if exist([Full_Pathname,'.elco_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting cohesive crack status of elements...'])
		ELCS= load([Full_Pathname,'.elco_',num2str(Num_Step_to_Plot)]);
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
			if ELCS(iElem) == 1   %粘聚单元
				xi_Elcs1(:,iElem) = Node_Coor(NN',1);                                     
				yi_Elcs1(:,iElem) = Node_Coor(NN',2); 	
				Yes_Contact1=1;
			end
			% if ELCS(iElem) == 2    %裂纹面由支撑剂支撑
				% xi_Elcs2(:,iElem) = Node_Coor(NN',1);                                     
				% yi_Elcs2(:,iElem) = Node_Coor(NN',2); 	
				% Yes_Contact2=1;
			% end
		end
		if Yes_Contact1==1
			patch(xi_Elcs1,yi_Elcs1,[1,192/255,203/255])               %用粉红色标记该单元
		end
		if Yes_Contact2==1
			patch(xi_Elcs2,yi_Elcs2,[160/255,102/255,211/255])         %用紫色标记该单元
		end
	end
end
            
% Plot the node numbers.
if Key_PLOT(1,1) ==1 && Key_PLOT(1,2) == 1
    disp(['      ----- Plotting node number...'])
    for iNode = 1:Num_Node
        text(Node_Coor(iNode,1)+0.05*delta,Node_Coor(iNode,2),1,num2str(iNode),...
		              'FontName',Title_Font,'FontSize',Node_Fontsize,'color',Node_Fontcolor)
    end
end

% Plot the element numbers.
if Key_PLOT(1,1) ==1 && Key_PLOT(1,3) == 1
    disp(['      ----- Plotting element number...'])
    for iElem = 1:Num_Elem
        NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
	          Elem_Node(iElem,3) Elem_Node(iElem,4)];
        XN = Node_Coor(NN,1);
        YN = Node_Coor(NN,2);
        text(mean(XN),mean(YN),1,num2str(iElem),'FontName',Title_Font,'FontSize',Elem_Fontsize,'color',Elem_Fontcolor)
    end
end

% 绘制天然裂缝.
% Na_Crack_X
if Key_PLOT(1,12) == 1
    disp(['      ----- Plotting natural crack line...'])
	if isempty(Na_Crack_X)==0
		for i = 1:num_Na_Crack
			nPt = size(Na_Crack_X{i},2);
			for iPt = 2:nPt
				x = [Na_Crack_X{i}(iPt-1) Na_Crack_X{i}(iPt)];
				y = [Na_Crack_Y{i}(iPt-1) Na_Crack_Y{i}(iPt)]; 
				% plot(x,y,'--','Color',Color_Crack,'LineWidth',Width_Crack)  
				plot(x,y,'w','Color','m','LineWidth',Width_Crack) 
			end
		end	
	end
end

% Plot cracks if necessary.
if Key_PLOT(1,5) == 1 || Key_PLOT(1,5) == 2
    disp(['      ----- Plotting crack line...'])
	if isempty(Crack_X)==0
		for i = 1:num_Crack(isub)
			nPt = size(Crack_X{i},2);
			for i_Seg = 1:nPt-1
				x = [Crack_X{i}(i_Seg) Crack_X{i}(i_Seg+1)];
				y = [Crack_Y{i}(i_Seg) Crack_Y{i}(i_Seg+1)];
				if Key_PLOT(1,5) == 2
					% Plot the kink of the crack lines.
					plot(x,y,'ko','MarkerSize',10,'Color','red','linewidth',1.5)%################
				end
				%----------------------------
				%如果是弧形裂缝,2017-07-17
				%----------------------------
				%Arc_Crack_Coor:x,y,r,Radian_Start,Radian_End,Radian,Point_Start_x,Point_Start_y,Point_End_x,Point_End_y
			    if abs(sum(Arc_Crack_Coor(i,i_Seg,1:11)))>=1.0e-10 
					c_R           =Arc_Crack_Coor(i,i_Seg,4);
					c_Direcction  =Arc_Crack_Coor(i,i_Seg,3);
					c_Radian_Start=Arc_Crack_Coor(i,i_Seg,5)*pi/180;
					c_Radian_End  =Arc_Crack_Coor(i,i_Seg,6)*pi/180;
					%c_Radian_End  =2.0*pi
					%**********************
					%   如果是逆时针圆弧
					%**********************
					if c_Direcction >0.5
						%若结束角度大于起始角度,则直接绘制圆弧
						if c_Radian_End>=c_Radian_Start 
							c_alpha=c_Radian_Start:pi/100:c_Radian_End;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						%若结束角度小于起始角度,则分成两部分绘制圆弧
						else
							%第1部分:Radian_Start到360
							c_alpha=c_Radian_Start:pi/100:2*pi;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
							%第2部分:0到Radian_End
							c_alpha=0.0:pi/100:c_Radian_End;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						end
					%**********************
					%   如果是顺时针圆弧
					%**********************
					elseif c_Direcction < (-0.5)
						%若结束角度大于起始角度,则分成两部分绘制圆弧
						if c_Radian_End>=c_Radian_Start 
							%第1部分:Radian_End到360
							c_alpha=c_Radian_End:pi/100:2*pi;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
							%第2部分:0到Radian_Start
							c_alpha=0.0:pi/100:c_Radian_Start;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						%若结束角度大于起始角度,则直接绘制圆弧
						else
							c_alpha=c_Radian_End:pi/100:c_Radian_Start;
							c_x=c_R*cos(c_alpha)+Arc_Crack_Coor(i,i_Seg,1);
							c_y=c_R*sin(c_alpha)+Arc_Crack_Coor(i,i_Seg,2);
							plot(c_x,c_y,'w','LineWidth',Width_Crack,'Color',Color_Crack)
						end

					end
				%如果是直线裂缝
				elseif abs(sum(Arc_Crack_Coor(i,i_Seg,1:11)))<1.0e-10 
					plot(x,y,'w','LineWidth',Width_Crack,'Color',Color_Crack)   
				end

			end
		end	
	end
end

% Plot circle holes.
if num_Hole ~=0
    disp(['      ----- Plotting circle holes...'])
	for i = 1:num_Hole
	    Coor_x  = Hole_Coor(i,1);
		Coor_y  = Hole_Coor(i,2);
		c_R  = Hole_Coor(i,3);
		alpha=0:pi/20:2*pi;
		x = Coor_x + c_R*cos(alpha);
		y = Coor_y + c_R*sin(alpha);
		plot(x,y,'-')
	end	
end

% Plot ellipse holes (2020-08-09).
if num_Ellipse_Hole ~=0
    disp(['      ----- Plotting ellipse holes...'])
	for i = 1:num_Ellipse_Hole
	    Coor_x  = Ellipse_Hole_Coor(i,1);
		Coor_y  = Ellipse_Hole_Coor(i,2);
		a_Hole  = Ellipse_Hole_Coor(i,3);
		b_Hole  = Ellipse_Hole_Coor(i,4);
		c_Hole  = sqrt(a_Hole^2-b_Hole^2);
		theta_Hole  = Ellipse_Hole_Coor(i,5);
		% Plot the ellipse.
		hold on;
		syms x y;
		theta_Hole = theta_Hole*pi/180.0;
		f=(a_Hole^2-c_Hole^2*cos(theta_Hole)^2)*(x-Coor_x).^2+(a_Hole^2-c_Hole^2*sin(theta_Hole)^2)*(y-Coor_y).^2-c_Hole^2*sin(2*theta_Hole)*(x-Coor_x).*(y-Coor_y)-a_Hole^2*b_Hole^2;
		h=ezplot(f,[-1.5*a_Hole+Coor_x,1.5*a_Hole+Coor_x,-1.5*a_Hole+Coor_y,1.5*a_Hole+Coor_y]);
		h.Color = 'black';
		title('Finite Element Mesh','FontName',Title_Font,'FontSize',Size_Font)
	end	
end

% 绘制圆形夹杂.
if num_Circ_Inclusion ~=0
    disp(['      ----- Plotting circle inclusions...'])
	for i = 1:num_Circ_Inclusion
	    Coor_x  = Circ_Inclu_Coor(i,1);
		Coor_y  = Circ_Inclu_Coor(i,2);
		c_R  = Circ_Inclu_Coor(i,3);
		alpha=0:pi/20:2*pi;
		x = Coor_x + c_R*cos(alpha);
		y = Coor_y + c_R*sin(alpha);
		%plot(x,y,'-')
		% plot(x_new,y_new,'-')
		patch(x,y,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
	end	
end


% 绘制多边形夹杂,2016-10-04
if num_Poly_Inclusion ~=0
    disp(['      ----- Plotting polygon inclusions...'])
	for i = 1:num_Poly_Inclusion
		nPt = size(Poly_Incl_Coor_x{i},2);
		for iPt = 1:nPt
			x(iPt) = Poly_Incl_Coor_x{i}(iPt);
			y(iPt) = Poly_Incl_Coor_y{i}(iPt);
		end
		patch(x,y,Color_Inclusion,'facealpha',0.3,'edgecolor','black','LineWidth',0.1)	 %透明度'facealpha'
	end	
end


% Plot enriched nodes of cracks.
if isempty(Post_Enriched_Nodes) ~= 1 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes...'])
    for i =1 :size(Post_Enriched_Nodes,2)
	    for j =1:Num_Node
		    x_node = Node_Coor(j,1);                                          
	        y_node = Node_Coor(j,2);                                          
		    if Post_Enriched_Nodes(j,i)==1     % Tip nodes
			    plot(x_node,y_node,'bs','Color','blue','MarkerSize',10)
			elseif Post_Enriched_Nodes(j,i)==2 % Heaviside nodes
			    plot(x_node,y_node,'bo','Color','blue','MarkerSize',8)
			elseif Post_Enriched_Nodes(j,i)==3 % Junction nodes
			    plot(x_node,y_node,'ko','MarkerSize',10+2,'Color','blue')
			elseif Post_Enriched_Nodes(j,i)==6 % Junction nodes of crack and hole
			    plot(x_node,y_node,'ko','MarkerSize',11,'Color','black')
			elseif Post_Enriched_Nodes(j,i)==4 % Cross nodes
			    plot(x_node,y_node,'^','MarkerSize',12,'Color','blue')
			end
		end
	end
end

% Plot enriched nodes of circle holes.
if num_Hole~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of hole...'])
	for i=1:num_Hole
	    for j =1:Num_Node
		    if Enriched_Node_Type_Hl(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',4,'Color','r')
			end
		end
		
	end
end

% Plot enriched nodes of ellipse holes(2020-08-09).
if num_Ellipse_Hole~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of ellipse hole...'])
	for i=1:num_Ellipse_Hole
	    for j =1:Num_Node
		    if Enriched_Node_Type_Hl(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',4,'Color','r')
			end
		end
		
	end
end

% Plot enriched nodes of crosses.
if num_Cross~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of cross...'])
	for i=1:num_Cross
	    for j =1:Num_Node
		    if Enriched_Node_Type_Cross(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',7,'Color','green')
			end
		end
		
	end
end

% Plot enriched nodes of crosses(五角星).
if num_HC~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of HC cross...'])
	for i=1:num_HC
	    for j =1:Num_Node
		    if Enriched_Node_Type_HC(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'p','MarkerSize',18,'Color','m')
			end
		end
		
	end
end

% 绘制圆形夹杂的增强节点.
if num_Circ_Inclusion~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of circle inclusions...'])
	for i=1:num_Circ_Inclusion
	    for j =1:Num_Node
		    if Enriched_Node_Type_Incl(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',4,'Color',[255/255,99/255,71/255])
			end
		end
		
	end
end

% 绘制多边形夹杂的增强节点.
if num_Poly_Inclusion~=0 && Key_PLOT(1,8) == 1
    disp(['      ----- Plotting enriched nodes of polygon inclusions...'])
	for i=1:num_Poly_Inclusion
	    for j =1:Num_Node
		    if Enriched_Node_Type_Incl(j,i) ==1
		        x_node = Node_Coor(j,1);                                          
	            y_node = Node_Coor(j,2);   
				plot(x_node,y_node,'bo','MarkerSize',6,'Color',[255/255,99/255,71/255])
			end
		end
		
	end
end

% Plot Gauss points.
if Key_PLOT(1,4) == 1
    if exist([Full_Pathname,'.gcor_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting Gauss points...'])
		% Read gauss point coordinates file.
		if Key_Data_Format==1 
			Gauss_Coor = load([Full_Pathname,'.gcor_',num2str(isub)]);
		elseif Key_Data_Format==2  %Binary
			c_file = fopen([Full_Pathname,'.gcor_',num2str(isub)],'rb');
			[cc_tem_Gauss_cor,cc_count]   = fread(c_file,inf,'double');
			fclose(c_file);
			%转换成Matlab中的数据格式
			for ccc_i=1:cc_count/2
				Gauss_Coor(ccc_i,1) = ccc_i;
				Gauss_Coor(ccc_i,2) = cc_tem_Gauss_cor(ccc_i*2-1);
				Gauss_Coor(ccc_i,3) = cc_tem_Gauss_cor(ccc_i*2);
			end
		end	
			
		
		plot(Gauss_Coor(:,2),Gauss_Coor(:,3),'bo','MarkerSize',1,'Color','black')

		clear Gauss_Coor
	end
end

% 绘制所有自由度的载荷.
if Key_PLOT(1,7)==1 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting forces of DOFs...'])
	%如果载荷文件存在则：
	if exist([Full_Pathname,'.fxdf_',num2str(Num_Step_to_Plot)], 'file') ==2  
		% Read x方向载荷.
		Force_x = load([Full_Pathname,'.fxdf_',num2str(isub)]);   
		% Read y方向载荷.
		Force_y = load([Full_Pathname,'.fydf_',num2str(isub)]); 
		Node_Force_x(1:Num_Node) = 0;
		Node_Force_y(1:Num_Node) = 0;
		% 节点循环,计算每个节点对应的增强节点号 
		for i_N =1:Num_Node
			for iCrack=1:num_Crack(isub)
				Num_Enriched_Node = POS(i_N,iCrack);
				if Num_Enriched_Node~=0 
					Node_Force_x(i_N)=Node_Force_x(i_N) + Force_x(Num_Enriched_Node);
					Node_Force_y(i_N)=Node_Force_y(i_N) + Force_y(Num_Enriched_Node);
				end
				% if Num_Enriched_Node~=0 %& (abs(Node_Force_x(i_N)) + abs(Node_Force_y(i_N))>10)
					% i_N,Force_x(Num_Enriched_Node),Force_y(Num_Enriched_Node)
				% end
			end
		end
		% Plot forces.
		Max_x_Force = max(abs(Node_Force_x));
		Max_y_Force = max(abs(Node_Force_y));
		%Min_x_Force = min(Node_Force_x)
		%Min_y_Force = min(Node_Force_y)
		Max_Force   = max(Max_x_Force,Max_y_Force);
		
		% Get the maximum and minimum value of the new coordinates of all nodes.
		Min_X_Coor = min(min(Node_Coor(:,1)));
		Max_X_Coor = max(max(Node_Coor(:,1)));
		Min_Y_Coor = min(min(Node_Coor(:,2)));
		Max_Y_Coor = max(max(Node_Coor(:,2)));	
		
		W = Max_X_Coor - Min_X_Coor;
		H = Max_Y_Coor - Min_Y_Coor;
			
		% length of force arrow         
		length_arrow = sqrt(aveg_area_ele)*1.8;
		% length_arrow = sqrt(aveg_area_ele)*0.8;
		% Loop through each node.
		for i = 1:Num_Node
			if Node_Force_x(i) ~=0 || Node_Force_y(i) ~=0          % If the nodes has force load, then:
				c_force_x   = Node_Force_x(i);
				c_force_y   = Node_Force_y(i);
				delta_L_x = c_force_x*length_arrow/Max_Force;
				delta_L_y = c_force_y*length_arrow/Max_Force;
				
				StartPoint = [Node_Coor(i,1)-delta_L_x   Node_Coor(i,2)-delta_L_y     0];
				EndPoint   = [Node_Coor(i,1)             Node_Coor(i,2)               0];

				line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','red','LineWidth',2)
				
				% The length of the head of the arrow.
				length_arrow_head = length_arrow/3;
				
				% Plot the head of the arrow.
				theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
				theta_1 = pi/2 - theta - pi/3;
				delta_x = -length_arrow_head*cos(theta_1);
				delta_y =  length_arrow_head*sin(theta_1);
				line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',1.5);
				theta_2 = 3*pi/2 - theta + pi/3;
				delta_x = -length_arrow_head*cos(theta_2);
				delta_y =  length_arrow_head*sin(theta_2);
				line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','red','LineWidth',1.5);
				
			end	
		end
	end
end

% 绘制计算点水压.
if Key_PLOT(1,7)==2 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting pressure of fluid nodes...'])
	%如果水压文件存在则：
	if exist([Full_Pathname,'.cpre_',num2str(Num_Step_to_Plot)], 'file') ==2
        %读取水压文件	
		namefile= [Full_Pathname,'.cpre_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Pressure(lineNum,1:c_num)  = str2num(TemData);   
			average_pres=sum(Pressure(lineNum,1:c_num))/c_num;  %平均水压
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%水压最大值
		Max_Pressure = max(max(abs(Pressure)));
		%读取计算点x坐标
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%读取计算点y坐标
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
		%读取计算点方位文件
		namefile= [Full_Pathname,'.cori_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Orient(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
        %绘制水压线		
		Max_Length = sqrt(aveg_area_ele)*4.5;  %水压线的最大长度
		for i_C = 1:num_Crack(Num_Step_to_Plot)
			for i_CalP=1:num_CalP_each_Crack(i_C)
                c_Pres_L = Pressure(i_C,i_CalP)/Max_Pressure*Max_Length;  %当前计算点水压线长度
				c_Pres_x = x(i_C,i_CalP);                                 %当前计算点x坐标
				c_Pres_y = y(i_C,i_CalP);                                 %当前计算点y坐标
				c_Orient = Orient(i_C,i_CalP);                            %当前计算点方位
				%若长度非0则绘制
				if c_Pres_L ~= 0
				    %得到起点坐标
					c_Start_x = c_Pres_x + c_Pres_L*0.5*sin(c_Orient);
					c_Start_y = c_Pres_y - c_Pres_L*0.5*cos(c_Orient);
					%得到终点坐标
					c_End_x   = c_Pres_x - c_Pres_L*0.5*sin(c_Orient);
					c_End_y   = c_Pres_y + c_Pres_L*0.5*cos(c_Orient);
					%绘制
					if c_Pres_L > 0     %正水压蓝色或翠绿色
					    %line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','blue','LineWidth',2.5);  %蓝色
                        line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color',[0/255,210/255,87/255],'LineWidth',2.5); %翠绿色					
					elseif c_Pres_L < 0 %负水压红色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','red','LineWidth',2.5);
					end
				end
			end
		end
	end
end

% 绘制计算点流量.
if Key_PLOT(1,7)==3 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting flux of fluid nodes...'])
	%如果流量文件存在则：
	if exist([Full_Pathname,'.cqua_',num2str(Num_Step_to_Plot)], 'file') ==2
        %读取流量文件	
		namefile= [Full_Pathname,'.cqua_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Quantity(lineNum,1:c_num)  = str2num(TemData);   
			average_pres=sum(Quantity(lineNum,1:c_num))/c_num;  %平均水压
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%流量最大值
		Max_Quantity = max(max(abs(Quantity)));
		%读取计算点x坐标
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%读取计算点y坐标
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
		%读取计算点方位文件
		namefile= [Full_Pathname,'.cori_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Orient(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
        %绘制流量线		
		Max_Length = sqrt(aveg_area_ele)*4.5;  %流量线的最大长度
		for i_C = 1:num_Crack(Num_Step_to_Plot)
			for i_CalP=1:num_CalP_each_Crack(i_C)
                c_Pres_L = Quantity(i_C,i_CalP)/Max_Quantity*Max_Length;  %当前计算点水压线长度
				c_Pres_x = x(i_C,i_CalP);                                 %当前计算点x坐标
				c_Pres_y = y(i_C,i_CalP);                                 %当前计算点y坐标
				c_Orient = Orient(i_C,i_CalP);                            %当前计算点方位
				%若长度非0则绘制
				if c_Pres_L ~= 0
				    %得到起点坐标
					c_Start_x = c_Pres_x + c_Pres_L*0.5*sin(c_Orient);
					c_Start_y = c_Pres_y - c_Pres_L*0.5*cos(c_Orient);
					%得到终点坐标
					c_End_x   = c_Pres_x - c_Pres_L*0.5*sin(c_Orient);
					c_End_y   = c_Pres_y + c_Pres_L*0.5*cos(c_Orient);
					%绘制
					if c_Pres_L > 0     %正流量蓝色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','blue','LineWidth',2.5);  %蓝色
                        %line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color',[0/255,210/255,87/255],'LineWidth',2.5); %翠绿色					
					elseif c_Pres_L < 0 %负流量红色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','red','LineWidth',2.5);
					end
				end
			end
		end
	end
end

% 绘制计算点开度.
if Key_PLOT(1,7)==4 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting aperture of calculation point...'])
	%如果开度文件存在则：
	if exist([Full_Pathname,'.cape_',num2str(Num_Step_to_Plot)], 'file') ==2
        %读取开度文件	
		namefile= [Full_Pathname,'.cape_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Aperture(lineNum,1:c_num)  = str2num(TemData);   
			average_Aper=sum(Aperture(lineNum,1:c_num))/c_num;  %平均开度
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%开度最大值
		Max_Aper = max(max(abs(Aperture)));
		% Max_Aper = 7.84*1e-3  %用于paper2-算例三
		% Max_Aper = 3.7237*1e-3%用于paper2-最后新增的验证算例4.3
		
		%读取计算点x坐标
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%读取计算点y坐标
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
		%读取计算点方位文件
		namefile= [Full_Pathname,'.cori_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Orient(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
        %绘制开度曲线			
		% Max_Length = sqrt(aveg_area_ele)*7.5;  %开度线的最大长度
		Max_Length = sqrt(aveg_area_ele)*3.5;  %开度线的最大长度
		for i_C = 1:num_Crack(Num_Step_to_Plot)
			for i_CalP=1:num_CalP_each_Crack(i_C)
                c_Aper_L = Aperture(i_C,i_CalP)/Max_Aper*Max_Length;  %当前计算点开度线长度
				c_Aper_x = x(i_C,i_CalP);                                 %当前计算点x坐标
				c_Aper_y = y(i_C,i_CalP);                                 %当前计算点y坐标
				c_Orient = Orient(i_C,i_CalP);                            %当前计算点方位
				%若长度非0则绘制
				if c_Aper_L ~= 0
				    %得到起点坐标
					c_Start_x = c_Aper_x + c_Aper_L*0.5*sin(c_Orient);
					c_Start_y = c_Aper_y - c_Aper_L*0.5*cos(c_Orient);
					%得到终点坐标
					c_End_x   = c_Aper_x - c_Aper_L*0.5*sin(c_Orient);
					c_End_y   = c_Aper_y + c_Aper_L*0.5*cos(c_Orient);
					%绘制
					if c_Aper_L > 0     %正开度蓝色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','blue','LineWidth',1.5);  %蓝色
                        %line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color',[0/255,210/255,87/255],'LineWidth',2.5); %翠绿色					
					elseif c_Aper_L < 0 %负开度红色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','red','LineWidth',1.5);
					end
				end
			end
		end
	end
end

% 绘制计算点粘聚裂缝粘聚力_x方向,2017-04-24.
if Key_PLOT(1,7)==5 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting traction x of calculation point of cohesive crack...'])
	%如果文件存在则：
	if exist([Full_Pathname,'.cohx_',num2str(Num_Step_to_Plot)], 'file') ==2
        %读取文件	
		namefile= [Full_Pathname,'.cohx_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Aperture(lineNum,1:c_num)  = str2num(TemData);   
			average_Aper=sum(Aperture(lineNum,1:c_num))/c_num;  %平均开度
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%牵引力最大值
		Max_Traction = max(max(abs(Aperture)));
		% Max_Aper = 7.84*1e-3  %用于paper2-算例三
		% Max_Aper = 3.7237*1e-3%用于paper2-最后新增的验证算例4.3
		
		%读取计算点x坐标
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%读取计算点y坐标
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
        %绘制牵引力			
		% Max_Length = sqrt(aveg_area_ele)*7.5;  %开度线的最大长度
		Max_Length = sqrt(aveg_area_ele)*3.5;  %开度线的最大长度
		for i_C = 1:num_Crack(Num_Step_to_Plot)
			for i_CalP=1:num_CalP_each_Crack(i_C)
                c_Aper_L = Aperture(i_C,i_CalP)/Max_Traction*Max_Length;  %当前计算点开度线长度
				c_Aper_x = x(i_C,i_CalP);                                 %当前计算点x坐标
				c_Aper_y = y(i_C,i_CalP);                                 %当前计算点y坐标
				%若长度非0则绘制
				if c_Aper_L ~= 0
				    %得到起点坐标
					c_Start_x = c_Aper_x+ c_Aper_L*0.5;
					c_Start_y = c_Aper_y;
					%得到终点坐标
					c_End_x   = c_Aper_x- c_Aper_L*0.5;
					c_End_y   = c_Aper_y ;
					%绘制
					if c_Aper_L > 0     %正值蓝色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','blue','LineWidth',1.5);  %蓝色
                        %line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color',[0/255,210/255,87/255],'LineWidth',2.5); %翠绿色					
					elseif c_Aper_L < 0 %负值红色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','red','LineWidth',1.5);
					end
				end
			end
		end
	end
end


% 绘制计算点粘聚裂缝粘聚力_y方向,2017-04-24.
if Key_PLOT(1,7)==6 && num_Crack(Num_Step_to_Plot)~=0
    disp(['      ----- Plotting traction y of calculation point of cohesive crack...'])
	%如果文件存在则：
	if exist([Full_Pathname,'.cohy_',num2str(Num_Step_to_Plot)], 'file') ==2
        %读取文件	
		namefile= [Full_Pathname,'.cohy_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			Aperture(lineNum,1:c_num)  = str2num(TemData);   
			average_Aper=sum(Aperture(lineNum,1:c_num))/c_num;  %平均开度
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%牵引力最大值
		Max_Traction = max(max(abs(Aperture)));
		% Max_Aper = 7.84*1e-3  %用于paper2-算例三
		% Max_Aper = 3.7237*1e-3%用于paper2-最后新增的验证算例4.3
		
		%读取计算点x坐标
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
		%读取计算点y坐标
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
        %绘制牵引力			
		% Max_Length = sqrt(aveg_area_ele)*7.5;  %开度线的最大长度
		Max_Length = sqrt(aveg_area_ele)*3.5;  %开度线的最大长度
		for i_C = 1:num_Crack(Num_Step_to_Plot)
			for i_CalP=1:num_CalP_each_Crack(i_C)
                c_Aper_L = Aperture(i_C,i_CalP)/Max_Traction*Max_Length;  %当前计算点开度线长度
				c_Aper_x = x(i_C,i_CalP);                                 %当前计算点x坐标
				c_Aper_y = y(i_C,i_CalP);                                 %当前计算点y坐标
				%若长度非0则绘制
				if c_Aper_L ~= 0
				    %得到起点坐标
					c_Start_x = c_Aper_x;
					c_Start_y = c_Aper_y + c_Aper_L*0.5;
					%得到终点坐标
					c_End_x   = c_Aper_x;
					c_End_y   = c_Aper_y - c_Aper_L*0.5;
					%绘制
					if c_Aper_L > 0     %正值蓝色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','blue','LineWidth',1.5);  %蓝色
                        %line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color',[0/255,210/255,87/255],'LineWidth',2.5); %翠绿色					
					elseif c_Aper_L < 0 %负值红色
					    line([c_Start_x,c_End_x],[c_Start_y,c_End_y],'color','red','LineWidth',1.5);
					end
				end
			end
		end
	end
end
		
% 绘制HF水力压裂分析计算点(流体节点).
if (Key_PLOT(1,5)==1 || Key_PLOT(1,5)==2) && num_Crack(Num_Step_to_Plot)~=0 && (Key_PLOT(1,6)==1||Key_PLOT(1,6)==2)
    disp(['      ----- Plotting calculation points of HF...'])
	%读取计算点x坐标
	if exist([Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)], 'file') ==2
		namefile= [Full_Pathname,'.apex_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			x(lineNum,1:c_num)  = str2num(TemData); 
			num_CalP_each_Crack(lineNum) = c_num;  %每条裂纹的计算点数目
		end
		fclose(data); 
	end

	%读取计算点y坐标
	if exist([Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)], 'file') ==2
		namefile= [Full_Pathname,'.apey_',num2str(Num_Step_to_Plot)];
		data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);              
			c_num   = size(str2num(TemData),2); 
			y(lineNum,1:c_num)  = str2num(TemData); 
		end
		fclose(data); 
	end
	
	%绘制计算点及其编号
	CalP_num = 0;
    for i=1:num_Crack(Num_Step_to_Plot)
	    for j=1:num_CalP_each_Crack(i)
		    CalP_num = CalP_num +1;
	        plot(x(i,j),y(i,j),'bo','Color',[218/255,112/255,214/255],'MarkerSize',6,'linewidth',1.5)
			if Key_PLOT(1,6)==2
			    text(x(i,j)+0.10*delta,y(i,j),1,num2str(CalP_num),...
		              'FontName',Title_Font,'FontSize',11,'color',[160/255,32/255,240/255])
			end
		end
	end
end

% 单元应力状态(Stress1-Stress3<Tol)
if Key_PLOT(1,11)==1 
    Yes_Fractured=0;  
	%如果单元应力状态文件存在则：
	if exist([Full_Pathname,'.elss_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting stress status of elements...'])
		ELSS= load([Full_Pathname,'.elss_',num2str(Num_Step_to_Plot)]);
		for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4)];                             % Nodes for current element
			if ELSS(iElem) == 1    %单元应力状态是否σ1-σ3 < Tol
				xi_Elcs1(:,iElem) = Node_Coor(NN',1);                                     
				yi_Elcs1(:,iElem) = Node_Coor(NN',2); 	
				Yes_Fractured=1;
			end
		end
		if Yes_Fractured==1
			patch(xi_Elcs1,yi_Elcs1,[255/255,255/255,0/255])       %用黄色标记该单元
		end
	end
end

% 绘制裂缝编号(不含天然裂缝).
if Key_PLOT(1,14) == 1
    disp(['      ----- Plotting crack number...'])
	if isempty(Crack_X)==0
		for i = 1:num_Crack(isub)
			nPt = size(Crack_X{i},2);
			%各裂缝的相对而言中心点的坐标
            ave_x = sum(Crack_X{i}(1:nPt))/nPt;
			ave_y = sum(Crack_Y{i}(1:nPt))/nPt;
			ave_x_offset = ave_x + delta*0.5;
			ave_y_offset = ave_y + delta*0.5;
	        plot(ave_x_offset,ave_y_offset,'bo','Color','black','MarkerSize',15,'linewidth',2)
			text(ave_x_offset,ave_y_offset,1,num2str(i),...
		              'FontName',Title_Font,'FontSize',10,'color','black','FontWeight','bold')
		end	
	end
end

% 绘制支撑剂的圆球
if Key_PLOT(1,10)==1 
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
			R = c_D_Coors(i,2)/2*Key_PLOT(2,6);
			old_Coor_x = c_D_Coors(i,3);
			old_Coor_y = c_D_Coors(i,4);
			Coor_x = old_Coor_x; 
			Coor_y = old_Coor_y; 
			
			x = Coor_x + R*cos(alpha);
			y = Coor_y + R*sin(alpha);
			
			plot(x,y,'-')
			%axis equal
			fill(x,y,[128/255,138/255,135/255])
		end
	end
end

% 绘制破裂区
if Key_PLOT(1,15)==1 
	%如果定义了破裂区
	if Yes_has_FZ ==1
	    disp(['      ----- Plotting fracture zone......'])  
        c_x = [frac_zone_min_x,frac_zone_max_x];
		c_y = [frac_zone_min_y,frac_zone_min_y];
		line(c_x,c_y,'LineWidth',1.5,'Color','green')
		% line(c_x,c_y,'LineWidth',3.0,'Color','black')
        c_x = [frac_zone_max_x,frac_zone_max_x];
		c_y = [frac_zone_min_y,frac_zone_max_y];
		line(c_x,c_y,'LineWidth',1.5,'Color','green')
		% line(c_x,c_y,'LineWidth',3.0,'Color','black')
        c_x = [frac_zone_max_x,frac_zone_min_x];
		c_y = [frac_zone_max_y,frac_zone_max_y];
		line(c_x,c_y,'LineWidth',1.5,'Color','green')
		% line(c_x,c_y,'LineWidth',3.0,'Color','black')
        c_x = [frac_zone_min_x,frac_zone_min_x];
		c_y = [frac_zone_max_y,frac_zone_min_y];
		line(c_x,c_y,'LineWidth',1.5,'Color','green')
		% line(c_x,c_y,'LineWidth',3.0,'Color','black')
	end
end

% 场问题相关绘制,2016-11-08
if Yes_Field_Problem==1 
    %边值边界条件
    for i=1:size(Field_Boundary_Value,1)
	    c_node = Field_Boundary_Value(i,1);
		x_node = Node_Coor(c_node,1);                                          
		y_node = Node_Coor(c_node,2);    
		plot(x_node,y_node,'ko','MarkerSize',10,'Color','blue')
	end
    %流量边界条件
    for i=1:size(Field_Boundary_Qn,1)
	    c_node = Field_Boundary_Qn(i,1);
		x_node = Node_Coor(c_node,1);                                          
		y_node = Node_Coor(c_node,2);    
		plot(x_node,y_node,'ko','MarkerSize',10,'Color','red')
	end
end

%测试用,2021-06-27
% TEST_NODES = load([Full_Pathname,'.TEST1']);
% TEST_ELES  = load([Full_Pathname,'.TEST2']);
% plot(TEST_NODES(:,1),TEST_NODES(:,2),'bo','MarkerSize',3,'Color','black')
% num_Added_Eles =  size(TEST_ELES,1);
% for i=1:num_Added_Eles
    % node_1 = TEST_ELES(i,1);
	% node_2 = TEST_ELES(i,2);
	% node_3 = TEST_ELES(i,3);
	% node_4 = TEST_ELES(i,4);
	% if node_1 <= Num_Node
	    % x_Coor(1) = Node_Coor(node_1,1);
		% y_Coor(1) = Node_Coor(node_1,2);
	% else
	    % x_Coor(1) = TEST_NODES(node_1,1);
		% y_Coor(1) = TEST_NODES(node_1,2);
	% end
	% if node_2 <= Num_Node
	    % x_Coor(2) = Node_Coor(node_2,1);
		% y_Coor(2) = Node_Coor(node_2,2);
	% else
	    % x_Coor(2) = TEST_NODES(node_2,1);
		% y_Coor(2) = TEST_NODES(node_2,2);
	% end
	% if node_3 <= Num_Node
	    % x_Coor(3) = Node_Coor(node_3,1);
		% y_Coor(3) = Node_Coor(node_3,2);
	% else
	    % x_Coor(3) = TEST_NODES(node_3,1);
		% y_Coor(3) = TEST_NODES(node_3,2);
	% end
	% if node_4 <= Num_Node
	    % x_Coor(4) = Node_Coor(node_4,1);
		% y_Coor(4) = Node_Coor(node_4,2);
	% else
	    % x_Coor(4) = TEST_NODES(node_4,1);
		% y_Coor(4) = TEST_NODES(node_4,2);
	% end	
    % patch(x_Coor,y_Coor,'red'); 
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
