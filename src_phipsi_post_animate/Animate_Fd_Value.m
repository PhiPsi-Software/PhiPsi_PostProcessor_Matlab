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

function Animate_Fd_Value(Real_num_iteration)
% 绘制场问题场变量云图.

global Node_Coor Elem_Node Outline Inline
global Num_Elem Version Num_Node
global Key_PLOT Full_Pathname Key_Dynamic
global Size_Font num_Crack
global aveg_area_ele Time_Delay delt_time_NewMark Num_Contourlevel
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor
global Key_Ani_Ave Width_Crack Color_Crack Key_Animation
global Num_Gauss_Points Key_Integral_Sol
global Num_Accuracy_Contour Key_Contour_Metd
global Output_Freq Color_Mesh
global Color_Contourlevel Key_Flipped_Gray Itera_Num Itera_HF_Time
global Na_Crack_X Na_Crack_Y num_Na_Crack Key_HF_Analysis
global frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global Num_Foc_x Num_Foc_y  Foc_x Foc_y
global num_Hole Hole_Coor Enriched_Node_Type_Hl POS_Hl Elem_Type_Hl
global num_Circ_Inclusion Circ_Inclu_Coor Enriched_Node_Type_Incl POS_Incl Elem_Type_Incl
global num_Poly_Inclusion Poly_Incl_Coor_x Poly_Incl_Coor_y
global Field_Flux_x Field_Flux_y
global max_area_ele
global Max_Flux_x Min_Flux_x Max_Flux_y Min_Flux_y Max_Flux Key_Time_String
global Key_Data_Format
global Title_Font Key_Figure_Control_Widget

disp('    Generating animations of field value......')

% Get the maximum and minimum value of displacements.
disp(['    > Calculating the range of the deformed body......']) 


% Get the max and min value of field value.
if (Key_Animation(4)==1 | Key_Animation(4)==3) & Key_Ani_Ave(4)==1
    i_output=0;
	for i=1:Real_num_iteration
	    if mod(i,Output_Freq)==0
		    i_output=i_output+1;
			if Key_Data_Format==1 
				Field_Value = load([Full_Pathname,'.fdvl_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.fdvl_',num2str(Itera_Num(i))],'rb');
				[Field_Value,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
			end
	
			if i_output==1
				MaxValue =max(Field_Value(1:Num_Node));
				MinValue =min(Field_Value(1:Num_Node));
			elseif i_output > 1
				if max(Field_Value(1:Num_Node)) > MaxValue, MaxValue=max(Field_Value(1:Num_Node)); end
				if min(Field_Value(1:Num_Node)) < MinValue, MinValue=min(Field_Value(1:Num_Node)); end
			end
			clear Field_Value
		end
	end
end

% 获得流量的最大值和最小值.
if (Key_Animation(4)==2 | Key_Animation(4)==3)  & Key_Ani_Ave(4)==1
    i_output=0;
	for i=1:Real_num_iteration
	    if mod(i,Output_Freq)==0
		    i_output=i_output+1;
			Field_Flux_x = load([Full_Pathname,'.fbfx_',num2str(Itera_Num(i))]);
			if i_output==1
				Max_Flux_x =max(Field_Flux_x(1:Num_Node));
				Min_Flux_x =min(Field_Flux_x(1:Num_Node));
			elseif i_output > 1
				if max(Field_Flux_x(1:Num_Node)) > Max_Flux_x, Max_Flux_x=max(Field_Flux_x(1:Num_Node)); end
				if min(Field_Flux_x(1:Num_Node)) < Min_Flux_x, Min_Flux_x=min(Field_Flux_x(1:Num_Node)); end
			end
			clear Field_Flux_x
		end
	end
    i_output=0;
	for i=1:Real_num_iteration
	    if mod(i,Output_Freq)==0
		    i_output=i_output+1;
			Field_Flux_y = load([Full_Pathname,'.fbfy_',num2str(Itera_Num(i))]);
			if i_output==1
				Max_Flux_y =max(Field_Flux_y(1:Num_Node));
				Min_Flux_y =min(Field_Flux_y(1:Num_Node));
			elseif i_output > 1
				if max(Field_Flux_y(1:Num_Node)) > Max_Flux_y, Max_Flux_y=max(Field_Flux_y(1:Num_Node)); end
				if min(Field_Flux_y(1:Num_Node)) < Min_Flux_y, Min_Flux_y=min(Field_Flux_y(1:Num_Node)); end
			end
			clear Field_Flux_x
		end
	end
	Max_Flux   = max(abs(Max_Flux_x),abs(Max_Flux_y));
end

% Loop through each frame to plot and store.
i_output=0;
for i=1:Real_num_iteration
    if mod(i,Output_Freq)==0
	    i_output = i_output + 1;
        % Read displacement file.
        if exist([Full_Pathname,'.fdvl_',num2str(Itera_Num(i))], 'file') ==2
			if Key_Data_Format==1 
				Field_Value = load([Full_Pathname,'.fdvl_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.fdvl_',num2str(Itera_Num(i))],'rb');
				[Field_Value,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
			end
			disp(['    > Plotting and saving field value contours of frame ',num2str(i),'......']) 
		else
		    disp(['    WARNING :: can not find *.fdcl file']) 
			break
		end
		scale = Key_PLOT(5,6);
	
		% Get resample coors.
		delta = sqrt(aveg_area_ele)/Num_Accuracy_Contour;
		gx    = Min_X_Coor:delta:Max_X_Coor; 
		gy    = Min_Y_Coor:delta:Max_Y_Coor;
		
		% Plot field value. 
		Tools_New_Figure
		hold on;
		if Key_Animation(4)==1 | Key_Animation(4)==3
			disp('      Resample field value....') 
			[c_Field_value,X,Y] = Tools_gridfit(Node_Coor(:,1),Node_Coor(:,2),Field_Value(1:Num_Node),gx,gy);
			contourf(X,Y,c_Field_value,Num_Contourlevel,'LineStyle','none')
			%Set colormap.
			%colormap(gray)
			%colormap(hot)
			colormap(cool)
			%colormap(jet)
		end
		
		%绘制流量矢量.
		if Key_Animation(4)==2 | Key_Animation(4)==3
		    Field_Flux_x = load([Full_Pathname,'.fbfx_',num2str(Itera_Num(i))]);
			Field_Flux_y = load([Full_Pathname,'.fbfy_',num2str(Itera_Num(i))]);
			W = Max_X_Coor - Min_X_Coor;
			H = Max_Y_Coor - Min_Y_Coor;
			% length of force arrow
			%length_arrow = sqrt(max_area_ele)*1.0;
			%length_arrow = max(W,H)/15.0;   
            length_arrow = (W+H)/2/20.0;     			
			% Loop through each node.
			for i = 1:Num_Node
				c_flux_x   = Field_Flux_x(i);
				c_flux_y   = Field_Flux_y(i);
				c_max_flux = max(abs(c_flux_x),abs(c_flux_y));
				delta_L_x  = c_flux_x*length_arrow/Max_Flux;
				delta_L_y  = c_flux_y*length_arrow/Max_Flux;
				StartPoint = [Node_Coor(i,1)-delta_L_x   Node_Coor(i,2)-delta_L_y     0];
				EndPoint   = [Node_Coor(i,1)             Node_Coor(i,2)               0];
                
				%流量达到一定范围才显示(超过最大值的千分之一)
				 c_Tol = 0.0;
				%c_Tol = Max_Flux*0.0001;
				if(c_max_flux > c_Tol)
				    line([StartPoint(1) EndPoint(1)],[StartPoint(2) EndPoint(2)],'color','black','LineWidth',1.5)
				end
				% The length of the head of the arrow.
				length_arrow_head = length_arrow/4;
				
				% Plot the head of the arrow.
				theta = atan2(EndPoint(2)-StartPoint(2),EndPoint(1)-StartPoint(1));
				theta_1 = pi/2 - theta - pi/3;
				delta_x = -length_arrow_head*cos(theta_1);
				delta_y =  length_arrow_head*sin(theta_1);
				%流量达到一定范围才显示(超过最大值的千分之一)
				if(c_max_flux > c_Tol)
				    line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','black','LineWidth',1.5);
				end
				theta_2 = 3*pi/2 - theta + pi/3;
				delta_x = -length_arrow_head*cos(theta_2);
				delta_y =  length_arrow_head*sin(theta_2);
				%流量达到一定范围才显示(超过最大值的千分之一)
				if(c_max_flux > c_Tol)
				    line([EndPoint(1) EndPoint(1)+delta_x],[EndPoint(2) EndPoint(2)+delta_y],'color','black','LineWidth',1.5);
				end
			end
			set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
		end
		
		% Set colormap.
		% if Key_Flipped_Gray==0
			% colormap(Color_Contourlevel)
		% elseif Key_Flipped_Gray==1
			% colormap(flipud(gray))
		% end
		title('Field value','FontName',Title_Font,'FontSize',Size_Font)
		clear c_Field_value
			
		% Plot the mesh.
		if Key_PLOT(5,9) ==1
			for iElem =1:Num_Elem
				NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				      Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)]; % Nodes for current element
				xi = Node_Coor(NN',1);                                           % Deformed x-coordinates of nodes
				yi = Node_Coor(NN',2);                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.1,'Color',Color_Mesh)
			end
		   end
			
		% The range of the plot.
		delta = sqrt(aveg_area_ele);
		min_x = Min_X_Coor-delta; max_x = Max_X_Coor+delta; 
		min_y = Min_Y_Coor-delta; max_y = Max_Y_Coor+delta;
		axis([min_x max_x min_y max_y]);
		
		% Fill the outside of the domain by white color.
        Tools_Fillout(Node_Coor(Outline(:,1),1),Node_Coor(Outline(:,1),2),[min_x max_x min_y max_y],'w'); 
		
		% Fill the inside of the domain by white color, holes, for example.
		if isempty(Inline)==0
            fill(Node_Coor(Inline(:,1),1),Node_Coor(Inline(:,1),2),'w');
		end	
	
		axis off; 
		axis equal;
		%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
		colorbar('FontName',Title_Font,'FontSize',Size_Font);
			
		% Control of the range of the colour bar.
		if Key_Ani_Ave(4)==1
            caxis([MinValue, MaxValue]); 
		end

		% Plot text on the left or the bottom of the figure.
		plot_string_MatMES = ['PhiPsi ',Version];
		plot_string_Frame  = ['Frame ',num2str(i_output),' / ',num2str(Real_num_iteration)];
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
		plot_string_Scale  = ['Scale factor: ',num2str(scale)];
		
		range_W = abs(Max_X_Coor-Min_X_Coor);
		range_H = abs(Max_Y_Coor-Min_Y_Coor);
		if range_H >= 0.75*range_W           % Left
			loc_x = -range_H/2+ Min_X_Coor;
			loc_y =  Max_Y_Coor-range_H*0.05;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_y =  Max_Y_Coor-range_H*0.15;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_y =  Max_Y_Coor-range_H*0.25;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			loc_y =  Max_Y_Coor-range_H*0.35;
			text(loc_x,loc_y,plot_string_Time, 'color','black');
		else                                  % Bottom
			loc_y =  Min_Y_Coor - range_H*0.2;
			loc_x =  Min_X_Coor;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_x =  Min_X_Coor + range_W*0.25;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_x =  Min_X_Coor + range_W*0.45;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			loc_x =  Min_X_Coor + range_W*0.60;
			text(loc_x,loc_y,plot_string_Time, 'color','black');
		end
			
	    % Save the current figure.
	    PLOT_Field_value(i_output) = getframe(gcf);       
		im=frame2im(PLOT_Field_value(i_output));         
		[I,map]=rgb2ind(im,256);
		k=i_output-0;
		str1='_field_value';
		str2=Full_Pathname;
		FileName1 =[str2,str1,'.gif'];
		if k==1;
			imwrite(I,map,FileName1,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
		else
			imwrite(I,map,FileName1,'gif','WriteMode','append','DelayTime',Time_Delay);
		end
		close
	end
end





