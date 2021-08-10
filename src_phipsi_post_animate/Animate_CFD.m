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

function Animate_CFD(Real_num_iteration)
% 绘制CFD云图.

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
global Mesh_Points number_mesh_points number_mesh_connect Mesh_Connects CFD_Contour_fine
global Title_Font Key_Figure_Control_Widget

disp('    Generating animations of CFD......')

% Get the maximum and minimum value of displacements.
disp(['    > Calculating the range of CFD value......']) 


% Get the max and min value of field value.
if Key_Animation(1)==1  & Key_Ani_Ave(1)==1
    i_output=0;
	for i=1:Real_num_iteration
	    if mod(i,Output_Freq)==0
		    i_output=i_output+1;
			if Key_Data_Format==1 
				CFD_Value = load([Full_Pathname,'.pres_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.pres_',num2str(Itera_Num(i))],'rb');
				[CFD_Value,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
			end
	
			if i_output==1
				MaxValue =max(CFD_Value(1:number_mesh_points));
				MinValue =min(CFD_Value(1:number_mesh_points));
			elseif i_output > 1
				if max(CFD_Value(1:number_mesh_points)) > MaxValue, MaxValue=max(CFD_Value(1:number_mesh_points)); end
				if min(CFD_Value(1:number_mesh_points)) < MinValue, MinValue=min(CFD_Value(1:number_mesh_points)); end
			end
			clear CFD_Value
		end
	end
end



% Loop through each frame to plot and store.
i_output=0;
for i=1:Real_num_iteration
    if mod(i,Output_Freq)==0
	    i_output = i_output + 1;
        % Read displacement file.
        if exist([Full_Pathname,'.pres_',num2str(Itera_Num(i))], 'file') ==2
			if Key_Data_Format==1 
				CFD_Value = load([Full_Pathname,'.pres_',num2str(Itera_Num(i))]);
			elseif Key_Data_Format==2  %Binary
				c_file = fopen([Full_Pathname,'.pres_',num2str(Itera_Num(i))],'rb');
				[CFD_Value,cc_count]   = fread(c_file,inf,'double');
				fclose(c_file);
			end
			disp(['    > Plotting and saving CFD value contours of frame ',num2str(i),'......']) 
		else
		    disp(['    WARNING :: can not find *.pres file']) 
			break
		end
		scale = 1.0;
		
		% Get resample coors.
		%delta = sqrt(aveg_area_ele)/Num_Accuracy_Contour;
		delta = ((Max_X_Coor - Min_X_Coor)+(Max_Y_Coor-Min_Y_Coor))/2/CFD_Contour_fine;
		gx    = Min_X_Coor:delta:Max_X_Coor; 
		gy    = Min_Y_Coor:delta:Max_Y_Coor;
		
		% Plot field value. 
		Tools_New_Figure
		hold on;

		disp('      Resample field value....') 
		[c_CFD_value,X,Y] = Tools_gridfit(Mesh_Points(:,2),Mesh_Points(:,3),CFD_Value(1:number_mesh_points),gx,gy);
		contourf(X,Y,c_CFD_value,Num_Contourlevel,'LineStyle','none')
		
		
		% Set colormap.
		if Key_Flipped_Gray==0
			colormap(Color_Contourlevel)
		elseif Key_Flipped_Gray==1
			colormap(flipud(gray))
		end
		title('CFD pressure','FontName',Title_Font,'FontSize',Size_Font)
		clear c_CFD_value
			
		% Plot the mesh		
		if Key_PLOT(2,3)==1		
			for iconnect =1:number_mesh_connect
				NN = [Mesh_Connects(iconnect,2) Mesh_Connects(iconnect,3) Mesh_Connects(iconnect,4) Mesh_Connects(iconnect,2)];   
				xi = Mesh_Points(NN',2);                                           % Deformed x-coordinates of nodes
				yi = Mesh_Points(NN',3);                                           % Deformed y-coordinates of nodes
				plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
			end
		end
			
		% The range of the plot.
		delta = 0;
		min_x = Min_X_Coor-delta; 
		max_x = Max_X_Coor+delta; 
		min_y = Min_Y_Coor-delta; 
		max_y = Max_Y_Coor+delta;
		axis([min_x max_x min_y max_y]);
		
		% Plot the outline of the mesh
		line([Mesh_Points(Outline(:,1),2) Mesh_Points(Outline(:,2),2)]', ...
			 [Mesh_Points(Outline(:,1),3) Mesh_Points(Outline(:,2),3)]','LineWidth',0.5,'Color','black')		
		%外边界填充
		Tools_Fillout(Mesh_Points(Outline(:,1),2),Mesh_Points(Outline(:,1),3),[Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor],'w'); 
		%内边界填充.
		if isempty(Inline)==0
			fill(Mesh_Points(Inline(:,1),2),Mesh_Points(Inline(:,1),3),'w');
		end	
				
	
		axis off; 
		axis equal;
		colorbar('FontName',Title_Font,'FontSize',Size_Font);
			
		% Control of the range of the colour bar.
		if Key_Ani_Ave(1)==1
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





