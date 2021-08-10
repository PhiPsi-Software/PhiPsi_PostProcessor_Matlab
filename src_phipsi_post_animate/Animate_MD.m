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

function Animate_MD
%分子动力学分析动画生成.

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
global Num_Step_to_Plot num_molecule
global Key_Data_Format
global Title_Font Key_Figure_Control_Widget

disp('    Generating MD animations......')

% 如果Num_Step_to_Plot = -999,则程序自动寻找最后一步的稳定计算结果并后处理
% if Num_Step_to_Plot == -999
    % for i_Check =1:5000
	    % if exist([Full_Pathname,'.mdco_',num2str(i_Check)], 'file') ==2 
	        % Num_Step_to_Plot = i_Check;
	    % end
	% end
% end
% if Num_Step_to_Plot==-999
	% disp(' >> Error :: No result files found, post-processor terminated.') 
	% Error_Message
% end




% 读取Box文件
if exist([Full_Pathname,'.mdbx'], 'file') ==2
    Yes_has_FZ = 1;
	disp('    > Reading mdbx file....') 
	namefile= [Full_Pathname,'.mdbx'];
	data=fopen(namefile,'r'); 
	lineNum =0;
	while ~feof(data)
		lineNum = lineNum+1;
		TemData = fgetl(data);    
		if lineNum==2   %第2行
			ttt_DATA(1:2)  = str2num(TemData);
		end
	end
	fclose(data); 
	box_zone_min_x = 0; 
	box_zone_max_x = ttt_DATA(1);
	box_zone_min_y = 0; 
	box_zone_max_y = ttt_DATA(2);
end
or_box_zone_min_x = box_zone_min_x; 
or_box_zone_max_x = box_zone_max_x;
or_box_zone_min_y = box_zone_min_y; 
or_box_zone_max_y = box_zone_max_y;

% 如果存在时间步文件
if  exist([Full_Pathname,'.hftm'], 'file') ==2 
	%读取水力压裂分析各破裂步对应的总的迭代次数
	disp('    > Reading hftm file....') 
	namefile= [Full_Pathname,'.hftm'];
	data=fopen(namefile,'r'); 
	lineNum = 0;
	num_Iter = 0;
	while ~feof(data)
		lineNum = lineNum+1;
		TemData = fgetl(data);    
		if lineNum>=2   %第一行是文件标识行,不予读取
			num_Iter = num_Iter+1;                     %总的迭代步数
			c_num   = size(str2num(TemData),2); 	   
			ttt_DATA(num_Iter,1:4)  = str2num(TemData);
		end
	end
	fclose(data); 
	%最大破裂步数
	Max_Frac = max(ttt_DATA(1:num_Iter,2));
	%提取每个破裂步对应的迭代步号
	for i_Fra = 1:Max_Frac
		Itera_Num(i_Fra) = 0;
		%所有破裂步间循环
		for i_ter = 1:num_Iter
			if ttt_DATA(i_ter,2)== i_Fra &&  ttt_DATA(i_ter,3) > Itera_Num(i_Fra)
				Itera_Num(i_Fra) = ttt_DATA(i_ter,3);
				Itera_HF_Time(i_Fra) = ttt_DATA(i_ter,4);
			end
		end
	end
% 如果不存在时间步文件
else
    for iii = 1:Num_Step_to_Plot
        Itera_Num(iii) = iii;
	end
end

num_Iter = num_Iter +1;
%Itera_Num
%-999
if Num_Step_to_Plot == -999
    Num_Step_to_Plot = Itera_Num(num_Iter);
end

%获取分子数目
MD_Coor = load([Full_Pathname,'.mdco_',num2str(Num_Step_to_Plot)]);
num_molecule = size(MD_Coor,1);
disp(['    Number of molecules:',num2str(num_molecule)]);

% 如果需要绘制轨迹线,则很耗时间,可以先把所有位置文件读入内存
if Key_PLOT(6,2) >= 1
    % all_MD_Coor=cell(num_Iter,)
    for ii=1:num_Iter
	    all_MD_Coor(ii,1:num_molecule,1:3) = load([Full_Pathname,'.mdco_',num2str(Itera_Num(ii))]);
    end
end
		
% Loop through each frame to plot and store.
i_output=0;
for c_i=1:num_Iter
	if mod(c_i,Output_Freq)==0
	    i_output = i_output + 1;

		scale = Key_PLOT(5,6);
		
		% New figure.
		Tools_New_Figure
		hold on;
		title('Position of particles','FontName',Title_Font,'FontSize',Size_Font)
		axis off; 
		axis equal;
        axis([box_zone_min_x box_zone_max_x box_zone_min_y box_zone_max_y]);

		%<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制分子位置及轨迹线.
		%<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(6,1) == 1
			% 轨迹线
			if Key_PLOT(6,2) >= 1
			    if c_i>=2
					disp(['      ----- Plotting moving paths of particles...'])
					tem1 = 1;
					tem2 = c_i-1;
					%如果Key_PLOT(6,2)=2,则只绘制最近 Key_PLOT(6,2)个计算步的轨迹
					if Key_PLOT(6,2) > 1
					    if (tem2 - Key_PLOT(6,2))>=1
						    tem1 = tem2 - Key_PLOT(6,2);
						end
					end
					for jjj =tem1:tem2
						% MD_Coor = load([Full_Pathname,'.mdco_',num2str(Itera_Num(jjj+1))]);
						% Old_MD_Coor = load([Full_Pathname,'.mdco_',num2str(Itera_Num(jjj))]);
						MD_Coor(1:num_molecule,1:3) = all_MD_Coor(jjj+1,1:num_molecule,1:3);
						Old_MD_Coor(1:num_molecule,1:3) = all_MD_Coor(jjj,1:num_molecule,1:3);
						all_c_x(1,1:num_molecule) = MD_Coor(1:num_molecule,2);
						all_c_x(2,1:num_molecule) = Old_MD_Coor(1:num_molecule,2);
						all_c_y(1,1:num_molecule) = MD_Coor(1:num_molecule,3);
						all_c_y(2,1:num_molecule) = Old_MD_Coor(1:num_molecule,3);
						
						% line(all_c_x,all_c_y,'LineWidth',1.0,'Color','blue');   %蓝色
						line(all_c_x,all_c_y,'LineWidth',1.0,'Color',[176/255,224/255,230/255]); %浅蓝色
					end
				end
			end
			disp(['      ----- Plotting positions of particles...'])
			% Read gauss point coordinates file.
			MD_Coor(1:num_molecule,1:3) = load([Full_Pathname,'.mdco_',num2str(Itera_Num(c_i))]);
			plot(MD_Coor(1:num_molecule,2),MD_Coor(1:num_molecule,3),'bo','MarkerSize',4,'Color','black','MarkerFaceColor','black')
			clear MD_Coor
		end
		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		%         绘制Box
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% if Key_PLOT(6,15)==1 
			% disp(['      ----- Plotting box......'])  
			% c_x = [box_zone_min_x,box_zone_max_x];
			% c_y = [box_zone_min_y,box_zone_min_y];
			% line(c_x,c_y,'LineWidth',1.5,'Color','green')
			% c_x = [box_zone_max_x,box_zone_max_x];
			% c_y = [box_zone_min_y,box_zone_max_y];
			% line(c_x,c_y,'LineWidth',1.5,'Color','green')
			% c_x = [box_zone_max_x,box_zone_min_x];
			% c_y = [box_zone_max_y,box_zone_max_y];
			% line(c_x,c_y,'LineWidth',1.5,'Color','green')
			% c_x = [box_zone_min_x,box_zone_min_x];
			% c_y = [box_zone_max_y,box_zone_min_y];
			% line(c_x,c_y,'LineWidth',1.5,'Color','green')
		% end
		if Key_PLOT(6,15)==1 
			disp(['      ----- Plotting box......'])  
			
			% 为了防止靠近边界的地方画不出来
			c_L = or_box_zone_max_x - or_box_zone_min_x;
			c_H = or_box_zone_max_y - or_box_zone_min_y;
			box_zone_min_x = or_box_zone_min_x + c_L*0.01;
			box_zone_max_x = or_box_zone_max_x - c_L*0.01;
			box_zone_min_y = or_box_zone_min_y + c_H*0.01;
			box_zone_max_y = or_box_zone_max_y - c_H*0.01;
			
			c_x = [box_zone_min_x,box_zone_max_x];
			c_y = [box_zone_min_y,box_zone_min_y];
			line(c_x,c_y,'LineWidth',1.5,'Color','black')
			c_x = [box_zone_max_x,box_zone_max_x];
			c_y = [box_zone_min_y,box_zone_max_y];
			line(c_x,c_y,'LineWidth',1.5,'Color','black')
			c_x = [box_zone_max_x,box_zone_min_x];
			c_y = [box_zone_max_y,box_zone_max_y];
			line(c_x,c_y,'LineWidth',1.5,'Color','black')
			c_x = [box_zone_min_x,box_zone_min_x];
			c_y = [box_zone_max_y,box_zone_min_y];
			line(c_x,c_y,'LineWidth',1.5,'Color','black')
			% 绘制5x5的标记线
			if Key_PLOT(6,9)==1 
				disp(['      ----- Plotting maker lines......'])  
				line_width_marker = 0.5;
				for i_marker_horizon =1:4
					delta_H = c_H/5.0;
					c_x = [box_zone_min_x,box_zone_max_x];
					c_y = [box_zone_min_y+delta_H*i_marker_horizon,box_zone_min_y+delta_H*i_marker_horizon];
					line(c_x,c_y,'LineWidth',line_width_marker,'Color',[192/255,192/255, 192/255])
				end
				for i_marker_vert =1:4
					delta_L = c_L/5.0;
					c_x = [box_zone_min_x+delta_L*i_marker_vert,box_zone_min_x+delta_L*i_marker_vert];
					c_y = [box_zone_min_y,box_zone_max_y];
					line(c_x,c_y,'LineWidth',line_width_marker,'Color',[192/255,192/255, 192/255])
				end		
			end	
		end
		
		
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Plot text on the left or the bottom of the figure.
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		plot_string_MatMES = ['PhiPsi  ',Version];
		plot_string_Frame  = ['Frame ',num2str(c_i),' / ',num2str(num_Iter)];
		plot_string_Time   = ['Time: ',num2str(c_i*delt_time_NewMark*1000),' ms'];
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
		range_W = abs(box_zone_max_x-box_zone_min_x);
		range_H = abs(box_zone_max_y-box_zone_min_y);
		if range_H >= 0.75*range_W      % Left
			loc_x = -range_H/2+ box_zone_min_x;
			loc_y =  box_zone_max_y-range_H*0.05;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_y =  box_zone_max_y-range_H*0.15;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_y =  box_zone_max_y-range_H*0.25;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_y =  box_zone_max_y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_y =  box_zone_max_y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
		else                            % Bottom
			loc_y =  box_zone_min_y-range_H*0.05;
			loc_x =  box_zone_min_x;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_x =  box_zone_min_x + range_W*0.25;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_x =  box_zone_min_x + range_W*0.45;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_x =  box_zone_min_x + range_W*0.60;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_x =  box_zone_min_x + range_W*0.60;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
		end
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		% Save the current figure
		%<<<<<<<<<<<<<<<<<<<<<<<<<<<
		deformation(i_output) = getframe(gcf);       
		im=frame2im(deformation(i_output));         
		[I,map]=rgb2ind(im,256);
		k=i_output-0;
		str1='_deformation';
		str2=Full_Pathname;
		FileName =[str2,str1,'.gif'];
		if k==1;
			imwrite(I,map,FileName,'gif','Loopcount',inf,'DelayTime',Time_Delay);    
		else
			imwrite(I,map,FileName,'gif','WriteMode','append','DelayTime',Time_Delay);
		end
		
		close
    end
end


% 输出完成信息
disp(' ')
disp('    Plot completed.')
disp(' ')

Cclock=clock;
% Display end time.
disp([' >> End time is ',num2str(Cclock(2)),'/',num2str(Cclock(3)),'/',num2str(Cclock(1))...
     ,' ',num2str(Cclock(4)),':',num2str(Cclock(5)),':',num2str(round(Cclock(6))),'.'])
	 
% Stop log file.
diary off;

