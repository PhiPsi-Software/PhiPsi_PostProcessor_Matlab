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

function Animate_PD
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

disp('    Generating PD animations......')


dscale = Key_PLOT(8,6);


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

%获取点的数目
PD_Results = load([Full_Pathname,'.pdrs_',num2str(Num_Step_to_Plot)]);
num_points = size(PD_Results,1);
disp(['    Number of points:',num2str(num_points)]);

%获取变形后的模型坐标范围
zone_min_x = min(PD_Results(1:num_points,2) + PD_Results(1:num_points,4)*dscale);
zone_max_x = max(PD_Results(1:num_points,2) + PD_Results(1:num_points,4)*dscale);
zone_min_y = min(PD_Results(1:num_points,3) + PD_Results(1:num_points,5)*dscale);
zone_max_y = max(PD_Results(1:num_points,3) + PD_Results(1:num_points,5)*dscale);

% 如果需要绘制轨迹线,则很耗时间,可以先把所有位置文件读入内存
if Key_PLOT(8,2) >= 1
    for ii=1:num_Iter
	    all_PD_RS(ii,1:num_points,1:6) = load([Full_Pathname,'.pdrs_',num2str(Itera_Num(ii))]);
    end
end
		
% Loop through each frame to plot and store.
i_output=0;
for c_i=1:num_Iter
	if mod(c_i,Output_Freq)==0
	    i_output = i_output + 1;

		scale = Key_PLOT(8,6);
		
		% New figure.
		Tools_New_Figure
		hold on;
		title('Position of points','FontName',Title_Font,'FontSize',Size_Font)
		axis off; 
		axis equal;
        axis([zone_min_x zone_max_x zone_min_y zone_max_y]);

		%<<<<<<<<<<<<<<<<<<<<<<<<
		% 绘制点的位置及轨迹线.
		%<<<<<<<<<<<<<<<<<<<<<<<<
		if Key_PLOT(8,1) == 1
			% 轨迹线
			if Key_PLOT(8,2) >= 1
			    if c_i>=2
					disp(['      ----- Plotting moving paths of points...'])
					tem1 = 1;
					tem2 = c_i-1;
					%如果Key_PLOT(6,2)=2,则只绘制最近 Key_PLOT(6,2)个计算步的轨迹
					if Key_PLOT(6,2) > 1
					    if (tem2 - Key_PLOT(8,2))>=1
						    tem1 = tem2 - Key_PLOT(8,2);
						end
					end
					for jjj =tem1:tem2
						PD_Results(1:num_points,1:6) = all_PD_RS(jjj+1,1:num_points,1:6);
						Old_PD_Results(1:num_points,1:6) = all_PD_RS(jjj,1:num_points,1:6);
						all_c_x(1,1:num_points) = PD_Results(1:num_points,2)+PD_Results(1:num_points,4)*dscale;
						all_c_x(2,1:num_points) = Old_PD_Results(1:num_points,2)+Old_PD_Results(1:num_points,4)*dscale;
						all_c_y(1,1:num_points) = PD_Results(1:num_points,3)+PD_Results(1:num_points,5)*dscale;
						all_c_y(2,1:num_points) = Old_PD_Results(1:num_points,3)+Old_PD_Results(1:num_points,5)*dscale;
						line(all_c_x,all_c_y,'LineWidth',1.0,'Color','blue');
					end
				end
			end
			disp(['      ----- Plotting positions of points...'])
			% Read gauss point coordinates file.
			PD_rs(1:num_points,1:6) = load([Full_Pathname,'.pdrs_',num2str(Itera_Num(c_i))]);
			plot(PD_rs(1:num_points,2)+dscale*PD_rs(1:num_points,4),PD_rs(1:num_points,3)+dscale*PD_rs(1:num_points,5),'bo','MarkerSize',4,'Color','black','MarkerFaceColor','black')
			clear PD_Results
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
		range_W = abs(zone_max_x-zone_min_x);
		range_H = abs(zone_max_y-zone_min_y);
		if range_H >= 0.75*range_W      % Left
			loc_x = -range_H/2+ zone_min_x;
			loc_y =  zone_max_y-range_H*0.05;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_y =  zone_max_y-range_H*0.15;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_y =  zone_max_y-range_H*0.25;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_y =  zone_max_y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_y =  zone_max_y-range_H*0.35;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
		else                            % Bottom
			loc_y =  zone_min_y-range_H*0.05;
			loc_x =  zone_min_x;
			text(loc_x,loc_y,plot_string_MatMES,'color','black');
			loc_x =  zone_min_x + range_W*0.25;
			text(loc_x,loc_y,plot_string_Scale,'color','black');
			loc_x =  zone_min_x + range_W*0.45;
			text(loc_x,loc_y,plot_string_Frame,'color','black');
			if Key_Dynamic == 1
				loc_x =  zone_min_x + range_W*0.60;
				text(loc_x,loc_y,plot_string_Time, 'color','black');
			end
			if  exist([Full_Pathname,'.hftm'], 'file') ==2 
				loc_x =  zone_min_x + range_W*0.60;
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

