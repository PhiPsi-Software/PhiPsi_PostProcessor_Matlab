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

function Plot_PD
% 近场动力学模拟相关的后处理.

global Key_PLOT
global Size_Font 
global Full_Pathname Num_Step_to_Plot
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global box_zone_min_x box_zone_max_x box_zone_min_y box_zone_max_y
global Title_Font Key_Figure_Control_Widget

disp(['      ----- Read files......'])

% 如果存在时间步文件
if  exist([Full_Pathname,'.hftm'], 'file') ==2 
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
	%最大计算步数
	Max_Steps = max(ttt_DATA(1:num_Iter,2));
	%提取每个破裂步对应的迭代步号
	for i_Fra = 1:Max_Steps
		Itera_Num(i_Fra) = 0;
		%所有计算步间循环
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

if Num_Step_to_Plot == -999
    Num_Step_to_Plot = Itera_Num(num_Iter);
end

%检查要绘制的是否存在
if exist([Full_Pathname,'.pdrs_',num2str(Num_Step_to_Plot)], 'file') ~=2  
	disp(' >> Error :: No result files found, post-processor terminated.') 
	Error_Message
end
dscale = Key_PLOT(8,6);

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
    % all_PD_RS=cell(num_Iter,)
    for ii=1:num_Iter
	    all_PD_RS(ii,1:num_points,1:6) = load([Full_Pathname,'.pdrs_',num2str(Itera_Num(ii))]);
    end
	%获取当前步对应的c_i号
	for uuu =1:num_Iter
	    if Itera_Num(uuu)==Num_Step_to_Plot
		    c_i = uuu;
		end
	end
end

% New figure.
Tools_New_Figure
hold on;
title('Position of points','FontName',Title_Font,'FontSize',Size_Font)
axis off; axis equal;
axis([zone_min_x zone_max_x zone_min_y zone_max_y]);

% 绘制轨迹线.
if Key_PLOT(8,2) >= 1
    disp(['      ----- Plotting moving path......'])  
	if c_i>=2
	    tem1 = 1;
	    tem2 = c_i-1;
		%%%%%如果Key_PLOT(8,2)=2,则只绘制最近 Key_PLOT(8,2)个计算步的轨迹
		if (tem2-tem1) > Key_PLOT(8,2)
		    tem1 = tem2 - Key_PLOT(8,2);
		end
		for jjj =tem1:tem2
			%%%%% PD_Results = load([Full_Pathname,'.mdco_',num2str(Itera_Num(jjj+1))]);
			%%%%% Old_PD_Results = load([Full_Pathname,'.mdco_',num2str(Itera_Num(jjj))]);
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


% 绘制点的变形后位置.
if Key_PLOT(8,1) == 1
    if exist([Full_Pathname,'.pdrs_',num2str(Num_Step_to_Plot)], 'file') ==2  
		disp(['      ----- Plotting new positions of points...'])
		PD_Results = load([Full_Pathname,'.pdrs_',num2str(Num_Step_to_Plot)]);
		plot(PD_Results(:,2)+dscale*PD_Results(:,4),PD_Results(:,3)+dscale*PD_Results(:,5),'ro','MarkerSize',Key_PLOT(8,7),'Color','black','MarkerFaceColor','black')
		% 绘制点的原始位置.
		if Key_PLOT(8,8) == 1
		    disp(['      ----- Plotting old positions of points...'])
			plot(PD_Results(:,2),PD_Results(:,3),'ro','MarkerSize',Key_PLOT(8,7),'Color','green','MarkerFaceColor','black')
		end
		clear PD_Results
	end
end

% 绘制Box
if Key_PLOT(8,15)==1 
	%如果定义了破裂区
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
end


% Save pictures.
Save_Picture(c_figure,Full_Pathname,'pdrs')

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
