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

global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack Defor_Factor
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves num_Na_Crack
global Plot_Aperture_Curves Plot_Pressure_Curves Plot_Velocity_Curves Num_Step_to_Plot
global Key_TipEnrich Plot_Quantity_Curves Plot_Concentr_Curves Itera_Num
global Na_Crack_X Na_Crack_Y num_Na_Crack Itera_HF_Time Key_HF_Analysis
global num_Hole Hole_Coor   
global num_Circ_Inclusion num_Poly_Inclusion
global Analysis_type
global Key_Heaviside_Value 
global Key_Hole_Value 
global Key_Data_Format


disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    Attention :: Results number to plot: ',num2str(Num_Step_to_Plot)])
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    ']) 

Analysis_type = 0;
% Get the full name of files.
Full_Pathname = [Work_Dirctory,'\',Filename];

%获得分析类型号
if exist([Full_Pathname,'.post'], 'file') ==2 
    file_post = fopen([Full_Pathname,'.post']);
	lineNum = 0;
	while ~feof(file_post)
		lineNum = lineNum+1;
		TemData = fgetl(file_post);
		if lineNum==2  %第一行   
			Tem_Info = str2num(TemData);
		end
	end
	fclose(file_post); 
	Analysis_type       = Tem_Info(1);   %分析类型
    Key_TipEnrich       = Tem_Info(2);   %裂尖增强类型
	Key_Data_Format     = Tem_Info(3);   %保存的数据的类型
	Key_Heaviside_Value = Tem_Info(4);   %Value keyword of Heaviside enrichmenet function:-1 (1 and -1) or 0 (1 and 0)
    Key_Hole_Value      = Tem_Info(5);   %Value keyword of Hole enrichmenet function:-1 (1 and -1) or 0 (1 and 0)
end
% 如果是分子动力学模拟(*.post文件中存储的分析类型号=21),则很简单,直接调用
if Analysis_type==21
	Animate_MD
	%退出该子程序
	return
end
% 如果是近场动力学模拟(*.post文件中存储的分析类型号=31),则很简单,直接调用
if Analysis_type==31
	Animate_PD
	%退出该子程序
	return
end
% Read input geometry files.
Read_Geo
disp(['  '])    

% 如果Num_Step_to_Plot = -999,则程序自动寻找最后一步的稳定计算结果并后处理

if Num_Step_to_Plot == -999
    for i_Check =1:10000
	    if exist([Full_Pathname,'.disn_',num2str(i_Check)], 'file') ==2 
	        Num_Step_to_Plot = i_Check;
	    end
	end
	if Key_PLOT(5,1) ~=0 
		for i_Check =1:10000
			if exist([Full_Pathname,'.fdvl_',num2str(i_Check)], 'file') ==2 
				Num_Step_to_Plot = i_Check;
			end
		end
	end
end

%如果Num_Step_to_Plot=1,则终止程序 (2021-08-02)
if Num_Step_to_Plot==1
	disp(' >> Error :: only one step to plot, post 2 processor terminated.') 
	Error_Message
end

%如果Num_Step_to_Plot=-999,且没有计算结果,则终止程序
if Num_Step_to_Plot==-999
	disp(' >> Error :: No result files found, post-processor terminated.') 
	Error_Message
end

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
		Itera_Num(i_Fra) = 1;
		%所有破裂步间循环
		for i_ter = 1:num_Iter
			if ttt_DATA(i_ter,2)== i_Fra &&  ttt_DATA(i_ter,3) > Itera_Num(i_Fra)
				Itera_Num(i_Fra) = i_ter;
				Itera_HF_Time(i_Fra) = ttt_DATA(i_ter,4);
			end
		end
	end
% 如果不存在时间步文件
else
    for i = 1:Num_Step_to_Plot
        Itera_Num(i) = i;
	end
end

% Check if result files exist.
if exist([Full_Pathname,'.disn_',num2str(Itera_Num(Num_Step_to_Plot))], 'file') ~=2 && exist([Full_Pathname,'.fdvl_',num2str(Itera_Num(Num_Step_to_Plot))], 'file') ~=2
    disp(' >> Error :: No result files found, post-processor terminated.') 
	Error_Message
end

% Get the number of cracks if crack file exists.
for i =1:Num_Step_to_Plot
	if exist([Full_Pathname,'.crax_',num2str(Itera_Num(i))], 'file') ==2  
		file_Crack_X = fopen([Full_Pathname,'.crax_',num2str(num2str(Itera_Num(i)))]);
		row=0;
		while ~feof(file_Crack_X)
			row=row+sum(fread(file_Crack_X,10000,'*char')==char(10));
		end
		fclose(file_Crack_X);
		num_Crack(i) = row;
	else
		num_Crack(i) = 0;
	end
end

% Get the number of holes if hole file exists.
if exist([Full_Pathname,'.hlcr'], 'file') ==2 
	file_Hole = fopen([Full_Pathname,'.hlcr']);
	row=0;
	while ~feof(file_Hole)
		row=row+sum(fread(file_Hole,10000,'*char')==char(10));
	end
	fclose(file_Hole);
	num_Hole = row;
else
    num_Hole = 0;
end

% 获得圆形夹杂数目.
if exist([Full_Pathname,'.jzcr'], 'file') ==2  
	file_Circ_Inclusion = fopen([Full_Pathname,'.jzcr']);
	row=0;
	while ~feof(file_Circ_Inclusion)
		row=row+sum(fread(file_Circ_Inclusion,10000,'*char')==char(10));
	end
	fclose(file_Circ_Inclusion);
	num_Circ_Inclusion = row;
else
    num_Circ_Inclusion = 0;
end

% 获得多边形夹杂数目.
if exist([Full_Pathname,'.jzpx'], 'file') ==2  
	file_Poly_Inclusion = fopen([Full_Pathname,'.jzpx']);
	row=0;
	while ~feof(file_Poly_Inclusion)
		row=row+sum(fread(file_Poly_Inclusion,10000,'*char')==char(10));
	end
	fclose(file_Poly_Inclusion);
	num_Poly_Inclusion = row;
else
    num_Poly_Inclusion = 0;
end
% Get the number of natural cracks if natural crack file exists.
for i =1:Num_Step_to_Plot
	if exist([Full_Pathname,'.ncrx'], 'file') ==2  
		file_Na_Crack_X = fopen([Full_Pathname,'.ncrx']);
		row=0;
		while ~feof(file_Na_Crack_X)
			row=row+sum(fread(file_Na_Crack_X,10000,'*char')==char(10));
		end
		fclose(file_Na_Crack_X);
		num_Na_Crack(i) = row;
	else
		num_Na_Crack(i) = 0;
	end
end

% Read coordinates of natural cracks if natural cracks exist(读取天然裂缝的坐标).
if num_Na_Crack(1) ~= 0
	disp('    > Reading coordinates of natural cracks....') 
	file_Na_Crack_X = fopen([Full_Pathname,'.ncrx']);
	file_Na_Crack_Y = fopen([Full_Pathname,'.ncry']);
	Na_Crack_X = cell(num_Na_Crack(1));
	Na_Crack_Y = cell(num_Na_Crack(1));
	for i=1:num_Na_Crack(1)
		Na_Crack_X{i} = str2num(fgetl(file_Na_Crack_X));
        Na_Crack_Y{i} = str2num(fgetl(file_Na_Crack_Y));
	end
	fclose(file_Na_Crack_X);
	fclose(file_Na_Crack_Y);
end

% Plot.
Animate_Main(Num_Step_to_Plot) 

% 输出完成信息
disp(' ')
disp('    Plot completed.')
disp(' ')

Cclock=clock;
% Display end time.
disp([' >> End time is ',num2str(Cclock(2)),'/',num2str(Cclock(3)),'/',num2str(Cclock(1))...
     ,' ',num2str(Cclock(4)),':',num2str(Cclock(5)),':',num2str(round(Cclock(6))),'.'])
disp([' >> Total elapsed time is ',num2str(toc),' s, i.e. ',num2str(toc/60),' mins.'])
% Stop log file.
diary off;
