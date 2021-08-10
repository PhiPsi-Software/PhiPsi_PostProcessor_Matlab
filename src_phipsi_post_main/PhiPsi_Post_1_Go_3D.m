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

% Start and define global variables
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot
global Key_TipEnrich
global Key_Data_Format Only_Plot_Mesh

% Get the full name of files.
Full_Pathname = [Work_Dirctory,'\',Filename];

% Only Plot Mesh (2021-08-07).
Only_Plot_Mesh =0;
if Key_PLOT(1,1)==1 & sum(Key_PLOT(2:5,1))==0 
    Only_Plot_Mesh =1;
end

%获得分析类型等数据
if exist([Full_Pathname,'.post'], 'file') ==2 
    file_post = fopen([Full_Pathname,'.post']);
	lineNum = 0;
	while ~feof(file_post)
		lineNum = lineNum+1;
		TemData = fgetl(file_post);
		if lineNum==2 
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

% 如果Num_Step_to_Plot = -999,则程序自动寻找最后一步的稳定计算结果并后处理
if Num_Step_to_Plot == -999
    for i_Check =1:1000
	    if exist([Full_Pathname,'.disn_',num2str(i_Check)], 'file') ==2 
	        Num_Step_to_Plot = i_Check;
	    end
	end
end

%如果Num_Step_to_Plot=-999,且没有计算结果,则终止程序
if Num_Step_to_Plot==-999
    if Only_Plot_Mesh ==0
		disp(' >> Error :: No result files found, post-processor terminated.') 
		Error_Message
	elseif Only_Plot_Mesh ==1
		disp(' >> WARNNING :: Plot only mesh.') 
		Num_Step_to_Plot = 1;
	end
end

disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    Attention :: Results number to plot: ',num2str(Num_Step_to_Plot)])
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    ']) 

% Get the full name of files.
Full_Pathname = [Work_Dirctory,'\',Filename];

% Read input geometry files.
Read_Geo3D
disp(['  '])    

% Check if result files exist.
if exist([Full_Pathname,'.disn_',num2str(Num_Step_to_Plot)], 'file') ~=2
    if Only_Plot_Mesh ==0
		disp(' >> Error :: No result files found, post-processor terminated.') 
		Error_Message
	end	
end

% Get the number of cracks if crack file exists.
if exist([Full_Pathname,'.crax_',num2str(Num_Step_to_Plot)], 'file') ==2  
	file_Crack_X = fopen([Full_Pathname,'.crax_',num2str(Num_Step_to_Plot)]);
	row=0;
	while ~feof(file_Crack_X)
		row=row+sum(fread(file_Crack_X,10000,'*char')==char(10));
	end
	fclose(file_Crack_X);
	num_Crack(Num_Step_to_Plot) = row;
else
    num_Crack(Num_Step_to_Plot) = 0;
end

% Plot.
Plot_Main3D(Num_Step_to_Plot) 

% 绘制HF曲线
% if Plot_Aperture_Curves == 1 || Plot_Pressure_Curves==1
    % Plot_curves(Num_Step_to_Plot) 
% end

Cclock=clock;
% Display end time.
disp([' >> End time is ',num2str(Cclock(2)),'/',num2str(Cclock(3)),'/',num2str(Cclock(1))...
     ,' ',num2str(Cclock(4)),':',num2str(Cclock(5)),':',num2str(round(Cclock(6))),'.'])
	 
disp([' >> Total elapsed time is ',num2str(toc),' s, i.e. ',num2str(toc/60),' mins.'])
% Stop log file.
diary off;

