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
%function PhiPsi_CFD_Post_1_Go

global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack Defor_Factor
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves num_Na_Crack
global Plot_Aperture_Curves Plot_Pressure_Curves Plot_Velocity_Curves Num_Step_to_Plot
global Plot_Tan_Aper_Curves
global Key_TipEnrich Plot_Quantity_Curves Plot_Concentr_Curves     
global num_Hole Plot_Wpnp_Curves Plot_Wphp_Curves
global num_Circ_Inclusion num_Poly_Inclusion
global Key_Gas_Prod_rate Key_Gas_Production Key_One_Node_Pres
global Analysis_type Key_Data_Format
global Key_Heaviside_Value 
global Key_Hole_Value Plot_Inj_Pres_Curves
global num_Cross Arc_Crack_Coor file_Arc Yes_Arc_Crack
global num_HC

% Get the full name of files.
Full_Pathname = [Work_Dirctory,'\',Filename];

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

%读取网格数据
Read_Geo_CFD

% 如果Num_Step_to_Plot = -999,则程序自动寻找最后一步的稳定计算结果并后处理
if Num_Step_to_Plot == -999
    for i_Check =1:10000
	    if exist([Full_Pathname,'.pres_',num2str(i_Check)], 'file') ==2 
	        Num_Step_to_Plot = i_Check;
	    end
	end
end

%如果Num_Step_to_Plot=-999,则终止程序
if Num_Step_to_Plot==-999
	disp(' >> Error :: No result files found, post-processor terminated.') 
	Error_Message
end


	
% GUI交互信息.
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    Attention :: Results number to plot: ',num2str(Num_Step_to_Plot)])
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    ']) 



% Plot mesh.
if Key_PLOT(1,1) ==1
    Plot_CFD_Mesh(Num_Step_to_Plot)
end

% Plot contour.
if Key_PLOT(2,1) ~=0
    Plot_CFD_Contours(Num_Step_to_Plot)
end

% 输出完成信息.
disp(' ')
disp('    Plot completed.')
disp(' ')

Cclock=clock;

% Display end time.
disp([' >> End time is ',num2str(Cclock(2)),'/',num2str(Cclock(3)),'/',num2str(Cclock(1))...
     ,' ',num2str(Cclock(4)),':',num2str(Cclock(5)),':',num2str(round(Cclock(6))),'.'])
	 
disp([' >> Total elapsed time is ',num2str(toc),' s, i.e. ',num2str(toc/60),' mins.'])

% Close multiple processors.
% matlabpool close	 
	 
% Stop log file.
diary off;

