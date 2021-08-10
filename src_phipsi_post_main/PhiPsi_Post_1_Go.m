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
%function PhiPsi_Post_1_Go

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
global num_HC num_Ellipse_Hole

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

% 如果是分子动力学模拟(*.post文件中存储的分析类型(Analysis_type=21)),则直接调用Plot_MD
% Analysis_type=21
if Analysis_type==21
    if Key_PLOT(6,1)~= 0
		Plot_MD
		return
	else
	    return
	end
end
% 如果是近场动力学模拟(*.post文件中存储的分析类型(Analysis_type=21)),则直接调用Plot_MD
% Analysis_type=31

if Analysis_type==31
    if Key_PLOT(8,1)~= 0
		Plot_PD
		return
	else
	    return
	end
end
% Read input geometry files.
Read_Geo

disp(['  '])    

% 如果Num_Step_to_Plot = -999,则程序自动寻找最后一步的稳定计算结果并后处理
if Num_Step_to_Plot == -999
    for i_Check =1:5000
	    if exist([Full_Pathname,'.disn_',num2str(i_Check)], 'file') ==2 
	        Num_Step_to_Plot = i_Check;
	    end
	end
    for i_Check =1:5000
	    if exist([Full_Pathname,'.fdvl_',num2str(i_Check)], 'file') ==2 
	        Num_Step_to_Plot = i_Check;
	    end
	end
end

%如果Num_Step_to_Plot=-999,则终止程序
if Num_Step_to_Plot==-999
	disp(' >> Error :: No result files found, post-processor terminated.') 
	Error_Message
end

%*************************
%曲线绘制相关设置的读取
%*************************
Plot_Pressure_Curves  = 0;                   %是否绘制裂缝内净水压分布曲线
Plot_Aperture_Curves  = 0;                   %是否绘制裂缝开度分布曲线
Plot_Tan_Aper_Curves  = 0;                   %是否绘制裂缝切向开度分布曲线
Plot_Velocity_Curves  = 0;                   %是否绘制裂缝计算点流速分布曲线
Plot_Quantity_Curves  = 0;                   %是否绘制裂缝计算点流量分布曲线
Plot_Concentr_Curves  = 0;                   %是否绘制裂缝计算点支撑剂浓度分布曲线
Plot_Wpnp_Curves      = 0;                   %是否绘制裂缝计算点支撑剂支撑裂缝开度(闭合压力为0),Width of Propped fracture when No Pressure
Num_Crack_HF_Curves   = [1];                 %绘制该裂纹对应的水力压裂分析相关曲线
Plot_SIF_KI_Curves    = 0;                   %是否绘制裂缝应力强度因子变化曲线
Num_Crack_SIF_Curves  = 1;                   %应力强度因子绘制裂缝号
Num_Tip_SIF_Curves    = 1;                   %应力强度因子绘制裂缝的裂尖号
Plot_Inj_Pres_Curves  = 0;                   %是否绘制裂缝注水点水压变化曲线
Key_InjP_Curves_x     = 2;                   %=1,以时间为x坐标轴;=2,以压裂步为x轴
Key_Gas_Prod_rate     = 0;                   %是否绘制页岩气产出率变化曲线
Key_Gas_Production    = 0;                   %是否绘制累积产量变化曲线
Key_One_Node_Pres     = 0;                   %是否绘制某一个点的压力变化


if Key_PLOT(7,1)==1
	%应力强度因子曲线的绘制
	Plot_SIF_KI_Curves    = Key_PLOT(7,2);    %是否绘制裂缝应力强度因子变化曲线
	Num_Crack_SIF_Curves  = Key_PLOT(7,3);    %应力强度因子绘制裂缝号码
	Num_Tip_SIF_Curves    = Key_PLOT(7,4);    %应力强度因子绘制裂缝的裂尖号
	%
	if Key_PLOT(7,5)==1
		Plot_Pressure_Curves  = 1;            %是否绘制裂缝内净水压分布曲线
		Plot_Aperture_Curves  = 0;            %是否绘制裂缝开度分布曲线
		Plot_Tan_Aper_Curves  = 0;            %是否绘制裂缝切向开度分布曲线
		Plot_Velocity_Curves  = 0;            %是否绘制裂缝计算点流速分布曲线
		Plot_Quantity_Curves  = 0;            %是否绘制裂缝计算点流量分布曲线
		Plot_Concentr_Curves  = 0;            %是否绘制裂缝计算点支撑剂浓度分布曲线
		Plot_Wpnp_Curves      = 0;            %是否绘制裂缝计算点支撑剂支撑裂缝开度(闭合压力�?),Width of Propped fracture without Pressure
	elseif Key_PLOT(7,5)==2
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 1;
		Plot_Tan_Aper_Curves  = 0;
		Plot_Velocity_Curves  = 0;
		Plot_Quantity_Curves  = 0;
		Plot_Concentr_Curves  = 0;
		Plot_Wpnp_Curves      = 0;
	elseif Key_PLOT(7,5)==3
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 0;
		Plot_Tan_Aper_Curves  = 1;
		Plot_Velocity_Curves  = 0;
		Plot_Quantity_Curves  = 0;
		Plot_Concentr_Curves  = 0;
		Plot_Wpnp_Curves      = 0;
	elseif Key_PLOT(7,5)==4
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 0;
		Plot_Tan_Aper_Curves  = 0;
		Plot_Velocity_Curves  = 1;
		Plot_Quantity_Curves  = 0;
		Plot_Concentr_Curves  = 0;
		Plot_Wpnp_Curves      = 0;
	elseif Key_PLOT(7,5)==5
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 0;
		Plot_Tan_Aper_Curves  = 0;
		Plot_Velocity_Curves  = 0;
		Plot_Quantity_Curves  = 1;
		Plot_Concentr_Curves  = 0;
		Plot_Wpnp_Curves      = 0;
	elseif Key_PLOT(7,5)==6
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 0;
		Plot_Tan_Aper_Curves  = 0;
		Plot_Velocity_Curves  = 0;
		Plot_Quantity_Curves  = 0;
		Plot_Concentr_Curves  = 1;
		Plot_Wpnp_Curves      = 0;
	elseif Key_PLOT(7,5)==7
		Plot_Pressure_Curves  = 0;
		Plot_Aperture_Curves  = 0;
		Plot_Tan_Aper_Curves  = 0;
		Plot_Velocity_Curves  = 0;
		Plot_Quantity_Curves  = 0;
		Plot_Concentr_Curves  = 0;
		Plot_Wpnp_Curves      = 1;
	end
		   
		
	Num_Crack_HF_Curves       = [Key_PLOT(7,6)];%绘制该裂纹对应的水力压裂分析相关曲线
	if Key_PLOT(7,7)==1
		Plot_Inj_Pres_Curves  = 1;              %是否绘制裂缝注水点水压变化曲线
		Key_InjP_Curves_x     = 1;              %=1,以时间为x坐标轴; =2,以压裂步为x坐标轴
	elseif Key_PLOT(7,7)==2
		Plot_Inj_Pres_Curves  = 1;              %是否绘制裂缝注水点水压变化曲线
		Key_InjP_Curves_x     = 2;              %=1,以时间为x坐标轴; =2,以压裂步为x坐标轴
	end
	
	%页岩气产量曲线.
	if Key_PLOT(7,8)==1
		Key_Gas_Prod_rate     = 1;              %是否绘制页岩气产出率变化曲线
	end
	if Key_PLOT(7,9)==1
		Key_Gas_Production    = 1;              %是否绘制累积产量变化曲线
	end
	if Key_PLOT(7,10)==1
		Key_One_Node_Pres     = 1;              %是否绘制某一个点的压力变化
	end
end
	
% GUI交互信息.
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    Attention :: Results number to plot: ',num2str(Num_Step_to_Plot)])
disp(['    ---#---#---#---#---#---#---#---#---#---#---']) 
disp(['    ']) 

% Check if result files exist.
if Key_PLOT(2,1) ~=0 || Key_PLOT(3,1) ~=0 || Key_PLOT(4,1) ~=0
	if exist([Full_Pathname,'.disn_',num2str(Num_Step_to_Plot)], 'file') ~=2 && ...    %弹性问题的节点位移
	   exist([Full_Pathname,'.fdvl_',num2str(Num_Step_to_Plot)], 'file') ~=2           %场问题的输出文件
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

% Get the Arc_Crack_Coor for arc crack if arcc file exists.
Yes_Arc_Crack = 0;
if exist([Full_Pathname,'.arcc_',num2str(Num_Step_to_Plot)], 'file') ==2  
	Yes_Arc_Crack = 1;
else
    % nothing
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

% Get the number of ellipse holes if ellipse hole file exists.
if exist([Full_Pathname,'.ehcr'], 'file') ==2  
	file_Ellipse_Hole = fopen([Full_Pathname,'.ehcr']);
	row=0;
	while ~feof(file_Ellipse_Hole)
		row=row+sum(fread(file_Ellipse_Hole,10000,'*char')==char(10));
	end
	fclose(file_Ellipse_Hole);
	num_Ellipse_Hole = row;
else
    num_Ellipse_Hole = 0;
end

% Get the number of crosses if hole file exists.
if exist([Full_Pathname,'.cscr'], 'file') ==2  
	file_Cross = fopen([Full_Pathname,'.cscr']);
	row=0;
	while ~feof(file_Cross)
		row=row+sum(fread(file_Cross,10000,'*char')==char(10));
	end
	fclose(file_Cross);
	num_Cross = row;
else
    num_Cross = 0;
end

% Get the number of HC if HC file exists.
if exist([Full_Pathname,'.hccr'], 'file') ==2  
	file_HC = fopen([Full_Pathname,'.hccr']);
	row=0;
	while ~feof(file_HC)
		row=row+sum(fread(file_HC,10000,'*char')==char(10));
	end
	fclose(file_HC);
	num_HC = row;
else
    num_HC = 0;
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
if exist([Full_Pathname,'.ncrx'], 'file') ==2  
	file_Na_Crack_X = fopen([Full_Pathname,'.ncrx']);
	row=0;
	while ~feof(file_Na_Crack_X)
		row=row+sum(fread(file_Na_Crack_X,10000,'*char')==char(10));
	end
	fclose(file_Na_Crack_X);
	num_Na_Crack = row;
else
    num_Na_Crack = 0;
end

% Plot.
Plot_Main(Num_Step_to_Plot) 

% 绘制HF曲线.
if Plot_Aperture_Curves == 1 | Plot_Pressure_Curves==1| Plot_Tan_Aper_Curves == 1 | ...
   Plot_Velocity_Curves==1   | Plot_Quantity_Curves==1| Plot_Concentr_Curves ==1 | ...
   Plot_Wpnp_Curves == 1     | Plot_Wphp_Curves ==1
    Plot_HF_curves(Num_Step_to_Plot) 
end

% 绘制应力强度因子曲线.
if Plot_SIF_KI_Curves==1
    Plot_SIF_curves(Num_Step_to_Plot,Num_Crack_SIF_Curves,Num_Tip_SIF_Curves) 
end

% 绘制页岩气产量曲线.
if Key_Gas_Prod_rate==1 | Key_Gas_Production==1
    Plot_Gas_Production_curves(Num_Step_to_Plot)
end

% 绘制某一个点的压力变化曲线.
if Key_One_Node_Pres==1 
    Plot_one_node_Pressure_curves(Num_Step_to_Plot)
end

if Plot_Inj_Pres_Curves==1
	if exist([Full_Pathname,'.injp'], 'file') ==2  
		Plot_Injection_Pressure_curves(Num_Step_to_Plot,Key_InjP_Curves_x)
    end
end

%绘制载荷位移曲线.
if Key_PLOT(7,1)==1
	if Key_PLOT(7,11)==1
	    Plot_FD_Curve
	end
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

