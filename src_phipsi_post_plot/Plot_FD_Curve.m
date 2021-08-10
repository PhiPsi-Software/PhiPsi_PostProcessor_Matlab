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

function Plot_FD_Curve

global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Title_Font Key_Figure_Control_Widget

disp('    > Plot F-D curves....') 
if exist([Full_Pathname,'.fdcu'], 'file') ==2
	namefile= [Full_Pathname,'.fdcu'];
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
	
	Max_FD_Point = max(ttt_DATA(1:num_Iter,2));
	
	disp(['    > 绘制载荷位移曲线, 节点号',num2str(ttt_DATA(1,1)),'...']) 
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	title('F-D curve','FontName',Title_Font,'FontSize',Size_Font)
	plot(ttt_DATA(1:Max_FD_Point,4)*1000,ttt_DATA(1:Max_FD_Point,3),'black-o','LineWidth',1)
	% set(gca,'xtick',1:1:Max_Frac)     
	xlabel('Displacement (mm)','FontName',Title_Font,'FontSize',Size_Font) 
	% xlabel('Time step','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('Force factor','FontName',Title_Font,'FontSize',Size_Font) 	
end

