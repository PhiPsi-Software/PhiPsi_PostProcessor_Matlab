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

function Plot_Gas_Production_curves(POST_Substep)

global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Key_Gas_Prod_rate Key_Gas_Production
global Title_Font Key_Figure_Control_Widget

%********************************
%读取各破裂步对应的总的迭代次数
%********************************
disp('    > 读取gasp文件....') 

if exist([Full_Pathname,'.gasp'], 'file') ==2  
	namefile= [Full_Pathname,'.gasp'];
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
else
    %若文件不存在,则直接退出
    return
end

%*************
%绘制曲线
%*************
if Key_Gas_Prod_rate==1
	disp(['    > 绘制产量变化曲线...']) 
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
    title('Production rates','FontName',Title_Font,'FontSize',Size_Font)
    plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,3),'black-o','LineWidth',1)
	% set(gca,'xtick',1:1:Max_Frac)     
    xlabel('Time (day)','FontName',Title_Font,'FontSize',Size_Font) 
	% xlabel('Time step','FontName',Title_Font,'FontSize',Size_Font) 
    ylabel('Gas production rate (MMscf/day)','FontName',Title_Font,'FontSize',Size_Font) 	
end
if Key_Gas_Production==1
	disp(['    > 绘制累积产量曲线...']) 
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
    title('Cumulative gas production','FontName',Title_Font,'FontSize',Size_Font)
    plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,4),'black-o','LineWidth',1)
	% set(gca,'xtick',1:1:Max_Frac)     
    xlabel('Time (day)','FontName',Title_Font,'FontSize',Size_Font) 
	% xlabel('Time step','FontName',Title_Font,'FontSize',Size_Font) 
    ylabel('Cumulative gas production (MMscf)','FontName',Title_Font,'FontSize',Size_Font) 	
end

