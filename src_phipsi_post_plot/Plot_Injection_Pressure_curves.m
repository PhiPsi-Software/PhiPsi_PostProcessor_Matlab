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
function Plot_Injection_Pressure_curves(POST_Substep,Key_InjP_Curves_x)
% ����ˮ��ѹ�ѷ���עˮ��ѹ���仯����.

global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Title_Font Key_Figure_Control_Widget

%***********
%��ȡ����
%***********
disp('    > Read *.injp file....') 
namefile= [Full_Pathname,'.injp'];
data=fopen(namefile,'r'); 
lineNum = 0;
num_frac = 0;
while ~feof(data)
	lineNum = lineNum+1;
	TemData = fgetl(data);    
	if lineNum>=2   %��һ�����ļ���ʶ��,�����ȡ
	    num_frac = num_frac+1;                     %�ܵ����Ѳ���	   
	    ttt_DATA(num_frac,1:3)  = str2num(TemData);
	end
end
fclose(data); 

%*************************************
% ���Ƹ������Ƶ�����(ˮ��ѹ�����)
%*************************************
disp(['    > ����עˮ��ˮѹ�仯����...']) 
c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
hold on;
title('Injection pressure (MPa)','FontName',Title_Font,'FontSize',Size_Font)
%��עˮʱ��Ϊx��
if Key_InjP_Curves_x==1
    plot(ttt_DATA(1:num_frac,2)/60.0,ttt_DATA(1:num_frac,3)/1.0e6,'black-o','LineWidth',1)
    xlabel('Injection time (min)','FontName',Title_Font,'FontSize',Size_Font) 
%��ѹ�Ѳ���Ϊx��
elseif Key_InjP_Curves_x==2
    plot(ttt_DATA(1:num_frac,1),ttt_DATA(1:num_frac,3)/1.0e6,'black-o','LineWidth',1)
    xlabel('Fracturing step','FontName',Title_Font,'FontSize',Size_Font)     
	set(gca,'xtick',1:1:num_frac)     
end
ylabel('Injection pressure (MPa)','FontName',Title_Font,'FontSize',Size_Font) 