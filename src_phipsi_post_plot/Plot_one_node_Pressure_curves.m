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

function Plot_one_node_Pressure_curves(POST_Substep)
% ����ĳһ�����ѹ���仯����
global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Key_Gas_Prod_rate Key_Gas_Production Key_One_Node_Pres
global Title_Font Key_Figure_Control_Widget

%��ȡ�����Ѳ���Ӧ���ܵĵ�������
disp('    > ��ȡonpr�ļ�....') 

if exist([Full_Pathname,'.onpr'], 'file') ==2  
	namefile= [Full_Pathname,'.onpr'];
	data=fopen(namefile,'r'); 
	lineNum = 0;
	num_Iter = 0;
	while ~feof(data)
		lineNum = lineNum+1;
		TemData = fgetl(data);    
		if lineNum>=2   %��һ�����ļ���ʶ��,�����ȡ
			num_Iter = num_Iter+1;                     %�ܵĵ�������
			c_num   = size(str2num(TemData),2); 	   
			ttt_DATA(num_Iter,1:12)  = str2num(TemData);
		end
	end
	fclose(data); 
else
    %���ļ�������,��ֱ���˳�
    return
end

%��������
if Key_One_Node_Pres==1
	disp(['    > ���Ƶ��ѹ���仯����...']) 
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
    title('Node pressure','FontName',Title_Font,'FontSize',Size_Font)
	if ttt_DATA(1,3) ~=0
		h1=plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,4),'black-o','LineWidth',1);
		str1 = ['Node number: ',num2str(ttt_DATA(1,3))];
		legend(h1,str1);
	end
	if ttt_DATA(1,5) ~=0
		h2=plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,6),'red-x','LineWidth',2);
		str2 = ['Node number: ',num2str(ttt_DATA(1,5))];
		legend([h1,h2],str1,str2);
	end
	if ttt_DATA(1,7) ~=0
		h3=plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,8),'blue-*','LineWidth',1);
		str3 = ['Node number: ',num2str(ttt_DATA(1,7))];
		legend([h1,h2,h3],str1,str2,str3);
	end
	if ttt_DATA(1,9) ~=0
		h4=plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,10),'green-.','LineWidth',2);
		str4 = ['Node number: ',num2str(ttt_DATA(1,9))];
		legend([h1,h2,h3,h4],str1,str2,str3,str4);
	end
	if ttt_DATA(1,11) ~=0
		h5=plot(ttt_DATA(1:num_Iter,2),ttt_DATA(1:num_Iter,12),'c-d','LineWidth',1);
		str5 = ['Node number: ',num2str(ttt_DATA(1,11))];
		legend([h1,h2,h3,h4,h5],str1,str2,str3,str4,str5);
	end	
	% set(gca,'xtick',1:1:Max_Frac)     
    xlabel('Time (months)','FontName',Title_Font,'FontSize',Size_Font) 
	% xlabel('Time step','FontName',Title_Font,'FontSize',Size_Font) 
    ylabel('Node pressure(MPa)','FontName',Title_Font,'FontSize',Size_Font) 	
end
