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

function Plot_SIF_curves(POST_Substep,num_Cr_to_Plot,num_Tip)
% ����Ӧ��ǿ����������.
global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Title_Font Key_Figure_Control_Widget

%********************************
%��ȡ�����Ѳ���Ӧ���ܵĵ�������
%********************************
disp('    > ��ȡ�����Ѳ���Ӧ�ĵ�������(hftm�ļ�)....') 

if exist([Full_Pathname,'.hftm'], 'file') ==2  
	namefile= [Full_Pathname,'.hftm'];
	data=fopen(namefile,'r'); 
	lineNum = 0;
	num_Iter = 0;
	while ~feof(data)
		lineNum = lineNum+1;
		TemData = fgetl(data);    
		if lineNum>=2   %��һ�����ļ���ʶ��,�����ȡ
			num_Iter = num_Iter+1;                     %�ܵĵ�������
			c_num   = size(str2num(TemData),2); 	   
			ttt_DATA(num_Iter,1:4)  = str2num(TemData);
		end
	end
	fclose(data); 
else
    %���ļ�������,��ֱ���˳�
    return
end

%������Ѳ���
Max_Frac = max(ttt_DATA(1:num_Iter,2));
%��ȡÿ�����Ѳ���Ӧ�ĵ�������
for i_Fra = 1:Max_Frac
    Itera_Num(i_Fra) = 1;
	%�������Ѳ���ѭ��
	for i_ter = 1:num_Iter
	    if ttt_DATA(i_ter,2)==i_Fra &&  ttt_DATA(i_ter,3) >Itera_Num(i_Fra)
		    Itera_Num(i_Fra) = i_ter;
		end
	end
end

%******************************************
%��ȡÿ�����Ѳ���Ӧ���ѷ�ŵ�Ӧ��ǿ������
%******************************************
for i_Fra = 1:Max_Frac
    c_Iter_num = Itera_Num(i_Fra);
	if exist([Full_Pathname,'.sifs_',num2str(c_Iter_num)], 'file') ==2 
	    namefile= [Full_Pathname,'.sifs_',num2str(c_Iter_num)];
	    data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);     %ÿ�й�4������: �Ѽ�1��KI���Ѽ�1��KII���Ѽ�2��KI���Ѽ�2��KII      
			if lineNum ==num_Cr_to_Plot
				Plot_SIFs(i_Fra,1:4)  = str2num(TemData);
			end
		end
		fclose(data); 
	else
	    Plot_SIFs(i_Fra,1:4)  = [0,0,0,0];
	end
end


%************************************
% ���Ƹ������Ƶ�Ӧ��ǿ����������
%************************************
disp(['    > �������� ',num2str(num_Cr_to_Plot),'�Ѽ� ',num2str(num_Tip),' I��Ӧ��ǿ����������...']) 
c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
hold on;
if num_Cr_to_Plot==1 
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==2
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==3
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==4
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==5
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==6
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==7
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==8
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==9
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
elseif num_Cr_to_Plot==10
	title('KI of crack 1 (MPa*m^1^/^2)','FontName',Title_Font,'FontSize',Size_Font)
end
Plot_x = 1:1:Max_Frac;
if num_Tip==1
	plot(Plot_x ,Plot_SIFs(1:Max_Frac,1)/1.0E6,'black-o','LineWidth',1)
	%��Ļ���Ӧ��ǿ������
	% Plot_SIFs(1:Max_Frac,1)/1.0E6
elseif num_Tip==2
	plot(Plot_x ,Plot_SIFs(1:Max_Frac,3)/1.0E6,'black-o','LineWidth',1)
	
end

%�����ļ�
% Plot_x=Plot_x/66;
% K = Plot_SIFs(1:Max_Frac,3);
% K=K/10e6/sqrt(2);
% Plot_x_T = Plot_x';
% K_T = K';
% save X:\dynamic_time_phipsi.txt Plot_x_T -ascii
% save X:\dynamic_K_phipsi.txt K -ascii

% set(gca,'xtick',1:1:Max_Frac)     
xlabel('Fracturing step','FontName',Title_Font,'FontSize',Size_Font) 
% xlabel('Time step','FontName',Title_Font,'FontSize',Size_Font) 
ylabel('KI','FontName',Title_Font,'FontSize',Size_Font) 	