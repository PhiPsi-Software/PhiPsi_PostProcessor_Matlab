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
% 绘制应力强度因子曲线.
global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Title_Font Key_Figure_Control_Widget

%********************************
%读取各破裂步对应的总的迭代次数
%********************************
disp('    > 读取各破裂步对应的迭代步号(hftm文件)....') 

if exist([Full_Pathname,'.hftm'], 'file') ==2  
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
else
    %若文件不存在,则直接退出
    return
end

%最大破裂步数
Max_Frac = max(ttt_DATA(1:num_Iter,2));
%提取每个破裂步对应的迭代步号
for i_Fra = 1:Max_Frac
    Itera_Num(i_Fra) = 1;
	%所有破裂步间循环
	for i_ter = 1:num_Iter
	    if ttt_DATA(i_ter,2)==i_Fra &&  ttt_DATA(i_ter,3) >Itera_Num(i_Fra)
		    Itera_Num(i_Fra) = i_ter;
		end
	end
end

%******************************************
%读取每个破裂步对应的裂缝号的应力强度因子
%******************************************
for i_Fra = 1:Max_Frac
    c_Iter_num = Itera_Num(i_Fra);
	if exist([Full_Pathname,'.sifs_',num2str(c_Iter_num)], 'file') ==2 
	    namefile= [Full_Pathname,'.sifs_',num2str(c_Iter_num)];
	    data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);     %每行共4个变量: 裂尖1的KI，裂尖1的KII，裂尖2的KI，裂尖2的KII      
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
% 绘制各条裂纹的应力强度因子曲线
%************************************
disp(['    > 绘制裂纹 ',num2str(num_Cr_to_Plot),'裂尖 ',num2str(num_Tip),' I型应力强度因子曲线...']) 
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
	%屏幕输出应力强度因子
	% Plot_SIFs(1:Max_Frac,1)/1.0E6
elseif num_Tip==2
	plot(Plot_x ,Plot_SIFs(1:Max_Frac,3)/1.0E6,'black-o','LineWidth',1)
	
end

%保存文件
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