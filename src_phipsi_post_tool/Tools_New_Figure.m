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

global Key_Figure_off
% New figure.

if Key_Figure_off == 0     % Visible on.
     % c_figure = figure('units','normalized','Visible','on','Renderer','OpenGL');   %采用OpenGL渲染器   
	 c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on','Renderer','OpenGL');   %采用OpenGL渲染器   
	 % c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on','Renderer','painters'); %默认的渲染器,最慢
	 % c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on','Renderer','zbuffer');    %zbuffer渲染器
	 % c_figure = figure('units','normalized','position',[0.38,0.38,0.42,0.42],'Visible','on','Renderer','zbuffer');    %zbuffer渲染器
	 
	 opengl hardware     %打开opengl硬件加速
	 % opengl software     %打开opengl软件加速
	 set(gcf,'renderer','opengl')
elseif Key_Figure_off == 1 % Visible off.
     c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','off');
	 set(gcf,'renderer','opengl')
end

