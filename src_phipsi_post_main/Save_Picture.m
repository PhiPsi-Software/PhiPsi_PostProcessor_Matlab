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

function Save_Picture(c_figure,Full_Pathname,type)
% This function save pictures.

global Key_Figure_off

% If visible of figure is on, then clear picture can be saved by the following codes:
if Key_Figure_off==0
	Picture = getframe(gcf);       
	im=frame2im(Picture);         
	[I,map]=rgb2ind(im,256);
	str1=['_',type];
	str2=Full_Pathname;
	FileName1 =[str2,str1,'.png'];
	imwrite(I,map,FileName1,'png'); 
% If visible of figure is off, then not very clear picture can be saved by the "saveas" function:
elseif Key_Figure_off==1
	str1=['_',type];
	str2=Full_Pathname;
	FileName1 =[str2,str1,'.png'];
	saveas(c_figure,FileName1);
end