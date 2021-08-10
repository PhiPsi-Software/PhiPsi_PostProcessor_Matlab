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

function [PS1,PS2,theta] = Cal_Principal_Stress(XX,YY,XY,Key_Matrix)
% This function calculate the principal stress of each Gauss point of entire element.

if Key_Matrix ==0               % Only calculate one point.
	PS1 = (XX+YY)/2+sqrt(((XX-YY)/2)^2+XY^2);
	PS2 = (XX+YY)/2-sqrt(((XX-YY)/2)^2+XY^2);
	if XX >= YY
		theta = atan(2*XY/(XX-YY))/2+pi/2;
	elseif XX < YY
		theta = atan(2*XY/(XX-YY))/2;
	end
elseif Key_Matrix ==1           % Calculate matrix.
	PS1   = (XX+YY)/2+(((XX-YY)/2).^2+XY.^2).^0.5;
	PS2   = (XX+YY)/2-(((XX-YY)/2).^2+XY.^2).^0.5;
	Flag  = XX>=YY;
	theta = atan(2*XY./(XX-YY))/2 + pi/2*Flag;
end
