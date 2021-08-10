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

function [Num] = Cal_Ele_Num_by_Coors(X,Y)
% This function calculate the element num by the coordinates of the point inside the element.

global G_NN G_X_NODES G_Y_NODES G_X_Min G_X_Max G_Y_Min G_Y_Max

Num       = 0;
tem     = intersect(intersect(find(G_X_Min<=X),find(G_X_Max>=X)), ...
                    intersect(find(G_Y_Min<=Y),find(G_Y_Max>=Y)));
	
for i=1:size(tem,2)
	[Yes_In Yes_On]= inpolygon(X,Y,G_X_NODES(:,tem(i)),G_Y_NODES(:,tem(i)));  
    if Yes_In ==1 || Yes_On==1
		Num = tem(i);
		break
	end
end	
