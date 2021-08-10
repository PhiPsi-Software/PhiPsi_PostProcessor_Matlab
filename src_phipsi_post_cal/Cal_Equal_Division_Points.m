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

function [Div_Points,Offsetted_D_P_Up,Offsetted_D_P_Down] = ...
        Cal_Equal_Division_Points(Num_Diversion,Line_AB,offset_delta,Key_Include_Endpoint)
% Get the equal diversion points of line AB, then offset them by a small delta.
% Diversion point are arranged from A to B.

a_x = Line_AB(1,1);
a_y = Line_AB(1,2);
b_x = Line_AB(2,1);
b_y = Line_AB(2,2);

for i = 1:Num_Diversion-1	
    Div_Points(i,1) = (i*b_x+(Num_Diversion-i)*a_x)/Num_Diversion;
	Div_Points(i,2) = (i*b_y+(Num_Diversion-i)*a_y)/Num_Diversion;
end

% If Key_Include_Endpoint==1 then add the Endpoints of Line_AB. 
if Key_Include_Endpoint==1
	ttt = Div_Points;
	num_ttt = size(ttt,1);
	Div_Points(1,:)= [a_x a_y];
	for i=1:size(Div_Points,1)
		Div_Points(i+1,:)   = ttt(i,:);
	end
	Div_Points(num_ttt+2,:) = [b_x b_y];
end

theta = atan2(b_y-a_y,b_x-a_x);

% Offset endpoints.
Div_Points(1,1) = Div_Points(1,1)+offset_delta*cos(theta);
Div_Points(1,2) = Div_Points(1,2)+offset_delta*sin(theta);

Div_Points(size(Div_Points,1),1) = Div_Points(size(Div_Points,1),1)-offset_delta*cos(theta);
Div_Points(size(Div_Points,1),2) = Div_Points(size(Div_Points,1),2)-offset_delta*sin(theta);

% Offset to left.
Offsetted_D_P_Up(:,1) = Div_Points(:,1) - offset_delta*sin(theta);
Offsetted_D_P_Up(:,2) = Div_Points(:,2) + offset_delta*cos(theta);
% Offset to right.
Offsetted_D_P_Down(:,1) = Div_Points(:,1) + offset_delta*sin(theta);
Offsetted_D_P_Down(:,2) = Div_Points(:,2) - offset_delta*cos(theta);