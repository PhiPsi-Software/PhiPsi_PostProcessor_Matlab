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

function [Offsetted_C_Up,Offsetted_C_Down] = Cal_Offseted_Single_Point(C,Line_AB,offset_delta)
% Offset the point c along the up and down normal of line_AB by the increment of offset_delta.
%
%              A                  B
%
%              ●-----------------●
%
%
%
%              ○
%
%              ●-----------------●
%              
%              ○

a_x = Line_AB(1,1);
a_y = Line_AB(1,2);
b_x = Line_AB(2,1);
b_y = Line_AB(2,2);
	

theta = atan2(b_y-a_y,b_x-a_x);

% Move C by offset_delta along the line_AB from B to A or from A to B.
% Case 1: C is at B point.
if (C(1)==b_x && C(2)==b_y)
	C(1)=C(1)-offset_delta*cos(theta);
	C(2)=C(2)-offset_delta*sin(theta);
% Case 2: C is at A point.
elseif  (C(1)==a_x && C(2)==a_y)
	C(1)=C(1)+offset_delta*cos(theta);
	C(2)=C(2)+offset_delta*sin(theta);
end

% Offset to left.
Offsetted_C_Up(1)  = C(1) - offset_delta*sin(theta);
Offsetted_C_Up(2)  = C(2) + offset_delta*cos(theta);
% Offset to right.
Offsetted_C_Down(1) = C(1) + offset_delta*sin(theta);
Offsetted_C_Down(2) = C(2) - offset_delta*cos(theta);