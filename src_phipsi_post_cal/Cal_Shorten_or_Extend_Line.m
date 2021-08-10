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

function [Line_AB,new_Point] = Cal_Shorten_or_Extend_Line(Line_AB,delta_L,Point_String)
% Shorten or extend line_AB at point a or b by the increment of offset_L.
% delta_L can be negative.
% Point_String ='A' or 'B'.
%                A                       B
%
%                ●----------------------●
%
% Point_String ='A',delta_L < 0:
%
%                |<---delta_L--->|
%                                ●------●
%              
%              
% Point_String ='B',delta_L > 0:
%
%                                        |<---delta_L--->|
%                ●--------------------------------------●

a_x = Line_AB(1,1);
a_y = Line_AB(1,2);
b_x = Line_AB(2,1);
b_y = Line_AB(2,2);
	

theta = atan2(b_y-a_y,b_x-a_x);

% Move C by offset_delta along the line_AB from B to A or from A to B.
% Case 1: Point_String ='A'
if lower(Point_String) =='a'
	Line_AB(1,1) = Line_AB(1,1)-delta_L*cos(theta);
	Line_AB(1,2) = Line_AB(1,2)-delta_L*sin(theta);
	new_Point = [Line_AB(1,1) Line_AB(1,2)];
% Case 2: Point_String ='B'
elseif  lower(Point_String) =='b'
	Line_AB(2,1) = Line_AB(2,1)+delta_L*cos(theta);
	Line_AB(2,2) = Line_AB(2,2)+delta_L*sin(theta);
    new_Point = [Line_AB(2,1) Line_AB(2,2)];
end
