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

function [angle_AB_BC] = Cal_Angle_of_AB_and_BC(Line_AB,Line_BC)
% Calculates the included angle of line AB and BC.
% If the direction is clockwise or anti-clockwise, controlled by Direc_String.
%
%                                         B                       A
%                                         ●----------------------●
%                                        -
%           clockwise                  -        angle_AB_B
%                \   angle_AB_B      -
%                 \                -    
%                  \             -
%                   \          -    
%                    \       ●    
%                    -/     C

if Line_AB(2,1) ~= Line_BC(1,1) | Line_AB(2,2) ~= Line_BC(1,2)
	disp(['*****************************************************'])
	disp(['Attention should be paid here!'])
	disp(['This situation seldom appears in Cal_Angle_of_AB_and_BC.m!'])
	disp(['*****************************************************'])
end   

a_x = Line_AB(1,1);
a_y = Line_AB(1,2);
b_x = Line_AB(2,1);
b_y = Line_AB(2,2);
c_x = Line_BC(2,1);
c_y = Line_BC(2,2);	

angle     = acos(dot([a_x-b_x,a_y-b_y],[c_x-b_x,c_y-b_y])/(norm([a_x-b_x,a_y-b_y])*norm([c_x-b_x,c_y-b_y])));

Direction = cross([a_x-b_x,a_y-b_y,0],[c_x-b_x,c_y-b_y,0]);

if Direction(3) >= 0
    angle_AB_BC = angle;
elseif Direction(3) < 0
    angle_AB_BC = 2*pi-angle;
end

% angle_AB_BC_degree = angle_AB_BC*180/pi;