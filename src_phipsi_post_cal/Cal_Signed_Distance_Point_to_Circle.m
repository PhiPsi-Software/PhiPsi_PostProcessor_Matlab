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
          
 
function [Signed_Distance] = Cal_Signed_Distance_Point_to_Circle(x0,y0,r,Point_C,S_Distance)
%This function calculates the signed distance from the Point_C to the circle.
%圆内为负,圆外为正,圆上为0

c_Dis = sqrt((x0-Point_C(1))^2 + (y0-Point_C(2))^2);

if c_Dis < r
  Signed_Distance = -abs(c_Dis-r);
elseif c_Dis > r
  Signed_Distance = abs(c_Dis-r);
end

if abs(c_Dis - r) <= 1.0e-8*r
  Signed_Distance = 0.0;
  %Signed_Distance = -1.0
  %Signed_Distance = abs(c_Dis-r);
end
