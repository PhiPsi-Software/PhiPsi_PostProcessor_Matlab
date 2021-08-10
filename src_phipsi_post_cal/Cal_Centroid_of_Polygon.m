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

function [center_x,center_y] = Cal_Centroid_of_Polygon(X,Y)
% Get the coordinates of the centroid of the polygon.
% The quadrangle is composed of four nodes whose coordinates are X and Y.
% X and Y should be arranged anticlockwise! 
% About the algorithm: http://en.wikipedia.org/wiki/Centroid

% For triangle, the centroid can be get by:
% center_x = (X(1)+X(2)+X(3)) / 3;
% center_y = (Y(1)+Y(2)+Y(3)) / 3;

n = size(X,1);
X(n+1) = X(1);                     % This is very important!
Y(n+1) = Y(1);

% Calculate the area of the polygon
A = 0;
for i = 1:n
    A = A + X(i)*Y(i+1)-X(i+1)*Y(i);
end
A = 0.5*A;

% X location of the centroid of the polygon
center_x = 0;
for i = 1:n
    center_x = center_x + (X(i)+X(i+1)) * (X(i)*Y(i+1)-X(i+1)*Y(i));
end
center_x = center_x/6/A;

% Y location of the centroid of the polygon
center_y = 0;
for i = 1:n
    center_y = center_y + (Y(i)+Y(i+1))*(X(i)*Y(i+1)-X(i+1)*Y(i));
end
center_y = center_y/6/A;

