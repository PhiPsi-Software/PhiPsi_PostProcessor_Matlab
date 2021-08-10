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

function [x,y] = Cal_Coor_by_KesiYita(kesi,yita,X_NODES,Y_NODES)
% Get the global coordinates of the point represented by the local kesi and yita.

N  = 1/4*[(1-kesi).*(1-yita) (1+kesi).*(1-yita) (1+kesi).*(1+yita) (1-kesi).*(1+yita)];

% JM=[0,            1-yita,      yita-kesi,        kesi-1    ;...
	% yita-1,       0,           1+kesi,          -kesi-yita ;...
	% kesi-yita,    -kesi-1,     0,                yita+1    ;...
	% 1-kesi,       yita+kesi,   -yita-1,          0         ];

% detJ=X_NODES*JM*Y_NODES'/8;

x  = N(:,1)*X_NODES(1)+N(:,2)*X_NODES(2)+N(:,3)*X_NODES(3)+N(:,4)*X_NODES(4);

y  = N(:,1)*Y_NODES(1)+N(:,2)*Y_NODES(2)+N(:,3)*Y_NODES(3)+N(:,4)*Y_NODES(4);
