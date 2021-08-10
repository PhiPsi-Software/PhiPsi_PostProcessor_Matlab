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

function [Kesi,Yita] = Cal_KesiYita_by_Coors_3D(X,Y)
% This function calculate the local coordinates kesi and yita by the coordinates of the point inside the element.

global G_X_NODES G_Y_NODES

% Firstly, get the number of the element.

[c_elem] = Cal_Ele_Num_by_Coors(X,Y);	

x1 = G_X_NODES(1,c_elem);
x2 = G_X_NODES(2,c_elem);
x3 = G_X_NODES(3,c_elem);
x4 = G_X_NODES(4,c_elem);

y1 = G_Y_NODES(1,c_elem);
y2 = G_Y_NODES(2,c_elem);
y3 = G_Y_NODES(3,c_elem);
y4 = G_Y_NODES(4,c_elem);

a1 = 0.25 * (-x1+x2+x3-x4);
a2 = 0.25 * ( x1-x2+x3-x4);
a3 = 0.25 * (-x1-x2+x3+x4);
a4 = 0.25 * ( x1+x2+x3+x4);

b1 = 0.25 * (-y1+y2+y3-y4);
b2 = 0.25 * ( y1-y2+y3-y4);
b3 = 0.25 * (-y1-y2+y3+y4);
b4 = 0.25 * ( y1+y2+y3+y4);

% ----------------------------------
% ----------- Option 1 -------------
% ----------------------------------
%options = optimoptions('fsolve','Display','off','TolFun',1e-10);  % For version higher than Matlab R2012
options = optimset('Display','off','TolFun',1e-10);       % For version lower than Matlab R2012

% Attention: no space are allowed in this line.
Eq = @(x) [a4+a3*x(2)+a1*x(1)+a2*x(1)*x(2)-X,b4+b3*x(2)+b1*x(1)+b2*x(1)*x(2)-Y];
pp=fsolve(Eq, [0 0],options);
Kesi=pp(1);
Yita=pp(2);

% ------------------------------------------------------------
% ----------- Option 2: can not be compiled into exe ---------
% ------------------------------------------------------------
% syms kesi yita;
% [kesi,yita] = solve(a4 + a3*yita + a1*kesi + a2*kesi*yita -X,b4 + b3*yita + b1*kesi + b2*kesi*yita -Y,kesi,yita);

% if size(kesi,1) ==2
	% kesi_1 = eval(kesi(1,1));
	% kesi_2 = eval(kesi(2,1));
	% yita_1 = eval(yita(1,1));
	% yita_2 = eval(yita(2,1));

	% if (kesi_1 >= -1) &&  (kesi_1 <=  1)  && (yita_1 >= -1) &&  (yita_1 <=  1) 
		% Kesi = kesi_1;
		% Yita = yita_1;
	% elseif (kesi_2 >= -1) &&  (kesi_2 <=  1)  && (yita_2 >= -1) &&  (yita_2 <=  1) 
		% Kesi = kesi_2;
		% Yita = yita_2;
	% else
		% Kesi = 0.0;
		% Yita = 0.0;
	% end
% elseif size(kesi,1) ==1
    % Kesi = eval(kesi);
	% Yita = eval(yita);
% end
% Kesi
% Yita