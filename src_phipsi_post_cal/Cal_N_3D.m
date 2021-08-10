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

function [N] = Cal_N_3D(kesi,yita,zeta)
% This function calculates N matrix for 3D problems.

N(1) = (1-kesi)*(1-yita)*(1-zeta)/8;    %kesi=-1,yita=-1,zeta=-1
N(2) = (1+kesi)*(1-yita)*(1-zeta)/8;    %kesi= 1,yita=-1,zeta=-1
N(3) = (1+kesi)*(1+yita)*(1-zeta)/8;    %kesi= 1,yita= 1,zeta=-1
N(4) = (1-kesi)*(1+yita)*(1-zeta)/8;    %kesi=-1,yita= 1,zeta=-1
N(5) = (1-kesi)*(1-yita)*(1+zeta)/8;    %kesi=-1,yita=-1,zeta= 1
N(6) = (1+kesi)*(1-yita)*(1+zeta)/8;    %kesi= 1,yita=-1,zeta= 1
N(7) = (1+kesi)*(1+yita)*(1+zeta)/8;    %kesi= 1,yita= 1,zeta= 1
N(8) = (1-kesi)*(1+yita)*(1+zeta)/8;    %kesi=-1,yita= 1,zeta= 1
	  




