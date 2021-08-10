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

function Animate_Main(Real_num_iteration)
% This function generates animations(gif files) of deformations, contours, vectors and so on.

global Key_Animation num_Hole Itera_Num


% Generate animations of nodal displacement contours
if Key_Animation(1)==1 | Key_Animation(1)==2
    Animate_Node_and_Gauss_Disp(Real_num_iteration)
end

% Generate animations of nodal stress contours
if Key_Animation(2)==1 | Key_Animation(2)==2
    Animate_Node_and_Gauss_Stress(Real_num_iteration)
end

% Generate deformation animations
if Key_Animation(3)==1
    Animate_Deformation(Real_num_iteration)
end

% 场变量动画绘制
if Key_Animation(4)~=0
    Animate_Fd_Value(Real_num_iteration)
end

if max(Key_Animation) ~= 0
    disp('    Animations generated.')
end
