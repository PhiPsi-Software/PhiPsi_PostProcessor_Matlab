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

function [F,dFdx,dFdy] = Cal_F_dFdx_dFdy(r,theta,omega,mu,c_mat_type)
% This function calculates the tip enrichment functions F and their derivative, dFdx and dFdy.

global aveg_area_ele 
global Key_TipEnrich

if r~=0
    r2 = sqrt(r);
else
    r2 = sqrt(aveg_area_ele)*0.1d-4;
	theta = 0.0d0;
end	

% ------------------------------
%         ISO material.
% ------------------------------
if c_mat_type ==1  
	fac = 0.5/r2;
	st2 = sin(theta/2.0);
	ct2 = cos(theta/2.0);
	s3t2= sin(1.5*theta);
	c3t2= cos(1.5*theta);
	st  = sin(theta);
	ct  = cos(theta);

	new_sin=sin(-3*theta/2);
	new_cos=cos(-3*theta/2);
	
	% Tip enrichment functions F.
	
	F(1) = r2*st2;
	F(2) = r2*ct2;
	F(3) = r2*st2*st;
	F(4) = r2*ct2*st;

	% dF1dx1 and dF1dx2
	dF1dx1 = -fac*st2;
	dF1dx2 =  fac*ct2;
	% dF2dx1 and dF2dx2
	dF2dx1 =  dF1dx2;
	dF2dx2 = -dF1dx1;
	% dF3dx1 and dF3dx2
	dF3dx1 = -fac*s3t2*st;
	dF3dx2 =  fac*(st2 + s3t2*ct);
	% dF4dx1 and dF4dx2
	dF4dx1 = -fac*c3t2*st;
	dF4dx2 =  fac*(ct2 + c3t2*ct);
 	
	% ------------------------
	% F(1) = r2*st2;
	% F(2) = r2*ct2;
	% F(3) = r2*new_sin;
	% F(4) = r2*new_cos;
	
	% dF1dx1 and dF1dx2
	% dF1dx1 = -fac*st2;
	% dF1dx2 =  fac*ct2;
	% dF2dx1 and dF2dx2
	% dF2dx1 =  dF1dx2;
	% dF2dx2 = -dF1dx1;
	% dF3dx1 and dF3dx2
	% dF3dx1 =  fac*new_sin*cos(theta)+3*fac*new_cos*sin(theta);
	% dF3dx2 =  fac*new_sin*sin(theta)-3*fac*new_cos*cos(theta);
	% dF4dx1 and dF4dx2
	% dF4dx1 =  fac*new_cos*cos(theta)-3*fac*new_sin*sin(theta);
	% dF4dx2 =  fac*new_cos*sin(theta)+3*fac*new_sin*cos(theta);

	% dx1dx
	dx1dx =  cos(omega);
	% dx2dx
	dx2dx = -sin(omega);
	% dx1dy
	dx1dy =  sin(omega);
	% dx2dy
	dx2dy =  cos(omega);

	% dFdx and dFdy
	dFdx(1) = dF1dx1*dx1dx + dF1dx2*dx2dx;
	dFdy(1) = dF1dx1*dx1dy + dF1dx2*dx2dy;	
	dFdx(2) = dF2dx1*dx1dx + dF2dx2*dx2dx;
	dFdy(2) = dF2dx1*dx1dy + dF2dx2*dx2dy;	
	dFdx(3) = dF3dx1*dx1dx + dF3dx2*dx2dx;
	dFdy(3) = dF3dx1*dx1dy + dF3dx2*dx2dy;	
	dFdx(4) = dF4dx1*dx1dx + dF4dx2*dx2dx;
	dFdy(4) = dF4dx1*dx1dy + dF4dx2*dx2dy;	
% ------------------------------
%       Orthotropic material.
% ------------------------------
elseif c_mat_type ==2 || c_mat_type ==3
    mu_1 = mu(1);
	mu_2 = mu(3);
	mu_1_real = real(mu_1); mu_1_imag = imag(mu_1);
	mu_2_real = real(mu_2); mu_2_imag = imag(mu_2);
	
	fac = 0.5/r2;
	st2 = sin(theta/2.0);
	ct2 = cos(theta/2.0);
	s3t2= sin(1.5*theta);
	c3t2= cos(1.5*theta);
	st  = sin(theta);
	ct  = cos(theta);
 	sa  = sin(omega);
	ca  = cos(omega);
	
	% theta_k
	Theta_1 = atan2(mu_1_imag*st, ct + mu_1_real*st);
	Theta_2 = atan2(mu_2_imag*st, ct + mu_2_real*st);
	% g_k
	g_1 = sqrt((ct + mu_1_real*st)^2 + (mu_1_imag*st)^2);
	g_2 = sqrt((ct + mu_2_real*st)^2 + (mu_2_imag*st)^2);

	dg1_dt = ((ct +mu_1_real*st)*(-st + mu_1_real*ct) + mu_1_imag^2*st*ct)/g_1; 
	dg2_dt = ((ct +mu_2_real*st)*(-st + mu_2_real*ct) + mu_2_imag^2*st*ct)/g_2; 
	
	U_1 = mu_1_imag*st/(ct + mu_1_real*st); 
	U_2 = mu_2_imag*st/(ct + mu_2_real*st); 
	
	dTheta1_dt = mu_1_imag/(ct + mu_1_real*st)^2/(1+U_1^2);
    dTheta2_dt = mu_2_imag/(ct + mu_2_real*st)^2/(1+U_2^2);
	
	G_1 = sqrt(g_1);   G_2 = sqrt(g_2);
	
	dG1_dt = dg1_dt/(2*G_1);
	dG2_dt = dg2_dt/(2*G_2);
	
	ct_1 = cos(Theta_1/2);   ct_2 = cos(Theta_2/2);
	st_1 = sin(Theta_1/2);   st_2 = sin(Theta_2/2);
	
	% Tip enrichment functions F.
	F(1) = r2*ct_1*G_1;    F(2) = r2*ct_2*G_2;
	F(3) = r2*st_1*G_1;    F(4) = r2*st_2*G_2;

	dct1_dt = -dTheta1_dt/2 * sin(Theta_1/2);
	dct2_dt = -dTheta2_dt/2 * sin(Theta_2/2);
	dst1_dt =  dTheta1_dt/2 * cos(Theta_1/2);
	dst2_dt =  dTheta2_dt/2 * cos(Theta_2/2);
	
	dF1_dr = fac * ct_1 * G_1; dF2_dr = fac * ct_2 * G_2;
	dF3_dr = fac * st_1 * G_1; dF4_dr = fac * st_2 * G_2;
	
	dF1_dt = r2 * (dct1_dt*G_1 + ct_1*dG1_dt);
	dF2_dt = r2 * (dct2_dt*G_2 + ct_2*dG2_dt);
	dF3_dt = r2 * (dst1_dt*G_1 + st_1*dG1_dt);
	dF4_dt = r2 * (dst2_dt*G_2 + st_2*dG2_dt);
	
	dF_dr  = [dF1_dr ,dF2_dr ,dF3_dr ,dF4_dr];
	dF_dt  = [dF1_dt ,dF2_dt ,dF3_dt ,dF4_dt];
	
    dFxy = [ca  -sa;sa  ca]*[ct  -st/r;st  ct/r]*[dF_dr;dF_dt];
	
	% dFdx and dFdy
	dFdx = dFxy(1,:);
	dFdy = dFxy(2,:); 
end
%如果是粘性裂尖
if Key_TipEnrich ==4
	st2 = sin(theta/2.0);
	ct2 = cos(theta/2.0);
    %Tip enrichment functions F.
	F(1) = r*st2;
	
	% F(1) = r2*st2;
	
    F(2) = 0.0;
    F(3) = 0.0;
    F(4) = 0.0;
	%dF1dx1 and dF1dx2
    dF1dx1 = st2;
    dF1dx2 = r/2.0*ct2;
	
	
	% dF1dx1 and dF1dx2
	% fac = 0.5/r2;
	% dF1dx1 = -fac*st2;
	% dF1dx2 =  fac*ct2;
	
	
    %dx1dx
    dx1dx =  cos(omega);
    %dx2dx
    dx2dx = -sin(omega);
    %dx1dy
    dx1dy =  sin(omega);
    %dx2dy
    dx2dy =  cos(omega);
    %dFdx and dFdy
    dFdx(1) = dF1dx1*dx1dx + dF1dx2*dx2dx;
    dFdy(1) = dF1dx1*dx1dy + dF1dx2*dx2dy;
    dFdx(2) = 0.0;
    dFdy(2) = 0.0;
    dFdx(3) = 0.0;
    dFdy(3) = 0.0;
    dFdx(4) = 0.0;
    dFdy(4) = 0.0;
end