% Written By: Shi Fang, 2014
% Website: phipsi.top
% Email: phipsi@sina.cn

function sum = hyp2f1(a,b,c,z)
%HYP2F1 Gauss hypergeometric function.
%   h = HYP2F1(a,b,c,z) evaluates the Gauss hypergeometric function 
%   2F1(a,b;c;z) for real value parameters a,b,c and argument z.
%

sum=1;
for n=0:90
    A=a; B=b;C=c;K=1;Z=z^(n+1);
    for k=1:n
        A=A*(a+k);B=B*(b+k);C=C*(c+k);K=K*(k+1);
    end
    %[n+1 A B C K];
    sum=sum+A*B*Z/C/K;
end % HYP2F1;