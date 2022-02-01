function [ Tmk ] = TTT( k,m,x )
% the k-th derivative of T_m(x), 
% the m-th chebyshev polynomial evaluated at x
% x should be column vector
% translated from fortran77 routine developed by Rainer Hollerbach
% Yufeng Lin at ETH Zurich
% Feb 8 2012

sx=length(x);
A=zeros(sx,m+1);
A(:,m+1)=1.0;
for j=1:k
    C=zeros(sx,m+2);
    for i=m:-1:1
        C(:,i)=C(:,i+2)+(2*i)*A(:,i+1); % maybe 2*i+0
    end
    C(:,1)=C(:,1)/2.0;
    A(:,1:m+1)=C(:,1:m+1);
end
B=zeros(sx,m+3);
for i=m+1:-1:1
    B(:,i)=2.0*x.*B(:,i+1)-B(:,i+2)+A(:,i);
end
Tmk=(A(:,1)+B(:,1)-B(:,3))/2.0;
end

