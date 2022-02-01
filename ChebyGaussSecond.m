function [x,w]=ChebyGaussSecond(N,a,b)

%%% calculate nodes and weights for Chebysheve-Gausse integral (second
%%% type)
%%% use formula 25.4.40-41, page 889 in A.S Handbook

nipi=((1:1:N)*pi)';

x=cos(nipi/(N+1));
x=(a+b)/2+(b-a)/2*x;

w=pi/(N+1)*sin(nipi/(N+1)).^2;
yab=sqrt((x-a).*(b-x));
w=(b-a)^2/4*w./yab;

