% Cardinal funciton of Chebyshev Polynomial
% Using formula F.44 on page 570 in Boyd J.P. @Chebyshev and Fourier Spectral Methods@
% Note there is a typo: cj=2 if j=N or 0;

function [ Card ] = Cardinal( j,N,x )

% Formula F.44(a)
% xj=cos(pi*j/N);
% if (j==N | j==0)
%     cj=2;
% else
%     cj=1;
% end
% 
% Card=(-1)^(j+1)*(1-x.^2).*TTT(1,N,x)./(cj*N^2*(x-xj));


% Formula F.44(b)
xj=cos(pi*j/N);
if (j==N | j==0)
    pj=2;
else
    pj=1;
end

Card=zeros(size(x));
for m=0:N
    if (m==N | m==0)
        pm=2;
    else
        pm=1;
    end
    Card=Card+1/pm*TTT(0,m,xj)*TTT(0,m,x);
end
Card=Card*2/N/pj;