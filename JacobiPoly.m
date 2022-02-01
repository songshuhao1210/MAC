%%% Jacobi Polynomial for n=0:NN
function y=JacobiPoly(NN,a,b,x)

x=reshape(x,1,length(x));

y(1,:)=1+0*x;

if NN>0
y(2,:)=0.5*(a-b+(a+b+2)*x);

for n=1:1:NN-1
    a1=2*(n+1)*(n+a+b+1)*(2*n+a+b);
    a2=(2*n+a+b+1)*(a^2-b^2);
    a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
    a4=2*(n+a)*(n+b)*(2*n+a+b+2);
    y(n+2,:)=((a2+a3*x).*y(n+1,:)-a4*y(n,:))/a1;
    
end
end
end
