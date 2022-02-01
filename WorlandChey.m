function [Wor,dWor,Wor_r]=WorlandChey(NN,l,r)

r=reshape(r,1,length(r));
Wor=zeros(NN+1,length(r));
Wor_r=zeros(NN+1,length(r));
dWor=zeros(NN+1,length(r));
if l==0
    lrl1=0*r;
else
    lrl1=l*r.^(l-1);
end
rl=r.^l;
rlp1=r.^(l+1);

x=2*r.^2-1;

al=-0.5;
be=l-0.5;
P0=1+0*r;
P1=0.5*(al-be+(al+be+2)*x);

al_1=0.5;
be_1=l+0.5;
P0_1=1+0*r;
P1_1=0.5*(al_1-be_1+(al_1+be_1+2)*x);

Pn(1,:)=P0;
Pn(2,:)=P1;
Wor=bsxfun(@times,Pn,rl);
Pn_1(1,:)=0*r;
Pn_1(2,:)=P0_1;
Pn_1(3,:)=P1_1;


for n=1:1:NN-1
    a1=2*(n+1)*(n+al+be+1)*(2*n+al+be);
    a2=(2*n+al+be+1)*(al^2-be^2);
    a3=(2*n+al+be)*(2*n+al+be+1)*(2*n+al+be+2);
    a4=2*(n+al)*(n+be)*(2*n+al+be+2);
    Pn(n+2,:)=((a2+a3*x).*Pn(n+1,:)-a4*Pn(n,:))/a1;
    
    if n<NN-1;
    a1=2*(n+1)*(n+al_1+be_1+1)*(2*n+al_1+be_1);
    a2=(2*n+al_1+be_1+1)*(al_1^2-be_1^2);
    a3=(2*n+al_1+be_1)*(2*n+al_1+be_1+1)*(2*n+al_1+be_1+2);
    a4=2*(n+al_1)*(n+be_1)*(2*n+al_1+be_1+2);
    Pn_1(n+3,:)=((a2+a3*x).*Pn_1(n+2,:)-a4*Pn_1(n+1,:))/a1;
    end
end

Wor=bsxfun(@times,Pn,rl);
Wor_r=bsxfun(@times,Pn,rlp1);

nl=(0:1:NN)'+l;
Pn_1=bsxfun(@times,Pn_1,nl);

dWor=bsxfun(@times,Pn,lrl1)+2*bsxfun(@times,Pn_1,rlp1);

end
%%%% Normalization factor
% for n=0:NN+1
%     n=double(n); l=double(l);
%     hn(n+1)=(n+l)/(2^(2*l+1)*(2*n+l))*factorial(2*n+2*l)*(factorial(n))^2/(factorial(2*n)*(factorial(n+l))^2);
%     hn(n+1)=hn(n+1)*2^l*gamma(n+1/2)*gamma(n+l+1/2)/((2*n+l)*gamma(n+1)*gamma(n+l));
% end
% hn=sqrt(hn);
function [test]=Jacobii(NN,l,r)

test=r;

end



