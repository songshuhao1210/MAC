function [KE,Diss]=Int_KE_Diss_Rieutord(ri,m,N,L,EqS,anl,cnl)
%% calculate total energy and total dissipation
%% RIEUTORD coefficients
%% 18 May 2017 by Yufeng Lin @ Cambridge

LL=L+m-1;
x=cos(pi*(0:N)/N)'; %[-1,1]
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
dxr=2/(1-ri);
DM1=chebdif(N+1,1)*dxr;

%% Using Gaussian-Chebysheve second type
[r_node,weight]=ChebyGaussSecond(N-1,ri,1.0);
weight_r2=weight.*r_node.^2;

%%%%% l*(l+1)
li_Pol=m+EqS:2:LL;
lliP=li_Pol.*(li_Pol+1); % l*(l+1) for an and bn
li_Tor=m+1-EqS:2:LL;
lliT=li_Tor.*(li_Tor+1); % l*(l+1) for cn

%%% bnl
bnl=2*anl+bsxfun(@times,DM1*anl,r);
bnl=bsxfun(@times,bnl,1./lliP);
if m+EqS==0
    bnl(:,1)=0;
end

Sum_an2=sum(abs(anl).^2, 2);
Sum_lli_bn2=(abs(bnl).^2)*lliP';
Sum_lli_cn2=(abs(cnl).^2)*lliT';

%% KE_r=r^2(an^2+l(l+1)(bn^2+cn^2))
KE_r=(Sum_an2+Sum_lli_bn2+Sum_lli_cn2).*r.^2;

Sum_Rn=(abs(cnl).^2)*(lliT'.^2);
Sn=bsxfun(@times,DM1*cnl,r)+cnl;
Sum_Sn=(abs(Sn).^2)*lliT';
Tn=anl-bnl-bsxfun(@times,DM1*bnl,r);
Sum_Tn=(abs(Tn).^2)*lliP';

%%% Diss_r=(l(l+1))^2 cn^2+l(l+1)(cn+r*cn').^2+l(l+1)(an-bn-r*bn')
Diss_r=Sum_Rn+Sum_Sn+Sum_Tn;

KE=0.25*weight'*KE_r(2:N); %% Totoal kinetic energy
Diss=0.5*weight'*Diss_r(2:N); %%% Viscous dissipation

end