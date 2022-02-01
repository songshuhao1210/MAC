function [KE,Diss]=Int_KE_Diss(ri,m,N,L,EqS,anl,cnl)
LL=L+m-1;

x=cos(pi*(0:N)/N)'; %[-1,1]
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
%r(1)=1-1e-6;
dxr=2/(1-ri);
DM1=chebdif(N+1,1)*dxr;

%% Using Gaussian-Chebysheve second type
[r_node,weight]=ChebyGaussSecond(N-1,ri,1.0);
weight_r2=weight.*r_node.^2;
weight_r4=weight.*r_node.^4;
li_Pol=m+EqS:2:LL;
lpol1=li_Pol.*(li_Pol+1);

bnl=DM1*anl+2*bsxfun(@times,anl,1./r);
bnl=bsxfun(@times,bnl,1./lpol1);
if m+EqS==0
    bnl(:,1)=0;
end

%%%%%%%%%%%% Int[r^2 a_n^2]
anl_r=anl(2:N,:);
Int_anl=weight_r2'*(abs(anl_r).^2);

%%%%%%%%%%%% Int[n(n+1)r^4 b_n^2]
bnl_r=bnl(2:N,:);
Int_bnl=weight_r4'*(abs(bnl_r).^2);
Int_bnl=Int_bnl.*lpol1;

%%%%%%%%%%%% Int[n(n+1)r^4 c_n^2]
cnl_r=cnl(2:N,:);
Int_cnl=weight_r4'*(abs(cnl_r).^2);
li=m+1-EqS:2:LL;
lli=li.*(li+1);
Int_cnl=Int_cnl.*lli;

%%%%%%%%%%%% Int[3 r^2 a'_n^2]
dr_anl=DM1(2:N,:)*anl;
Int_3dan=3*weight_r2'*(abs(dr_anl).^2);

%%%%%%%%%%%% Int[n(n+1)(a_n+r^2 b'_n)^2]
dr_bnl=DM1(2:N,:)*bnl;
r2_dr_bnl=bsxfun(@times, dr_bnl, r_node.^2);
Int_anl_dbnl=weight'*(abs(anl_r+r2_dr_bnl).^2);
li=m+EqS:2:LL;
lli=li.*(li+1);
Int_anl_dbnl=Int_anl_dbnl.*lli;

%%%%%%%%%%%% Int[(n-1)n(n+1)(n+2)r^2 b_n^2]
Int_n4_bnl=weight_r2'*(abs(bnl_r).^2);
li=m+EqS:2:LL;
lli=li.*(li+1).*(li+2).*(li-1);
Int_n4_bnl=Int_n4_bnl.*lli;

%%%%%%%%%%%% Int[n(n+1)r^4 c'_n^2]
dr_cnl=DM1(2:N,:)*cnl;
Int_dr_cnl=weight_r4'*(abs(dr_cnl).^2);
li=m+1-EqS:2:LL;
lli=li.*(li+1);
Int_dr_cnl=Int_dr_cnl.*lli;

%%%%%%%%%%%% Int[(n-1)n(n+1)(n+2)r^2 c_n^2]
Int_n4_cnl=weight_r2'*(abs(cnl_r).^2);
li=m+1-EqS:2:LL;
lli=li.*(li+1).*(li+2).*(li-1);
Int_n4_cnl=Int_n4_cnl.*lli;


KE=0.25*(sum(Int_anl+Int_bnl)+sum(Int_cnl)); %% Totoal kinetic energy
Diss=0.5*(sum(Int_3dan+Int_anl_dbnl+Int_n4_bnl)+sum(Int_dr_cnl+Int_n4_cnl)); %%% Viscous dissipation