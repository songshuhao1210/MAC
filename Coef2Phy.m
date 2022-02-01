function [Vel,Vor,Helicity]=Coef2Phy(ri,m,N,L,EqS,anl,cnl,r,theta)
ic=sqrt(-1);
LL=L+m-1;
dxr=2/(1-ri);
DM1=chebdif(N+1,1)*dxr;

li_Pol=m+EqS:2:LL;
lpol1=li_Pol.*(li_Pol+1);
bnl=DM1*anl+2*bsxfun(@times,anl,1./r);
bnl=bsxfun(@times,bnl,1./lpol1);
if m+EqS==0
    bnl(:,1)=0;
end

Vor_anbn=bsxfun(@times,anl,1./r)-bsxfun(@times,DM1*bnl,r)-2*bnl;
Vor_cn=bsxfun(@times,DM1*cnl,r)+2*cnl;

for li=m:LL
    tmp_legendre=zeros(li+2,length(theta));
    tmp_legendre(1:li+1,:)=legendre(li,cos(theta),'norm')/sqrt(2*pi);
    Plm=(-1)^m*tmp_legendre(m+1,:);
    Norm_nm=sqrt((li+m+1)*(li-m));
    Plm1=(-1)^(m+1)*tmp_legendre(m+2,:);

    if (mod(li-m,2)==EqS)
        a_ind=(li-m+2-EqS)/2;
        Ylm_an(a_ind,:)=Plm;
        icm_Ylm_bn(a_ind,:)=ic*m*Plm./sin(theta);
        Ylm_dth_bn(a_ind,:)=m*cot(theta).*Plm+Norm_nm*Plm1;
    else
        c_ind=(li-m+1+EqS)/2;           
        Ylm_cn(c_ind,:)=Plm*li*(li+1);
        icm_Ylm_cn(c_ind,:)=ic*m*Plm./sin(theta);
        Ylm_dth_cn(c_ind,:)=m*cot(theta).*Plm+Norm_nm*Plm1;
    end
end

u_r=anl*Ylm_an;
u_t=(bnl*Ylm_dth_bn+cnl*icm_Ylm_cn);
u_t=bsxfun(@times,u_t,r);
u_p=(bnl*icm_Ylm_bn-cnl*Ylm_dth_cn);
u_p=bsxfun(@times,u_p,r);

Vel(:,:,1)=u_r;
Vel(:,:,2)=u_t;
Vel(:,:,3)=u_p;

Vor_r=cnl*Ylm_cn;
Vor_t=Vor_anbn*icm_Ylm_bn+Vor_cn*Ylm_dth_cn;
Vor_p=-Vor_anbn*Ylm_dth_bn+Vor_cn*icm_Ylm_cn;

Vor(:,:,1)=Vor_r;
Vor(:,:,2)=Vor_t;
Vor(:,:,3)=Vor_p;


Helicity=(u_r).*conj(Vor_r)+(u_t).*conj(Vor_t)+(u_p).*conj(Vor_p);
%Helicity_R=real(u_r).*real(Vor_r)+real(u_t).*real(Vor_t)+real(u_p).*real(Vor_p);
%Helicity_I=imag(u_r).*imag(Vor_r)+imag(u_t).*imag(Vor_t)+imag(u_p).*imag(Vor_p);


