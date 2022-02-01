function [AA,BB]=GenMat(ri,m,L,N,Omega,Ekman,Le,N2,En,Em,EqS,MagF,BCuo,BCui,BCbi,N2r)
% Base on fomulas in Ricon & Rieutord 2003

% ri-- radius ratio
% m -- Azimuthal wave number
% L -- Resolution in latitude, should be even?
% N -- Resolution in radius
% Omega -- Mean rotation rate
% Ekman -- Ekman number
% Le -- Lehner

%number: B0/(sqrt(rho mu) \Omega R)
% Pm -- Magnetic Prandtl number
% EqS -- 0-- Symmetric 1--AntiSymmetric around the equator
% MagF -0 -- Dipolar field 1-- uniform vertical field.
% BCuo -0 -- noslip 1 --stress free on the outer boundary condition
% BCui -0 -- noslip 1 --stress free on the inner boundary condition
% BCbi -0 --- insulating 1-- conducting for the magnetic on the inner
% boundary, the outer is always insulating

LL=L+m-1; % Spherical harmonics degree
N1=N+1;
ic=sqrt(-1);
%% Define size of AA,BB
total_index = find_index(L,m,EqS);


%% Gauss-Lobatto colcation points
x=cos(pi*(0:N)/N)'; %[-1,1]vstack((ar,ai)).T
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
dxr=2/(1-ri);

%%%% Radial dependence of the background magnetic field
if MagF==0
    Br=r.^(-3);
    dBr=-3*r.^(-4);
    Bt=0.5*Br;
    dBt=0.5*dBr;
elseif MagF==1
    Br=1+0*r;
    dBr=0*r;
    Bt=-1+0*r;
    dBt=0*r;
else
    warning('Wrong Magnetic Flag');
end
   


%% Differention Matrix of Chebyshev Polynomial in radial direction
DM0=eye(N1);
DM=chebdif(N1,4);
DM1=DM(:,:,1)*dxr;
DM2=DM(:,:,2)*dxr^2;
DM3=DM(:,:,3)*dxr^3;
DM4=DM(:,:,4)*dxr^4;
clear DM;


LmN=L*N1;

BB1=sparse(total_index*N1,total_index*N1);

countA = 1;

for li=m:LL
    qlm=sqrt((li^2-m^2)/(4*li^2-1));
    qlm1=sqrt(((li+1)^2-m^2)/(4*(li+1)^2-1));
     
    Alm=1/li^2*qlm;
    Alm1=1/(li+1)^2*qlm1;
    Blm=(li^2-1)*qlm;
    Blm1=(li^2+2*li)*qlm1;
       
    ll1=li*(li+1);

%%%%%%%%%%% Matrix block
%     Iner=zeros(N1,N1);
%     Visc=zeros(N1,N1);
%     Lor=zeros(N1,N1);
%     lmdB=zeros(N1,N1);
%     DiffB=zeros(N1,N1);
%     Adve=zeros(N1,N1);
%     CorLeft=zeros(N1,N1);
%     LorLeft=zeros(N1,N1);
%     AdveLeft=zeros(N1,N1);
%     CorRight=zeros(N1,N1);
%     LorRight=zeros(N1,N1);
%     AdveRight=zeros(N1,N1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mod(li-m,2)==EqS)
      %%%%% Poloidal velocity u_n, toroidal magnetic field, c_n
      Iner=bsxfun(@times,DM2,r)+4*DM1-(ll1-2)*bsxfun(@times,DM0,1./r); %u_n
      Visc=bsxfun(@times,DM4,r)+8*DM3-2*(ll1-6)*bsxfun(@times,DM2,1./r) ...
            -4*ll1*bsxfun(@times,DM1,1./r.^2)-(2*ll1-ll1^2)*bsxfun(@times,DM0,1./r.^3); % u_n
      Visc=Ekman*Visc;
      Lor=ic*m*(bsxfun(@times,DM2,Br)+bsxfun(@times,DM1,2*Br./r+dBr)-bsxfun(@times,DM0,ll1*dBt./r-dBr./r-ll1*Bt./r.^2)); %c_n
      Lor=Le^2*Lor;
      
      lmdB=DM0; % c_n
      DiffB=DM2+bsxfun(@times,DM1,2./r)-ll1*bsxfun(@times,DM0,1./r.^2); %c_n
      DiffB=Em*DiffB;
      Adve=ic*m*(-1/ll1^2*bsxfun(@times,DM2,Br.*r)-1/ll1^2*bsxfun(@times,DM1,4*Br+dBr.*r) ...
          +bsxfun(@times,DM0,dBt/ll1-Bt./r/ll1-2*dBr/ll1^2-2*Br./r/ll1^2)); %u_n
       
      Iner(1:2,:)=0; Iner(N:N1,:)=0;
      Visc(1:2,:)=0; Visc(N:N1,:)=0;
      Lor(1:2,:)=0; Lor(N:N1,:)=0;
      
      lmdB(1,:)=0;  lmdB(N1,:)=0;
      DiffB(1,:)=0;  DiffB(N1,:)=0;
      Adve(1,:)=0;  Adve(N1,:)=0;
      
      %addition
      Fory=-ll1*N2*DM0;%t_n ..
      lmdC=DM0;% t_n ..
      %Ad1=-DM0;% u_n.. %£¿£¿£¿ÊÇ¼ÓºÅÂð
      Ad1 = bsxfun(@times,DM0,N2r);
      Ad2=En*(bsxfun(@times,DM2,1) +bsxfun(@times,DM1,2./r) -ll1*bsxfun(@times,DM0,1./r.^2));
      
      Fory(1:2,:)=0; Fory(N:N1,:)=0;
      lmdC(1,:)=0; lmdC(N1,:)=0;
      Ad1(1,:)=0;Ad1(N1,:)=0;
      Ad2(1,:)=0;Ad2(N1,:)=0;
      
      %Ad2(1,:)=DM0(1,:);Ad2(N1,:)=DM0(N1,:); ???DM1
      Ad2(1,:)=DM1(1,:);Ad2(N1,:)=DM1(N1,:);
       
      
      if li>m
          CorLeft=Omega*Blm*(DM1-(li-1)*bsxfun(@times,DM0,1./r)); % w_n-1
          LorLeft=ll1*Alm*(bsxfun(@times,DM3,Br.*r)+bsxfun(@times,DM2,li*Bt+dBr.*r+6*Br)...
              -bsxfun(@times,DM1,(li+2)*(li-3)*Br./r-4*dBr-4*li*Bt./r) ...
              -(li+1)*(li-2)*bsxfun(@times,DM0,li*Bt./r.^2+dBr./r)); %a_n-1
          LorLeft=Le^2*LorLeft;
          AdveLeft=Alm*li*(li-1)*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,li*Bt./r+Br./r+dBr)); % w_n-1
          
          %%% reserve zeros for the B.C.
          CorLeft(1:2,:)=0; CorLeft(N:N1,:)=0;
          LorLeft(1:2,:)=0; LorLeft(N:N1,:)=0;
          AdveLeft(1,:)=0; AdveLeft(N1,:)=0;
      end
      if li<LL
          CorRight=Omega*Blm1*(DM1+(li+2)*bsxfun(@times,DM0,1./r)); % w_n+1
          LorRight=ll1*Alm1*(bsxfun(@times,DM3,Br.*r)-bsxfun(@times,DM2,(li+1)*Bt-dBr.*r-6*Br) ...
              -bsxfun(@times,DM1,(li+4)*(li-1)*Br./r-4*dBr+4*(li+1)*Bt./r) ...
              +li*(li+3)*bsxfun(@times,DM0,(li+1)*Bt./r.^2-dBr./r)); % a_n+1
          LorRight=Le^2*LorRight;
          AdveRight=(li+1)*(li+2)*Alm1*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,Br./r+dBr-(li+1)*Bt./r)); % w_n+1
          
          CorRight(1:2,:)=0; CorRight(N:N1,:)=0;
          LorRight(1:2,:)=0; LorRight(N:N1,:)=0;
          AdveRight(1,:)=0; AdveRight(N1,:)=0;          
      end
      
      if BCuo==0 
          %%% no-slip on the outer 
          Visc(1,:)=DM0(1,:);    % u_n=0
          Visc(2,:)=r(1)*DM1(1,:)+2*DM0(1,:);  % v_n=0-> r u_n'+2u_n=0
      elseif BCuo==1
          %%% stress-free on the outer
          Visc(1,:)=DM0(1,:);    % u_n=0
          Visc(2,:)=r(1)^2*DM2(1,:)+2*r(1)*DM1(1,:)+(ll1-2)*DM0(1,:); 
          % (v_n/r)'+u_n/r^2=0-> u_n''+2u_n'/r-(2+n(n+1))u_n/r^2=0 
      else
          warning('Please check the boundary condition');
      end
        
      
      if BCui==0 
          %%% no-slip on the inner boundary
          Visc(N1,:)=DM0(N1,:); % u_n=0 on r=r_i
          Visc(N,:)=r(N1)*DM1(N1,:)+2*DM0(N1,:);  % v_n=0-> r u_n'+2u_n=0
      elseif BCui==1
          %%% stree free on the inner boundary
          Visc(N1,:)=DM0(N1,:); % u_n=0 on r=r_i
          Visc(N,:)=r(N1)^2*DM2(N1,:)+2*r(N1)*DM1(N1,:)+(ll1-2)*DM0(N1,:);  
          % (v_n/r)'+u_n/r^2=0-> u_n''+2u_n'/r-(2+n(n+1))u_n/r^2=0 
      else
          warning('Please check the boundary condition');
      end
      
      if BCbi==0
          %% Insulating boundary on the inner boundary
          DiffB(N1,:)=DM0(N1,:);
      elseif BCbi==1
          %% Conducting boundary condition for B rcn'+cn=0 at r=ri
          DiffB(N1,:)=r(N1)*DM1(N1,:)+DM0(N1,:);
      else
         warning('Please check the boundary condition');
      end
                 
      %%%%%%% Insulating boundary condition for B cn=0 at r=1;
      DiffB(1,:)=DM0(1,:);
    else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      Iner=DM0; %w_n
      Visc=DM2+bsxfun(@times,DM1,2./r)-ll1*bsxfun(@times,DM0,1./r.^2); % w_n
      Visc=Ekman*Visc;
      Lor=ic*m/ll1^2*(-bsxfun(@times,DM2,Br.*r)-4*bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,(ll1-2)*Br./r)); %a_n
      Lor=Le^2*Lor;
      
      lmdB=DM0; % a_n
      DiffB=DM2+bsxfun(@times,DM1,4./r)+(2-ll1)*bsxfun(@times,DM0,1./r.^2); %a_n
      DiffB=Em*DiffB;
      Adve=ic*m*bsxfun(@times,DM0,Br./r); %w_n
      
      Iner(1,:)=0; Iner(N1,:)=0;
      Visc(1,:)=0; Visc(N1,:)=0;
      Lor(1,:)=0; Lor(N1,:)=0;
      lmdB(1,:)=0;  lmdB(N1,:)=0;
      DiffB(1,:)=0;  DiffB(N1,:)=0;
      Adve(1,:)=0;  Adve(N1,:)=0;
      
      
      
      if li>m
          CorLeft=-Omega*Alm*(bsxfun(@times,DM1,r)-(li-2)*DM0); %u_n-1
          LorLeft=li*(li-1)*Alm*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,(Br+li*Bt)./r)); %c_n-1
          LorLeft=Le^2*LorLeft;
          AdveLeft=ll1*Alm*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,(2*Br+li*Bt)./r)); % u_n-1
                    %%% reserve zeros for the B.C.
          CorLeft(1,:)=0; CorLeft(N1,:)=0;
          LorLeft(1,:)=0; LorLeft(N1,:)=0;
          AdveLeft(1,:)=0; AdveLeft(N1,:)=0;
      end
      if li<LL
          CorRight=-Omega*Alm1*((li+3)*DM0+bsxfun(@times,DM1,r)); % u_n+1
          LorRight=(li+1)*(li+2)*Alm1*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,(Br-(li+1)*Bt)./r)); % c_n+1
          LorRight=Le^2*LorRight;
          AdveRight=ll1*Alm1*(bsxfun(@times,DM1,Br)+bsxfun(@times,DM0,(2*Br-(li+1)*Bt)./r)); % u_n+1
          CorRight(1,:)=0; CorRight(N1,:)=0;
          LorRight(1,:)=0; LorRight(N1,:)=0;
          AdveRight(1,:)=0;AdveRight(N1,:)=0;          
      end
      
          
      
      %%%%%%%%%%% Set B.C. for w_n in Vis matrix rw'-w=0 
      if BCuo==0
          %%% no-slip on the outer 
          Visc(1,:)=DM0(1,:); % w_n=0
      elseif BCuo==1
          %%% stress free on the outer bounda
          Visc(1,:)=r(1)*DM1(1,:)-DM0(1,:); %(w_n/r)'=0  
      else
          warning('Please check the boundary conditions')
      end
      
      if BCui==0
          %%% no-slip on the inner boundary
          Visc(N1,:)=DM0(N1,:);  %% w_n=0
      elseif BCui==1
          %%% stress-free on the inner boundary
          Visc(N1,:)=r(N1)*DM1(N1,:)-DM0(N1,:);  %(w_n/r)'=0  
      else
          warning('Please check the boundary conditions')
      end
      
      
      if BCbi==0
          %% Insulating boundary on the inner ra'-(l-1)a=0
          DiffB(N1,:)=r(N1)*DM1(N1,:)-(li-1)*DM0(N1,:);
      elseif BCbi==1
          %%%%% Conducting boundary on the inner
          DiffB(N1,:)=DM0(N1,:); %a_n=0
          DiffB(N,:)=r(N1)*DM2(N1,:)+4*DM1(N1,:); %(rb_n=0)->r^2a_n''+4ra_n'+2a_n
          lmdB(N,:)=0;Adve(N,:)=0;
          AdveLeft(N,:)=0;AdveRight(N,:)=0;          
      else
          warning('Please check the boundary conditions')
      end      
      %%%%%%% Insulating boundary condition for B ra'+(l+2)a=0
      DiffB(1,:)=r(1)*DM1(1,:)+(li+2)*DM0(1,:);
                      
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    A{2*(li-m)+1}=sparse(Iner);
    A{2*(li-m+1)}=sparse(lmdB);

    
    Visc=ic*m*Omega/ll1*Iner+Visc; %!
    Bdig=[Visc, Lor; Adve,DiffB];
    B{li-m+1}=sparse(Bdig);
    
    if (mod(li-m,2)==EqS)
        A{2*LmN+countA} = sparse(lmdC);
        BB1(2*LmN+(countA-1)*N1+1:2*LmN+countA*N1,2*(li-m)*N1+1:(2*(li-m)+1)*N1) = Ad1;
        BB1(2*(li-m)*N1+1:(2*(li-m)+1)*N1,2*LmN+(countA-1)*N1+1:2*LmN+countA*N1) = Fory;
        B{2*LmN+countA} = sparse(Ad2);
        countA = countA + 1;
    end
    
    

    
      
    if ((li-m)>0)
        Bdown=[-CorLeft, LorLeft; AdveLeft,zeros(N1,N1)];   
        BdiagLow{li-m}=sparse(Bdown);    
    end
    if (li<LL)
        Bup=[-CorRight, LorRight; AdveRight,zeros(N1,N1)];  
        BdiagUp{li-m+1}=sparse(Bup);
    end
end
display('Calculating matrix values')



AA=blkdiag(A{:});
BB=blkdiag(B{:});
BdownT=blkdiag(BdiagLow{:});
BupT=blkdiag(BdiagUp{:});

B1 = sparse(2*N1,total_index*N1);
B2 = sparse(2*LmN-2*N1,(total_index+2)*N1-2*LmN);
B3 = sparse(total_index*N1-2*LmN,total_index*N1);
BBdown = [B1;BdownT,B2;B3];

B1 = sparse(2*LmN-2*N1,2*N1);
B2 = sparse(2*LmN-2*N1,total_index*N1-2*LmN);
B3 = sparse((total_index+2)*N1-2*LmN,total_index*N1);
BBup = [B1,BupT,B2;B3];

BB = BB + BBdown + BBup + BB1;





end
