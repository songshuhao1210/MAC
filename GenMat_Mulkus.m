function [AA,BB]=GenMat_Mulkus(ri,m,L,N,Omega,Ekman,Le,Pm,EqS,BCuo,BCui,BCbo,BCbi)
% Base on fomulas in Ricon & Rieutord 2003

% ri-- radius ratio
% m -- Azimuthal wave number
% L -- Resolution in latitude, should be even?
% N -- Resolution in radius
% Omega -- Mean rotation rate
% Ekman -- Ekman number
% Le -- Lehnert number: B0/(sqrt(rho mu) \Omega R)
% Pm -- Magnetic Prandtl number
% EqS -- 0-- Symmetric 1--AntiSymmetric around the equator
% BCuo -0 -- noslip 1 --stress free on the outer boundary condition
% BCui -0 -- noslip 1 --stress free on the inner boundary condition
% BCbi -0 --- insulating 1-- conducting for the magnetic on the inner
% boundary, the outer is always insulating

LL=L+m-1; % Spherical harmonics degree
N1=N+1;
ic=sqrt(-1);

%% Gauss-Lobatto colcation points
x=cos(pi*(0:N)/N)'; %[-1,1]
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
dxr=2/(1-ri);

% Differention Matrix of Chebyshev Polynomial in radial direction
DM0=eye(N1);
DM=chebdif(N1,4);
DM1=DM(:,:,1)*dxr;
DM2=DM(:,:,2)*dxr^2;
DM3=DM(:,:,3)*dxr^3;
DM4=DM(:,:,4)*dxr^4;
clear DM;


LmN=L*N1;

BB=sparse(2*LmN,2*LmN);
AA=sparse(2*LmN,2*LmN);

for li=m:LL
    qlm=sqrt((li^2-m^2)/(4*li^2-1));
    qlm1=sqrt(((li+1)^2-m^2)/(4*(li+1)^2-1));
     
    Alm=1/li^2*qlm;
    Alm1=1/(li+1)^2*qlm1;
    Blm=(li^2-1)*qlm;
    Blm1=(li^2+2*li)*qlm1;
       
    ll1=li*(li+1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mod(li-m,2)==EqS)
      %%%%% Poloidal velocity u_n, toroidal magnetic field, c_n
      Iner=bsxfun(@times,DM2,r)+4*DM1-(ll1-2)*bsxfun(@times,DM0,1./r); %u_n
      Visc=bsxfun(@times,DM4,r)+8*DM3-2*(ll1-6)*bsxfun(@times,DM2,1./r) ...
            -4*ll1*bsxfun(@times,DM1,1./r.^2)-(2*ll1-ll1^2)*bsxfun(@times,DM0,1./r.^3); % u_n
      Visc=Ekman*Visc;
      Lor=Le^2*ic*m*(1-2/ll1)*Iner; %a_n
         
      
      lmdB=DM0; % a_n
      DiffB=DM2+bsxfun(@times,DM1,4./r)+(2-ll1)*bsxfun(@times,DM0,1./r.^2); %a_n
      DiffB=Ekman/Pm*DiffB;
      
      Adve=ic*m*DM0; %u_n
       
      Iner(1:2,:)=0; Iner(N:N1,:)=0;
      Visc(1:2,:)=0; Visc(N:N1,:)=0;
      Lor(1:2,:)=0; Lor(N:N1,:)=0;
      
      lmdB(1,:)=0;  lmdB(N1,:)=0;
      DiffB(1,:)=0;  DiffB(N1,:)=0;
      Adve(1,:)=0;  Adve(N1,:)=0;
      
      
      if li>m
          zLeft=Blm*(DM1-(li-1)*bsxfun(@times,DM0,1./r)); 
          CorLeft=2*Omega*zLeft; % w_n-1
          LorLeft=-2*Le^2*zLeft; % c_n-1
          
          %%% reserve zeros for the B.C.
          CorLeft(1:2,:)=0; CorLeft(N:N1,:)=0;
          LorLeft(1:2,:)=0; LorLeft(N:N1,:)=0;
      end
      if li<LL
          zRight=Blm1*(DM1+(li+2)*bsxfun(@times,DM0,1./r));
          CorRight=2*Omega*zRight; % w_n+1
          LorRight=-2*Le^2*zRight; % c_n+1
                    
          CorRight(1:2,:)=0; CorRight(N:N1,:)=0;
          LorRight(1:2,:)=0; LorRight(N:N1,:)=0;
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
          %% Insulating boundary on the inner ra'-(l-1)a=0
          DiffB(N1,:)=r(N1)*DM1(N1,:)-(li-1)*DM0(N1,:);
      elseif BCbi==1
          %%%%% Conducting boundary on the inner
          DiffB(N1,:)=DM0(N1,:); %a_n=0
          DiffB(N,:)=r(N1)*DM2(N1,:)+4*DM1(N1,:); %(rb_n=0)->r^2a_n''+4ra_n'+2a_n
          lmdB(N,:)=0;   
      else
          warning('Please check the boundary conditions')
      end      
      
      if BCbo==0
          % Insulating boundary condition for B ra'+(l+2)a=0
         DiffB(1,:)=r(1)*DM1(1,:)+(li+2)*DM0(1,:);
      elseif BCbo==1
         %%%%%%% Conducting boundary condition for B a_n=0
         DiffB(1,:)=DM0(1,:); %a_n=0
         DiffB(2,:)=r(1)*DM2(1,:)+4*DM1(1,:); %(rb_n=0)->r^2a_n''+4ra_n'+2a_n
         lmdB(2,:)=0;  
      else
         warning('Please check the boundary conditions') 
      end
      
   
    else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      Iner=DM0; %w_n
      Visc=DM2+bsxfun(@times,DM1,2./r)-ll1*bsxfun(@times,DM0,1./r.^2); % w_n
      Visc=Ekman*Visc;
      Lor=Le^2*ic*m*(1-2/ll1)*Iner; %c_n
     
      
      lmdB=DM0; % c_n
      DiffB=DM2+bsxfun(@times,DM1,2./r)-ll1*bsxfun(@times,DM0,1./r.^2); %c_n
      DiffB=Ekman/Pm*DiffB;
      Adve=ic*m*DM0; %w_n
      
      Iner(1,:)=0; Iner(N1,:)=0;
      Visc(1,:)=0; Visc(N1,:)=0;
      Lor(1,:)=0; Lor(N1,:)=0;
      lmdB(1,:)=0;  lmdB(N1,:)=0;
      DiffB(1,:)=0;  DiffB(N1,:)=0;
      Adve(1,:)=0;  Adve(N:N1,:)=0;
      
      
      
      if li>m
          zLeft=Alm*(bsxfun(@times,DM1,r)-(li-2)*DM0);
          CorLeft=-2*Omega*zLeft; %u_n-1
          LorLeft=2*Le^2*zLeft; %a_n-1
          
          %%% reserve zeros for the B.C.
          CorLeft(1,:)=0; CorLeft(N1,:)=0;
          LorLeft(1,:)=0; LorLeft(N1,:)=0;
      end
      if li<LL
          zRight=Alm1*((li+3)*DM0+bsxfun(@times,DM1,r));
          CorRight=-2*Omega*zRight; % u_n+1
          LorRight=2*Le^2*zRight; % a_n+1
          CorRight(1,:)=0; CorRight(N1,:)=0;
          LorRight(1,:)=0; LorRight(N1,:)=0;
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
          %% Insulating boundary on the inner boundary
          DiffB(N1,:)=DM0(N1,:);
      elseif BCbi==1
          %% Conducting boundary condition for B rcn'+cn=0 at r=ri
          DiffB(N1,:)=r(N1)*DM1(N1,:)+DM0(N1,:);
      else
         warning('Please check the boundary condition');
      end
      
      if BCbo==0
          %%%%%%% Insulating boundary condition for B cn=0 at r=1;
          DiffB(1,:)=DM0(1,:);
      elseif BCbo==1
          %% Conducting boundary condition for B rcn'+cn=0 at r=ri
          DiffB(1,:)=r(1)*DM1(1,:)+DM0(1,:);
      else
         warning('Please check the boundary condition');
      end
          
                      
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    A{2*(li-m)+1}=sparse(Iner);
    A{2*(li-m+1)}=sparse(lmdB);
    
    Visc=ic*2*m*Omega/ll1*Iner+Visc;
    
    
    Bdig=[Visc, Lor; Adve,DiffB];
    
    B{li-m+1}=sparse(Bdig);
    
      
    if ((li-m)>0)
        Bdown=[-CorLeft, -LorLeft; zeros(N1,2*N1)];   
        BdiagLow{li-m}=sparse(Bdown);    
    end
    if (li<LL)
        Bup=[-CorRight, -LorRight; zeros(N1,2*N1)];  
        BdiagUp{li-m+1}=sparse(Bup);
    end
end
display('Calculating matrix values')



AA=blkdiag(A{:});
BB=blkdiag(B{:});
BdownT=blkdiag(BdiagLow{:});
BupT=blkdiag(BdiagUp{:});
B1=sparse(2*N1,2*LmN);
B2=sparse(2*LmN-2*N1,2*N1);
BBdown=[B1; BdownT, B2];
BBup=[B2, BupT; B1];
BB=BB+BBdown+BBup;

end
