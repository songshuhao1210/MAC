% Main Program
% by Yufeng Lin 11 Nov 2016
% based on Ricon & Rieturd,2003

clear all
close all


L=20 % resolution in latitude, Should be even
N=10 %N>4 Chebyshev polymonical numbers，collocation related
N1=N+1;
LmN=L*(N+1);
ic=sqrt(-1);

Amp=1;
Omega=2;%(科氏力归一化参数)
Le=1e-2;
Ekman=1e-6;%Ek
Pm=0.01;% Pm

Em=Ekman/Pm;%Em

Pr=1;
En=Ekman/Pr; %Eman/Pr;%En ..
N2=-2;%N**2 ..

ri=0.35;
omg=1.1;
m=2;% Azimuthal wavenumber
EqS=0; % 0-- Symmetric 1--AntiSymmetric around the equator

MagF=0; % 0-- dipolar magnetic field 1-- Uniform vertical magnetic field  2-- Malkus Toroidal field
BCuo=1;BCui=1; %0--no-slip 1--stree free
BCbo=0;BCbi=0; %0-- Insulating inner core 1-- conducting inner core



F_flag = 2;% 0--forced problem for  MC, 1 -- forced problem for MAC , 2 -- eignen  problem for MAC
f_scan=0;
omgg=linspace(-2,0,101);
%% Define size of AA,BB
total_index = find_index(L,m,EqS);




%% Gauss-Lobatto colcation points
x=cos(pi*(0:N)/N)'; %[-1,1]
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
dxr=2/(1-ri);

%%% N2 profile
N2r=r; % 锘緿intrans et al., 1999 JFM
%N2r=(1-0.5*(1-tanh((r-0.96)/0.005)));  %锘縑idal & Schaeffer, 2015 GJI
%N2r=(r-0.96)/(1-0.96)+abs((r-0.96)/(1-0.96)); %锘縑idal & Schaeffer, 2015 GJI
%N2r=(1-0.5*(1-tanh((0.5-r)/0.001))); %%

Coef=sparse(total_index*N1,1);
Force=sparse(total_index*N1,1);
%%%% Buffet 2010 Nature  u_2^1=Amp on the outer boundary
tic
if MagF==2
    [AA,BB]=GenMat_Mulkus(ri,m,L,N,Omega,Ekman,Le,Pm,EqS,BCuo,BCui,BCbo,BCbi);
else
    [AA,BB]=GenMat(ri,m,L,N,Omega,Ekman,Le,N2,En,Em,EqS,MagF,BCuo,BCui,BCbi,N2r);
end
toc

%%%%%%%%%% add the constraint to remove the spin-over mode
if F_flag == 0
    AA=ic*omg*AA+BB;
    Force(3*N1)=Amp*2*(1-0.01)*0.0025*2*sqrt(2*pi/15)*ri; % Correction one
    spparms('spumoni',1);
    Coef=-AA\Force;
    res=full([real(Coef),imag(Coef)]);
elseif F_flag == 1
    m=2;
    EqS=0;
    Force(1)=Amp;
    AA1=ic*omg*AA+BB;
    Coef=-AA1\Force;
    res=full([real(Coef),imag(Coef)]);
elseif F_flag == 2
    tic
    display('Solve the eigen-value problem (matlab solver)');
    lmd0=-ic*omg;
    if m==0
        if EqS==0
            [EigVec,EigVal]=eigs(BB(2*N1+1:1.5*LmN,2*N1+1:1.5*LmN),AA(2*N1+1:1.5*LmN,2*N1+1:1.5*LmN),1,lmd0);
            Coef(2*N1+1:1.5*LmN)=EigVec;
        else
            [EigVec,EigVal]=eigs(BB(N1+1:1.5*LmN,N1+1:1.5*LmN),AA(N1+1:1.5*LmN,N1+1:1.5*LmN),1,lmd0);
            Coef(N1+1:1.5*LmN)=EigVec;
        end
        EigVal
    else
        if f_scan==0
            [EigVec,EigVal]=eigs(BB,AA,1,lmd0);
            Coef=EigVec;
        else
            
            for nn=1:length(omgg)
                lmd0=-ic*omgg(nn);
                [EigVec,EigVal]=eigs(BB,AA,1,lmd0);
                Coef=EigVec;
                scanf(nn,:)=[omgg(nn),real(EigVal),imag(EigVal)]
            end
            save('Scanf.dat','scanf','-ascii');
        end
    end
    toc
    GrowthRate=real(EigVal)
    Frequency=imag(EigVal)
    res=full([real(Coef),imag(Coef)]);
else
    warning('Please check the forcing')
end

%filename=['CoefF3B',num2str(MagF),'f',num2str(omg),'r',num2str(ri),'Ek',num2str(Ekman,'%.1e'),...
%    'Le',num2str(Le,'%.1e'),'Em',num2str(Em,'%.1e'),'L',num2str(L),'N',num2str(N),'.dat'];
%save(filename,'res','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Extract coeffecients
CoefMat=reshape(Coef,N1,2.5*L);

if MagF==2
    unl=CoefMat(:,2*EqS+1:4:2*L);
    wnl=CoefMat(:,3-2*EqS:4:2*L);
    cnl=CoefMat(:,4-2*EqS:4:2*L);
    anl=CoefMat(:,2*EqS+2:4:2*L);
else
    unl=CoefMat(:,2*EqS+1:4:2*L);
    wnl=CoefMat(:,3-2*EqS:4:2*L);
    anl=CoefMat(:,4-2*EqS:4:2*L);
    cnl=CoefMat(:,2*EqS+2:4:2*L);
end


%%%%% Integrals of the total energy and dissipation and power input
%% based on Rieutord coefficients
[KEU,Diss]=Int_KE_Diss(ri,m,N,L,EqS,unl,bsxfun(@times,wnl,1./r));
DissU=Diss*Ekman;

if MagF==2
    [KEBR,DissR]=Int_KE_Diss_Rieutord(ri,m,N,L,EqS,anl,cnl);
else
    [KEBR,DissR]=Int_KE_Diss_Rieutord(ri,m,N,L,1-EqS,anl,cnl);
end

%%% Calculate integral quantities

KEBR=KEBR*Le^2;
DissBR=DissR*Em*Le^2
DissTot=DissU+DissBR

res=[Ekman,Le,Em,KEU,DissU,KEBR,DissBR,L,N];


%% Spectral coefficients to Physical space
theta=linspace(0.0001,pi-0.0001,L+1);
[Vel,Vor,Helicity]=Coef2Phy(ri,m,N,L,EqS,unl,wnl,r,theta);
[MagB,CurrJ,HelB]=Coef2Phy(ri,m,N,L,1-EqS,anl,cnl,r,theta);


%Mode=['m',num2str(m),'Fre',num2str(-imag(EigVal),'%.2e')]

%%%% Equatorial plane
phi=linspace(0,2*pi,L);
Temp_eq_r=Vel(:,L/2+1,1);
Temp_eq=real(Temp_eq_r)*cos(m*phi)-imag(Temp_eq_r)*sin(m*phi);
xx=1.0*r*cos(phi);
yy=1.0*r*sin(phi);
figure
hold on
axis equal
axis off
pcolor(xx,yy,Temp_eq); shading interp;
colormap('jet')
xb=cos(linspace(0, 2*pi, 100));
yb=sin(linspace(0, 2*pi, 100));
plot(xb,yb,'k-','LineWidth',1.0);

%fileformat = "D:/MAC_test/omg=%.4f/%s.png";
%ddle = Le;
%ddEkman = Ekman;
%dd = Frequency;
%ss1 = 'EquatorialPlane_0';
%str = sprintf(fileformat,dd,ss1);
%saveas(gcf,str)

%%% CMB
Data_r0=MagB(1,:,1);
Data_CMB=real(Data_r0')*cos(m*phi)-imag(Data_r0')*sin(m*phi);
lat=90-theta*180/pi;
lon=phi*180/pi;
figure('Position',[100 100 700 400],'color','w')
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on');
%axesm ('hammer', 'Frame', 'on', 'Grid', 'on');
pcolorm(lat,lon,Data_CMB);
%colormap(mycolor);
axis off
hold on
%h=colorbar;
%title(['$\omega$=' num2str(-imag(EigVal),'%.2e')],'Interpreter','latex','FontSize',20)
%export_fig([Mode,'CMB.png'],'-m 2')


[tt,rr]=meshgrid(theta,r);
x=rr.*sin(tt);
y=rr.*cos(tt);

%ss1 = 'CMB_0';
%str = sprintf(fileformat,dd,ss1);
%saveas(gcf,str)


figure('position',[100 100 400 400],'Color','w');
subplot('Position',[0.1,0.02,0.4,0.95])
axis off
axis equal
hold on
pcolor(-x,y,real(Vel(:,:,3))); shading interp;
%colorbar('WestOutside')
ax1=gca;
colormap(ax1,'jet');
%plot([-ri, -1],[0,0],'k-');
plot(-ri*sin(theta),ri*cos(theta),'k-');
plot(-sin(theta),cos(theta),'k-');
%xlim([-1,0])
subplot('Position',[0.5,0.02,0.4,0.95])
axis off
axis equal
hold on
h2=pcolor(x,y,real(MagB(:,:,1))); shading interp;
%colorbar('EastOutside')
ax2=gca;
colormap(ax2);
%colormap(ax2,mycolor);
hold on
plot([0, 0],[-ri,-1],'k-');
plot([0, 0],[ri,1],'k-');
plot(sin(theta),cos(theta),'k-');
plot(ri*sin(theta),ri*cos(theta),'k-');
%ss1 = 'merid_0';
%str = sprintf(fileformat,dd,ss1);
%saveas(gcf,str)
%export_fig([Mode,'merid.png'],'-m 2')
% figure('position',[100 100 300 600],'Color','w');
% axis off
% axis equal
% hold on
% pcolor(-x,y,real(Vel(:,:,3))); shading interp;
% %colorbar('WestOutside')
% ax1=gca;
% colormap(mycolor);
% export_fig Vel.png -transparent

% figure('position',[100 100 300 600],'Color','w');
% axis off
% axis equal
% hold on
% pcolor(x,y,real(MagB(:,:,1))); shading interp;
% %colorbar('WestOutside')
% ax1=gca;
% colormap(ax1,'jet');
%export_fig Mag.png -transparent