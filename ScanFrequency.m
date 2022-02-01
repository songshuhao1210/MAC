function[]=ScanFrequency(Le,Ekman,N2,F_flag,Plot_flag,path,N2r_mode,m,omg0,omg_max,step0,iteration_flag,L,N,ri)
    % Main Program
    % by Yufeng Lin 11 Nov 2016
    % based on Ricon & Rieturd,2003

    
    %% Define funtion parameters

    
    %% Define computation parameters
    
    
    Pm=0.01;% Em=Ekman/Pm
    Pr=1;%En=Ekman/Pr
    
     % Azimuthal wavenumber
    EqS=0; % 0-- Symmetric 1--AntiSymmetric around the equator

    MagF=0; % 0-- dipolar magnetic field 1-- Uniform vertical magnetic field  2-- Malkus Toroidal field
    BCuo=1;BCui=1; %0--no-slip 1--stree free
    BCbo=0;BCbi=0; %0-- Insulating inner core 1-- conducting inner core
    
    Amp=1;
    Omega=2;%(科氏力归一化参数，勿动)
    
    
    N1=N+1;
    LmN=L*(N+1);
    ic=sqrt(-1);

    Em=Ekman/Pm;
    En=Ekman/Pr; 
    
    f_scan=0;
    omgg=linspace(-2,0,101);
    %% Define size of AA,BB
    total_index = find_index(L,m,EqS);

    %% Gauss-Lobatto colcation points
    x=cos(pi*(0:N)/N)'; %[-1,1]
    r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
    dxr=2/(1-ri);

    %%% N2 profile
    if N2r_mode == 1
        N2r = r; % 锘緿intrans et al., 1999 JFM
    elseif N2r_mode == 2
        N2r = (1-0.5*(1-tanh((r-0.96)/0.005)));
    elseif N2r_mode == 3
        N2r = (r-0.96)/(1-0.96)+abs((r-0.96)/(1-0.96));
    else
        N2r = (1-0.5*(1-tanh((0.5-r)/0.001)));
    end


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


     %% scaning different omg
    omg = omg0;
    count_same = 0;
    count_loop = 0;
    step = step0;
    AA0 = AA;
    BB0 = BB;
    
    while omg <= omg_max 
        AA = AA0;
        BB = BB0;
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
            AA1=ic*omg*AA0+BB0;
            Coef=-AA1\Force;
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
            %   'Le',num2str(Le,'%.1e'),'Em',num2str(Em,'%.1e'),'L',num2str(L),'N',num2str(N),'.dat'];
            %save(filename,'res','-ascii');

        if omg == omg0
            freqs = zeros(1,2);
            freqs(1,1) = omg;
            freqs(1,2) = Frequency;
        else
            freqs = [freqs;omg,Frequency];
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Extract coeffecients
        CoefMat=reshape(Coef,N1,total_index);

        if MagF==2
            unl=CoefMat(:,2*EqS+1:5:total_index);
            wnl=CoefMat(:,3-2*EqS:5:total_index);
            cnl=CoefMat(:,4-2*EqS:5:total_index);
            anl=CoefMat(:,2*EqS+2:5:total_index);
        else
            unl=CoefMat(:,2*EqS+1:4:2*L);
            wnl=CoefMat(:,3-2*EqS:4:2*L);
            anl=CoefMat(:,4-2*EqS:4:2*L);
            cnl=CoefMat(:,2*EqS+2:4:2*L);
            tnl=CoefMat(:,2*L+1:total_index);
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
        DissBR=DissR*Em*Le^2;
        DissTot=DissU+DissBR;

        res=[Ekman,Le,Em,KEU,DissU,KEBR,DissBR,L,N];

        if omg == omg0
            %plot 
            myPlot(Plot_flag,Frequency,ri,m,N,L,EqS,unl,wnl,anl,cnl,r,path)
            %record dataset into datarecord.xlsx
            titlex = [{'Frequency'},{'GrowthRate'},{'DissTot'},{'DissBR'},{'Omg'},{'Times'}];
            xlfileformat = path;
            xlswrite([xlfileformat,'datarecord.xlsx'],titlex,'sheet1','A1');
            strFormat = "A%d";

            id = 2;
            strId = sprintf( strFormat,id);
            dataset = [Frequency,GrowthRate,DissBR,DissTot,omg];
            xlswrite([xlfileformat,'datarecord.xlsx'],dataset,'sheet1',strId);
            %record dataset into para with single files
            filename=['CoefF3B',num2str(MagF),'f',num2str(omg),'r',num2str(ri),'Ek',num2str(Ekman,'%.1e'),...
                'Le',num2str(Le,'%.1e'),'Em',num2str(Em,'%.1e'),'L',num2str(L),'N',num2str(N),'.dat'];
            filepath = [path,'\para\'];
            save([filepath, filename],'res','-ascii');
            
        else
            if abs(abs(freqs(count_loop+1,2)) -abs( freqs(count_loop,2))) > 0.0001
                myPlot(Plot_flag,Frequency,ri,m,N,L,EqS,unl,wnl,anl,cnl,r,path)
                
                
                strFormat_2 = "F%d";
                strID_2 = sprintf(strFormat_2,id);
                xlswrite([xlfileformat,'datarecord.xlsx'],count_same+1,'sheet1',strID_2);
                
                id = id+1;
                strId = sprintf( strFormat,id);
                dataset = [Frequency,GrowthRate,DissBR,DissTot,omg];

                xlfileformat = path;
                xlswrite([xlfileformat,'datarecord.xlsx'],dataset,'sheet1',strId);

                filename=['CoefF3B',num2str(MagF),'f',num2str(omg),'r',num2str(ri),'Ek',num2str(Ekman,'%.1e'),...
                    'Le',num2str(Le,'%.1e'),'Em',num2str(Em,'%.1e'),'L',num2str(L),'N',num2str(N),'.dat'];
                filepath = [path,'\para\'];
                save([filepath, filename],'res','-ascii');

                count_same = 0
            else
                count_same = count_same + 1

            end
        end

        close all
        
        if iteration_flag == 1
            if count_same == 5
                step = step+step0;
            elseif count_same == 10
                step = step+5*step0;
            elseif count_same == 20
                step = step*5;
            elseif count_same == 40
                step = step*10;
            elseif abs(abs(omg) - abs(Frequency)) > abs(30*Frequency)
                step = abs(abs(omg) - abs(Frequency))
                omg = omg+step;
                xlswrite([xlfileformat,'datarecord.xlsx'],dataset,'sheet2',strId)
            end       
        end

        

        step
        omg = omg+step

        count_loop = count_loop + 1

    end
    plot_dec(m,Le,Ekman,N2,N2r_mode,path);
    clear all
    clear all
end