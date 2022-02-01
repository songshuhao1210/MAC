LeS = [1e-4];
EkmanS = [1e-6];
N2S = -4:-1;
MS = 2:7;
N2r_modeS = [3];%1--N2r=r; 2--N2r=(1-0.5*(1-tanh((r-0.96)/0.005))); 3--N2r=(r-0.96)/(1-0.96)+abs((r-0.96)/(1-0.96)); 4--N2r=(1-0.5*(1-tanh((0.5-r)/0.001)));

L=200 % resolution in latitude, Should be even
N=100%N>4 Chebyshev polymonical numbers£¬collocation related
ri=0.35;%ÆðÊ¼°ë¾¶

omg0 = -0.1;
omg_max = 0.1;
step0 = 0.001;
iteration_flag = 0;%0--no self adaptive; 1--self adaptive

test_flag = 0;

F_flag = 2;% 0--forced problem for  MC, 1 -- forced problem for MAC , 2 -- eignen  problem for MAC


Plot_flag = 0;%0--plot 3 single plots ; 1--plot together

titlex = [{'m'},{'Le'},{'Ekman'},{'N2'},{'N2r_mode'}];
path0 = ['D:/MAC/'];
xlswrite([path0,'error_data.xlsx'],titlex,'sheet1','A1');
id_0 = 2;

for m = MS(1):MS(end)
    for i = 1:length(LeS)
        for j = 1:length(EkmanS)
            for k = 1:length(N2S)
                for l = 1:length(N2r_modeS)
                    %if i==1 && j==1 && k==1 && l==1 && m==2
                    %    continue;
                    %end
                    
                    if test_flag == 0
                        Le = LeS(i);
                        Ekman = EkmanS(j);
                        N2 = N2S(k);
                        N2r_mode = N2r_modeS(l);
                        
                        plot_mode(N2r_mode,N,ri,path0);
                        
                        path = ['D:/MAC/','m=',num2str(m,'%d'),'_Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'_mode=',num2str(N2r_mode,'%d'),'/'];
                        new_folder = path;
                        mkdir(new_folder);
                        mkdir([new_folder,'/figure']);
                        mkdir([new_folder,'/para']);
                        ScanFrequency(LeS(i),EkmanS(j),N2S(k),F_flag,Plot_flag,path,N2r_mode,m,omg0,omg_max,step0,iteration_flag,L,N,ri);
                    else
                    
                        try
                            Le = LeS(i);
                            Ekman = EkmanS(j);
                            N2 = N2S(k);
                            N2r_mode = N2r_modeS(l);
                            
                            path = ['D:/MAC/','m=',num2str(m,'%d'),'_Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'_mode=',num2str(N2r_mode,'%d'),'/'];
                            new_folder = path;
                            mkdir(new_folder);
                            mkdir([new_folder,'/figure']);
                            mkdir([new_folder,'/para']);
                            ScanFrequency(LeS(i),EkmanS(j),N2S(k),F_flag,Plot_flag,path,N2r_mode,m.omg0,omg_max,iteration_flag);
                        catch
                            
                            strFormat_0 = "A%d";
                            strId_0 = sprintf(strFormat_0,id_0);
                            dataset_0 = [Le,Ekman,N2,N2r_mode,m];
                            xlswrite([path_0,'error_data.xlsx'],dataset_0,'sheet1',strId_0);
                            id_0 = id_0+1;
                            
                            plot_dec(m,Le,Ekman,N2,N2r_mode,path);
                            continue;
                        end
                        
                    end
                end
            end
        end
    end
end
