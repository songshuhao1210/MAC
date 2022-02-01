%Le = 1e-2;
%Ekman = 1e-6;
%N2 = -3;
%N2r_mode = 3 ;
%m = 2;
%path = ['D:/MAC/','m=',num2str(m,'%d'),'_Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'_mode=',num2str(N2r_mode,'%d'),'/'];

%plot_flag_1 = 1;% 1--single plot, 2--reduce first value , 3--different plot

LeS = [1e-4];
EkmanS = [1e-6];
N2S = -4:-1;
MS = 5:7;
N2r_modeS = [3];
plot_flag_1 = 3;%plot_flag_1 = 1;% 1--single plot , 2--different plots 




if plot_flag_1 == 1
    plot_dec(m,Le,Ekman,N2,N2r_mode,path);
elseif plot_flag_1 == 2
    plot_dec_1(m,Le,Ekman,N2,N2r_mode,path);
else
    for m = MS(1):MS(end)
        for i = 1:length(LeS)
            for j = 1:length(EkmanS)
                for k = 1:length(N2S)
                    %if i==1 && j==1 && k==1
                    %    continue;
                    %end
                    Le = LeS(i);
                    Ekman = EkmanS(j);
                    N2 = N2S(k);
                    path = ['D:/MAC/','m=',num2str(m,'%d'),'_Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'_mode=',num2str(N2r_mode,'%d'),'/'];
                    plot_dec_1(m,Le,Ekman,N2,N2r_mode,path);
                end
            end
        end
    end
    
    
end