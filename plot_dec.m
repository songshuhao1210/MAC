function[] = plot_dec(m,Le,Ekman,N2,N2r_mode,path)
%% read and plot
%path = ['D:/MAC/','Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'/','_mode=',num2str(N2r_modeS,'%d'),'/'];
xlfileformat = path;
%xlfileformat =  ['D:\MAC\','Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'\'];

freqs = xlsread([xlfileformat,'datarecord.xlsx'],'A:A');
dec =abs( xlsread([xlfileformat,'datarecord.xlsx'],'B:B'));
flag_this = 0; %0 -- normal; 1--too little values; 2--first too large;

max_index = find(dec == min(dec));

if length(freqs) <= 4
    flag_this = 1;
elseif abs(1/min(dec)) > abs(40*1/dec(max_index+1))
    plot_dec_1(m,Le,Ekman,N2,N2r_mode,path);
    flag_this = 2;
end

%close 

if length(dec) > 1
    y_max = max(1./dec);
    y_min = min(1./dec);
    x_max = max(freqs);
    x_min = min(freqs);
    
    ylim([y_min,y_max])
    xlim([x_min,x_max])
end

%hold on

%line([min(freqs),max(freqs)],[1/mean_dec,1/mean_dec]);

%% text
if flag_this == 0
    grid on
    dec_trans = 1./dec;
    flag = 0.5*(max(dec_trans) - min(dec_trans));
    
    stem(freqs(1./dec>= flag),1./dec(1./dec>= flag),'red-o','MarkerFaceColor','red','MarkerEdgeColor','green');
    hold on
    stem(freqs(1./dec< flag),1./dec(1./dec< flag),'blue-.o');
    hold on
    line([0,0],[y_min,y_max],'color','green','Linestyle','--')
    
    f_text = freqs(dec_trans>flag);
    dec_text = dec_trans(dec_trans>flag);
    
    for i = 1:length(f_text)
        if length(f_text) >= 10
            text(f_text(i),dec_text(i),['(',num2str(roundn(f_text(i),-4)),',',num2str(int32(dec_text(i))),')'],'color','black','FontSize',7)
        else
            text(f_text(i),dec_text(i),['(',num2str(roundn(f_text(i),-4)),',',num2str(int32(dec_text(i))),')'],'color','black')
        end
    end
    
else
    %close Figure 1
    grid off
    stem(freqs,1./dec,'black-o');
end

%% labels
title(['    Decay rate; ','m=',num2str(m,'%d'),' Le=',num2str(Le,'%.1e'),'  Ekman=',num2str(Ekman,'%.1e'),' N2=',num2str(N2,'%d'),' mode=',num2str(N2r_mode,'%d')])
xlabel('freqency')
ylabel('1/decay rate')

%% saves
pngfileformat = "D:/MAC/m=%d_Le=%.1e_Ekman=%.1e_N2=%d_mode=%d/TotalDecayRate_comparison.png";
str = sprintf(pngfileformat,m,Le,Ekman,N2,N2r_mode);
pause(0.1)
saveas(gcf,str)



end