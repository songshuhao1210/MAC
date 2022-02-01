function[]=plot_mode(N2r_mode,N,ri,path0)
    x=cos(pi*(0:N)/N)'; %[-1,1]
    r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);
    N2r_modeS = zeros(size(r,1),1);
    
    if N2r_mode == 1
        N2r_modeS = r; 
    elseif N2r_mode == 2
        N2r_modeS = (1-0.5*(1-tanh((r-0.96)/0.005)));
    elseif N2r_mode == 3
        N2r_modeS = (r-0.96)/(1-0.96)+abs((r-0.96)/(1-0.96));
    else
        N2r_modeS = (1-0.5*(1-tanh((0.5-r)/0.001)));
    end
    
    figure
    plot(r,N2r_modeS');
    xlabel('r');
    ylabel('N2r');
    title(['mode',num2str(N2r_mode)]);
    
    fileformat = strcat(path0,"mode %d.png");
    str = sprintf(fileformat,N2r_mode);
    pause(0.1)
    saveas(gcf,str)

end