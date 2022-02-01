N = 100;
ri = 0.35;
x=cos(pi*(0:N)/N)'; %[-1,1]
r=(1-ri)/2*x+(1+ri)/2; % x=2*r/(1-ri)-(1+ri)/(1-ri);

N2r_modeS = zeros(size(r,1),4);

N2r_modeS(:,1) = r;
N2r_modeS(:,2) = (1-0.5*(1-tanh((r-0.96)/0.005))); 
N2r_modeS(:,3) = (r-0.96)/(1-0.96)+abs((r-0.96)/(1-0.96)); 
N2r_modeS(:,4)=(1-0.5*(1-tanh((0.5-r)/0.001)));

for i = 1:size(N2r_modeS,2)
    figure
    plot(r,N2r_modeS(:,i)');
    xlabel('r');
    ylabel('N2r');
    title(['mode',num2str(i)]);
end