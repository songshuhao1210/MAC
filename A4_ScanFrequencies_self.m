Le = 1e-2;
Ekman = 1e-2;
N2 = -2;
N2r_mode = 1;

F_flag = 2;% 0--forced problem for  MC, 1 -- forced problem for MAC , 2 -- eignen  problem for MAC
Plot_flag = 0;%0--plot 3 single plots ; 1--plot together
new_folder = ['D:/MAC/','Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'_mode=',num2str(N2r_mode,'%d'),'/'];
mkdir(new_folder);
mkdir([new_folder,'\figure']);
mkdir([new_folder,'\para']);

ScanFrequency(Le,Ekman,N2,F_flag,Plot_flag,new_folder);