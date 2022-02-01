Le = 0.01;
N2 = -2;
Ekman = 1e-6;
new_folder = ['D:\MAC\','Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d')];
mkdir(new_folder);
mkdir([new_folder,'\figure']);
mkdir([new_folder,'\para']);

titlex = [{'Frequency'},{'GrowthRate'},{'DissTot'},{'DissBR'},{'Omg'}];

xlfileformat = ['D:\MAC\','Le=',num2str(Le,'%.1e'),'_Ekman=',num2str(Ekman,'%.1e'),'_N2=',num2str(N2,'%d'),'\'];
xlswrite([xlfileformat,'datarecord.xlsx'],titlex,'sheet1','A1');
strFormat = "A%d";

id = 2;
strId = sprintf( strFormat,id);
dataset = [Frequency,GrowthRate,DissBR,DissTot,omg];
xlswrite([xlfileformat,'datarecord.xlsx'],dataset,'sheet1',strId);