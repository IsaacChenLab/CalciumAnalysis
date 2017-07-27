function [dispData,xt,Hz]=MC_spGetElectrodesDisplay(file,ch)

global CAT;
global DATE;
global DATA_DIR;

cd(DATA_DIR);
cd('DisplayFiles');

s=sprintf('load D%s_%s_%d-%d dispData_%d xt Hz;',CAT,DATE,file(1),file(end),ch);
eval(s);
s=sprintf('dispData=dispData_%d;',ch);
eval(s);
clear dispData_*;

return;


