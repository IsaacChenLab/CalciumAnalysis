function cat(catname,datee)

global CAT;
global DATE;
CAT=catname;
if nargin>1
    DATE=datee;
else
    DATE='';
    disp('note, date was not assigned');
end


 
 
 

