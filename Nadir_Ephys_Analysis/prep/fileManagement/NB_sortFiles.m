% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

function [] = NB_sortFiles(numRecordings)

for rec = 0:numRecordings-1;
    
    
    if rec == 0;
        if exist('Events.nev')
            movefile('Events.nev','./Recording_0000/');
        end
        
        for i = 1:32;
            if exist(strcat('CSC',num2str(i),'.ncs'))
                movefile(strcat('CSC',num2str(i),'.ncs'),'./Recording_0000/');
            end
        end
        for i = 1:8;
            if exist(strcat('TT',num2str(i),'.ntt'))
                movefile(strcat('TT',num2str(i),'.ntt'),'./Recording_0000/');
            end
        end
    
    elseif rec < 10;
        if exist(strcat('Events_000',num2str(rec),'.nev'))
            movefile(strcat('Events_000',num2str(rec),'.nev'),strcat('./Recording_000',num2str(rec),'/'));
        end
        for i = 1:32;
            if exist(strcat('CSC',num2str(i),'_000',num2str(rec),'.ncs'))
                movefile(strcat('CSC',num2str(i),'_000',num2str(rec),'.ncs'),strcat('./Recording_000',num2str(rec),'/'));
            end
        end
        for i = 1:8;
            if exist(strcat('TT',num2str(i),'_000',num2str(rec),'.ntt'))
                movefile(strcat('TT',num2str(i),'_000',num2str(rec),'.ntt'),strcat('./Recording_000',num2str(rec),'/'));
            end
        end
 
    elseif rec < 100;
        if exist(strcat('Events_00',num2str(rec),'.nev'))
            movefile(strcat('Events_00',num2str(rec),'.nev'),strcat('./Recording_00',num2str(rec),'/'));
        end
        
        for i = 1:32;
            if exist(strcat('CSC',num2str(i),'_00',num2str(rec),'.ncs')) 
                movefile(strcat('CSC',num2str(i),'_00',num2str(rec),'.ncs'),strcat('./Recording_00',num2str(rec),'/'));
            end           
        end
        for i = 1:8;
            if exist(strcat('TT',num2str(i),'_00',num2str(rec),'.ntt'))
                movefile(strcat('TT',num2str(i),'_00',num2str(rec),'.ntt'),strcat('./Recording_00',num2str(rec),'/'));
            end           
        end
    
    else
        if exist(strcat('Events_0',num2str(rec),'.nev'))
            movefile(strcat('Events_0',num2str(rec),'.nev'),strcat('./Recording_0',num2str(rec),'/'));
        end

        
        for i = 1:32;
            if exist(strcat('CSC',num2str(i),'_0',num2str(rec),'.ncs'))
               movefile(strcat('CSC',num2str(i),'_0',num2str(rec),'.ncs'),strcat('./Recording_0',num2str(rec),'/'));
            end            
        end
        for i = 1:8;
            if exist(strcat('TT',num2str(i),'_0',num2str(rec),'.ntt'))
              movefile(strcat('TT',num2str(i),'_0',num2str(rec),'.ntt'),strcat('./Recording_0',num2str(rec),'/'));
            end                        
        end        
    end
end