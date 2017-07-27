%Create timestamp arrays for different TTL signals.

Timestamp0=[];
Timestamp1=[];
Timestamp2=[];
Timestamp3=[];
Timestamp4=[];
Timestamp12=[];
Timestamp13=[];
Timestamp14=[];
Timestamp23=[];
Timestamp24=[];
Timestamp34=[];
Timestamp123=[];
Timestamp124=[];
Timestamp134=[];
Timestamp234=[];
Timestamp1234=[];

for b=1:numel(TTLs)
    if TTLs(b)==0
        Timestamp0=[Timestamp0 TimestampEvents(b)];
    else if TTLs(b)==1
            Timestamp2=[Timestamp2 TimestampEvents(b)];
        else if TTLs(b)==2
                Timestamp1=[Timestamp1 TimestampEvents(b)];
            else if TTLs(b)==3
                    Timestamp12=[Timestamp12 TimestampEvents(b)];
                else if TTLs(b)==4
                        Timestamp3=[Timestamp3 TimestampEvents(b)];
                    else if TTLs(b)==5
                            Timestamp23=[Timestamp23 TimestampEvents(b)];
                        else if TTLs(b)==6
                                Timestamp13=[Timestamp13 TimestampEvents(b)];
                            else if TTLs(b)==7
                                    Timestamp123=[Timestamp123 TimestampEvents(b)];
                                else if TTLs(b)==8
                                        Timestamp4=[Timestamp4 TimestampEvents(b)];
                                    else if TTLs(b)==9
                                            Timestamp24=[Timestamp24 TimestampEvents(b)];
                                        else if TTLs(b)==10
                                                Timestamp14=[Timestamp14 TimestampEvents(b)];
                                            else if TTLs(b)==11
                                                    Timestamp124=[Timestamp124 TimestampEvents(b)];
                                                else if TTLs(b)==12
                                                        Timestamp34=[Timestamp34 TimestampEvents(b)];
                                                    else if TTLs(b)==13
                                                            Timestamp234=[Timestamp234 TimestampEvents(b)];
                                                        else if TTLs(b)==14
                                                                Timestamp134=[Timestamp134 TimestampEvents(b)];
                                                            else if TTLs(b)==15
                                                                    Timestamp1234=[Timestamp1234 TimestampEvents(b)];
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Change timestamps to units of seconds.

Timestamp0=Timestamp0 / 1000000;
Timestamp1=Timestamp1 / 1000000;
Timestamp2=Timestamp2 / 1000000;
Timestamp3=Timestamp3 / 1000000;
Timestamp4=Timestamp4 / 1000000;
Timestamp12=Timestamp12 / 1000000;
Timestamp13=Timestamp13 / 1000000;
Timestamp14=Timestamp14 / 1000000;
Timestamp23=Timestamp23 / 1000000;
Timestamp24=Timestamp24 / 1000000;
Timestamp34=Timestamp34 / 1000000;
Timestamp123=Timestamp123 / 1000000;
Timestamp124=Timestamp124 / 1000000;
Timestamp134=Timestamp134 / 1000000;
Timestamp234=Timestamp234 / 1000000;
Timestamp1234=Timestamp1234 / 1000000;

%% second data set
% Timestamp12342=[];
% 
% for b=1:numel(TTLs2)
%     if TTLs2(b)==15
%         Timestamp12342=[Timestamp12342 TimestampEvents2(b)];
%     end
% end

% Timestamp12342=Timestamp12342 / 1000000;
% 
% %% third data set
% Timestamp12343=[];
% 
% for b=1:numel(TTLs3)
%     if TTLs3(b)==15
%         Timestamp12343=[Timestamp12343 TimestampEvents3(b)];
%     end
% end
% 
% Timestamp12343=Timestamp12343 / 1000000;
% 
% %% fourth data set
% Timestamp12344=[];
% 
% for b=1:numel(TTLs4)
%     if TTLs4(b)==15
%         Timestamp12344=[Timestamp12344 TimestampEvents4(b)];
%     end
% end
% 
% Timestamp12344=Timestamp12344 / 1000000;
% 
% %% fifth data set
% Timestamp12345=[];
% 
% for b=1:numel(TTLs5)
%     if TTLs5(b)==15
%         Timestamp12345=[Timestamp12345 TimestampEvents5(b)];
%     end
% end
% 
% Timestamp12345=Timestamp12345 / 1000000;