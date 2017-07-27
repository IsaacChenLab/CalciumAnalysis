%Create timestamp arrays for different cells within tetrode file.

% TT1
a=size(CellNumbersTT1);

Cell_1_0=[];
Cell_1_1=[];
Cell_1_2=[];
Cell_1_3=[];
Cell_1_4=[];
Cell_1_5=[];
Cell_1_6=[];
Cell_1_7=[];
Cell_1_8=[];
Cell_1_9=[];
Cell_1_10=[];

for b=1:a(2)
    if CellNumbersTT1(b)==0
        Cell_1_0=[Cell_1_0 TimestampsTT1(b)];
    else if CellNumbersTT1(b)==1
            Cell_1_1=[Cell_1_1 TimestampsTT1(b)];
        else if CellNumbersTT1(b)==2
                Cell_1_2=[Cell_1_2 TimestampsTT1(b)];
            else if CellNumbersTT1(b)==3
                    Cell_1_3=[Cell_1_3 TimestampsTT1(b)];
                else if CellNumbersTT1(b)==4
                        Cell_1_4=[Cell_1_4 TimestampsTT1(b)];
                    else if CellNumbersTT1(b)==5
                            Cell_1_5=[Cell_1_5 TimestampsTT1(b)];
                        else if CellNumbersTT1(b)==6
                                Cell_1_6=[Cell_1_6 TimestampsTT1(b)];
                            else if CellNumbersTT1(b)==7
                                    Cell_1_7=[Cell_1_7 TimestampsTT1(b)];
                                else if CellNumbersTT1(b)==8
                                        Cell_1_8=[Cell_1_8 TimestampsTT1(b)];
                                    else if CellNumbersTT1(b)==9
                                            Cell_1_9=[Cell_1_9 TimestampsTT1(b)];
                                        else if CellNumbersTT1(b)==10
                                                Cell_1_10=[Cell_1_10 TimestampsTT1(b)];
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

% Changes timestamps to units of seconds.

Cell_1_0=Cell_1_0 / 1000000;
Cell_1_1=Cell_1_1 / 1000000;
Cell_1_2=Cell_1_2 / 1000000;
Cell_1_3=Cell_1_3 / 1000000;
Cell_1_4=Cell_1_4 / 1000000;
Cell_1_5=Cell_1_5 / 1000000;
Cell_1_6=Cell_1_6 / 1000000;
Cell_1_7=Cell_1_7 / 1000000;
Cell_1_8=Cell_1_8 / 1000000;
Cell_1_9=Cell_1_9 / 1000000;
Cell_1_10=Cell_1_10 / 1000000;

% % TT2
% a=size(CellNumbersTT2);
% 
% Cell_2_0=[];
% Cell_2_1=[];
% Cell_2_2=[];
% Cell_2_3=[];
% Cell_2_4=[];
% Cell_2_5=[];
% Cell_2_6=[];
% Cell_2_7=[];
% Cell_2_8=[];
% Cell_2_9=[];
% Cell_2_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT2(b)==0
%         Cell_2_0=[Cell_2_0 TimestampsTT2(b)];
%     else if CellNumbersTT2(b)==1
%             Cell_2_1=[Cell_2_1 TimestampsTT2(b)];
%         else if CellNumbersTT2(b)==2
%                 Cell_2_2=[Cell_2_2 TimestampsTT2(b)];
%             else if CellNumbersTT2(b)==3
%                     Cell_2_3=[Cell_2_3 TimestampsTT2(b)];
%                 else if CellNumbersTT2(b)==4
%                         Cell_2_4=[Cell_2_4 TimestampsTT2(b)];
%                     else if CellNumbersTT2(b)==5
%                             Cell_2_5=[Cell_2_5 TimestampsTT2(b)];
%                         else if CellNumbersTT2(b)==6
%                                 Cell_2_6=[Cell_2_6 TimestampsTT2(b)];
%                             else if CellNumbersTT2(b)==7
%                                     Cell_2_7=[Cell_2_7 TimestampsTT2(b)];
%                                 else if CellNumbersTT2(b)==8
%                                         Cell_2_8=[Cell_2_8 TimestampsTT2(b)];
%                                     else if CellNumbersTT2(b)==9
%                                             Cell_2_9=[Cell_2_9 TimestampsTT2(b)];
%                                         else if CellNumbersTT2(b)==10
%                                                 Cell_2_10=[Cell_2_10 TimestampsTT2(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_2_0=Cell_2_0 / 1000000;
% Cell_2_1=Cell_2_1 / 1000000;
% Cell_2_2=Cell_2_2 / 1000000;
% Cell_2_3=Cell_2_3 / 1000000;
% Cell_2_4=Cell_2_4 / 1000000;
% Cell_2_5=Cell_2_5 / 1000000;
% Cell_2_6=Cell_2_6 / 1000000;
% Cell_2_7=Cell_2_7 / 1000000;
% Cell_2_8=Cell_2_8 / 1000000;
% Cell_2_9=Cell_2_9 / 1000000;
% Cell_2_10=Cell_2_10 / 1000000;
% 
% %TT3
% a=size(CellNumbersTT3);
% 
% Cell_3_0=[];
% Cell_3_1=[];
% Cell_3_2=[];
% Cell_3_3=[];
% Cell_3_4=[];
% Cell_3_5=[];
% Cell_3_6=[];
% Cell_3_7=[];
% Cell_3_8=[];
% Cell_3_9=[];
% Cell_3_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT3(b)==0
%         Cell_3_0=[Cell_3_0 TimestampsTT3(b)];
%     else if CellNumbersTT3(b)==1
%             Cell_3_1=[Cell_3_1 TimestampsTT3(b)];
%         else if CellNumbersTT3(b)==2
%                 Cell_3_2=[Cell_3_2 TimestampsTT3(b)];
%             else if CellNumbersTT3(b)==3
%                     Cell_3_3=[Cell_3_3 TimestampsTT3(b)];
%                 else if CellNumbersTT3(b)==4
%                         Cell_3_4=[Cell_3_4 TimestampsTT3(b)];
%                     else if CellNumbersTT3(b)==5
%                             Cell_3_5=[Cell_3_5 TimestampsTT3(b)];
%                         else if CellNumbersTT3(b)==6
%                                 Cell_3_6=[Cell_3_6 TimestampsTT3(b)];
%                             else if CellNumbersTT3(b)==7
%                                     Cell_3_7=[Cell_3_7 TimestampsTT3(b)];
%                                 else if CellNumbersTT3(b)==8
%                                         Cell_3_8=[Cell_3_8 TimestampsTT3(b)];
%                                     else if CellNumbersTT3(b)==9
%                                             Cell_3_9=[Cell_3_9 TimestampsTT3(b)];
%                                         else if CellNumbersTT3(b)==10
%                                                 Cell_3_10=[Cell_3_10 TimestampsTT3(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_3_0=Cell_3_0 / 1000000;
% Cell_3_1=Cell_3_1 / 1000000;
% Cell_3_2=Cell_3_2 / 1000000;
% Cell_3_3=Cell_3_3 / 1000000;
% Cell_3_4=Cell_3_4 / 1000000;
% Cell_3_5=Cell_3_5 / 1000000;
% Cell_3_6=Cell_3_6 / 1000000;
% Cell_3_7=Cell_3_7 / 1000000;
% Cell_3_8=Cell_3_8 / 1000000;
% Cell_3_9=Cell_3_9 / 1000000;
% Cell_3_10=Cell_3_10 / 1000000;
% 
% % TT4
% a=size(CellNumbersTT4);
% 
% Cell_4_0=[];
% Cell_4_1=[];
% Cell_4_2=[];
% Cell_4_3=[];
% Cell_4_4=[];
% Cell_4_5=[];
% Cell_4_6=[];
% Cell_4_7=[];
% Cell_4_8=[];
% Cell_4_9=[];
% Cell_4_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT4(b)==0
%         Cell_4_0=[Cell_4_0 TimestampsTT4(b)];
%     else if CellNumbersTT4(b)==1
%             Cell_4_1=[Cell_4_1 TimestampsTT4(b)];
%         else if CellNumbersTT4(b)==2
%                 Cell_4_2=[Cell_4_2 TimestampsTT4(b)];
%             else if CellNumbersTT4(b)==3
%                     Cell_4_3=[Cell_4_3 TimestampsTT4(b)];
%                 else if CellNumbersTT4(b)==4
%                         Cell_4_4=[Cell_4_4 TimestampsTT4(b)];
%                     else if CellNumbersTT4(b)==5
%                             Cell_4_5=[Cell_4_5 TimestampsTT4(b)];
%                         else if CellNumbersTT4(b)==6
%                                 Cell_4_6=[Cell_4_6 TimestampsTT4(b)];
%                             else if CellNumbersTT4(b)==7
%                                     Cell_4_7=[Cell_4_7 TimestampsTT4(b)];
%                                 else if CellNumbersTT4(b)==8
%                                         Cell_4_8=[Cell_4_8 TimestampsTT4(b)];
%                                     else if CellNumbersTT4(b)==9
%                                             Cell_4_9=[Cell_4_9 TimestampsTT4(b)];
%                                         else if CellNumbersTT4(b)==10
%                                                 Cell_4_10=[Cell_4_10 TimestampsTT4(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_4_0=Cell_4_0 / 1000000;
% Cell_4_1=Cell_4_1 / 1000000;
% Cell_4_2=Cell_4_2 / 1000000;
% Cell_4_3=Cell_4_3 / 1000000;
% Cell_4_4=Cell_4_4 / 1000000;
% Cell_4_5=Cell_4_5 / 1000000;
% Cell_4_6=Cell_4_6 / 1000000;
% Cell_4_7=Cell_4_7 / 1000000;
% Cell_4_8=Cell_4_8 / 1000000;
% Cell_4_9=Cell_4_9 / 1000000;
% Cell_4_10=Cell_4_10 / 1000000;
% 
% % TT5
% a=size(CellNumbersTT5);
% 
% Cell_5_0=[];
% Cell_5_1=[];
% Cell_5_2=[];
% Cell_5_3=[];
% Cell_5_4=[];
% Cell_5_5=[];
% Cell_5_6=[];
% Cell_5_7=[];
% Cell_5_8=[];
% Cell_5_9=[];
% Cell_5_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT5(b)==0
%         Cell_5_0=[Cell_5_0 TimestampsTT5(b)];
%     else if CellNumbersTT5(b)==1
%             Cell_5_1=[Cell_5_1 TimestampsTT5(b)];
%         else if CellNumbersTT5(b)==2
%                 Cell_5_2=[Cell_5_2 TimestampsTT5(b)];
%             else if CellNumbersTT5(b)==3
%                     Cell_5_3=[Cell_5_3 TimestampsTT5(b)];
%                 else if CellNumbersTT5(b)==4
%                         Cell_5_4=[Cell_5_4 TimestampsTT5(b)];
%                     else if CellNumbersTT5(b)==5
%                             Cell_5_5=[Cell_5_5 TimestampsTT5(b)];
%                         else if CellNumbersTT5(b)==6
%                                 Cell_5_6=[Cell_5_6 TimestampsTT5(b)];
%                             else if CellNumbersTT5(b)==7
%                                     Cell_5_7=[Cell_5_7 TimestampsTT5(b)];
%                                 else if CellNumbersTT5(b)==8
%                                         Cell_5_8=[Cell_5_8 TimestampsTT5(b)];
%                                     else if CellNumbersTT5(b)==9
%                                             Cell_5_9=[Cell_5_9 TimestampsTT5(b)];
%                                         else if CellNumbersTT5(b)==10
%                                                 Cell_5_10=[Cell_5_10 TimestampsTT5(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_5_0=Cell_5_0 / 1000000;
% Cell_5_1=Cell_5_1 / 1000000;
% Cell_5_2=Cell_5_2 / 1000000;
% Cell_5_3=Cell_5_3 / 1000000;
% Cell_5_4=Cell_5_4 / 1000000;
% Cell_5_5=Cell_5_5 / 1000000;
% Cell_5_6=Cell_5_6 / 1000000;
% Cell_5_7=Cell_5_7 / 1000000;
% Cell_5_8=Cell_5_8 / 1000000;
% Cell_5_9=Cell_5_9 / 1000000;
% Cell_5_10=Cell_5_10 / 1000000;
% 
% % TT6
% a=size(CellNumbersTT6);
% 
% Cell_6_0=[];
% Cell_6_1=[];
% Cell_6_2=[];
% Cell_6_3=[];
% Cell_6_4=[];
% Cell_6_5=[];
% Cell_6_6=[];
% Cell_6_7=[];
% Cell_6_8=[];
% Cell_6_9=[];
% Cell_6_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT6(b)==0
%         Cell_6_0=[Cell_6_0 TimestampsTT6(b)];
%     else if CellNumbersTT6(b)==1
%             Cell_6_1=[Cell_6_1 TimestampsTT6(b)];
%         else if CellNumbersTT6(b)==2
%                 Cell_6_2=[Cell_6_2 TimestampsTT6(b)];
%             else if CellNumbersTT6(b)==3
%                     Cell_6_3=[Cell_6_3 TimestampsTT6(b)];
%                 else if CellNumbersTT6(b)==4
%                         Cell_6_4=[Cell_6_4 TimestampsTT6(b)];
%                     else if CellNumbersTT6(b)==5
%                             Cell_6_5=[Cell_6_5 TimestampsTT6(b)];
%                         else if CellNumbersTT6(b)==6
%                                 Cell_6_6=[Cell_6_6 TimestampsTT6(b)];
%                             else if CellNumbersTT6(b)==7
%                                     Cell_6_7=[Cell_6_7 TimestampsTT6(b)];
%                                 else if CellNumbersTT6(b)==8
%                                         Cell_6_8=[Cell_6_8 TimestampsTT6(b)];
%                                     else if CellNumbersTT6(b)==9
%                                             Cell_6_9=[Cell_6_9 TimestampsTT6(b)];
%                                         else if CellNumbersTT6(b)==10
%                                                 Cell_6_10=[Cell_6_10 TimestampsTT6(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_6_0=Cell_6_0 / 1000000;
% Cell_6_1=Cell_6_1 / 1000000;
% Cell_6_2=Cell_6_2 / 1000000;
% Cell_6_3=Cell_6_3 / 1000000;
% Cell_6_4=Cell_6_4 / 1000000;
% Cell_6_5=Cell_6_5 / 1000000;
% Cell_6_6=Cell_6_6 / 1000000;
% Cell_6_7=Cell_6_7 / 1000000;
% Cell_6_8=Cell_6_8 / 1000000;
% Cell_6_9=Cell_6_9 / 1000000;
% Cell_6_10=Cell_6_10 / 1000000;
% 
% % TT7
% a=size(CellNumbersTT7);
% 
% Cell_7_0=[];
% Cell_7_1=[];
% Cell_7_2=[];
% Cell_7_3=[];
% Cell_7_4=[];
% Cell_7_5=[];
% Cell_7_6=[];
% Cell_7_7=[];
% Cell_7_8=[];
% Cell_7_9=[];
% Cell_7_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT7(b)==0
%         Cell_7_0=[Cell_7_0 TimestampsTT7(b)];
%     else if CellNumbersTT7(b)==1
%             Cell_7_1=[Cell_7_1 TimestampsTT7(b)];
%         else if CellNumbersTT7(b)==2
%                 Cell_7_2=[Cell_7_2 TimestampsTT7(b)];
%             else if CellNumbersTT7(b)==3
%                     Cell_7_3=[Cell_7_3 TimestampsTT7(b)];
%                 else if CellNumbersTT7(b)==4
%                         Cell_7_4=[Cell_7_4 TimestampsTT7(b)];
%                     else if CellNumbersTT7(b)==5
%                             Cell_7_5=[Cell_7_5 TimestampsTT7(b)];
%                         else if CellNumbersTT7(b)==6
%                                 Cell_7_6=[Cell_7_6 TimestampsTT7(b)];
%                             else if CellNumbersTT7(b)==7
%                                     Cell_7_7=[Cell_7_7 TimestampsTT7(b)];
%                                 else if CellNumbersTT7(b)==8
%                                         Cell_7_8=[Cell_7_8 TimestampsTT7(b)];
%                                     else if CellNumbersTT7(b)==9
%                                             Cell_7_9=[Cell_7_9 TimestampsTT7(b)];
%                                         else if CellNumbersTT7(b)==10
%                                                 Cell_7_10=[Cell_7_10 TimestampsTT7(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_7_0=Cell_7_0 / 1000000;
% Cell_7_1=Cell_7_1 / 1000000;
% Cell_7_2=Cell_7_2 / 1000000;
% Cell_7_3=Cell_7_3 / 1000000;
% Cell_7_4=Cell_7_4 / 1000000;
% Cell_7_5=Cell_7_5 / 1000000;
% Cell_7_6=Cell_7_6 / 1000000;
% Cell_7_7=Cell_7_7 / 1000000;
% Cell_7_8=Cell_7_8 / 1000000;
% Cell_7_9=Cell_7_9 / 1000000;
% Cell_7_10=Cell_7_10 / 1000000;
% 
% % TT8
% a=size(CellNumbersTT8);
% 
% Cell_8_0=[];
% Cell_8_1=[];
% Cell_8_2=[];
% Cell_8_3=[];
% Cell_8_4=[];
% Cell_8_5=[];
% Cell_8_6=[];
% Cell_8_7=[];
% Cell_8_8=[];
% Cell_8_9=[];
% Cell_8_10=[];
% 
% for b=1:a(2)
%     if CellNumbersTT8(b)==0
%         Cell_8_0=[Cell_8_0 TimestampsTT8(b)];
%     else if CellNumbersTT8(b)==1
%             Cell_8_1=[Cell_8_1 TimestampsTT8(b)];
%         else if CellNumbersTT8(b)==2
%                 Cell_8_2=[Cell_8_2 TimestampsTT8(b)];
%             else if CellNumbersTT8(b)==3
%                     Cell_8_3=[Cell_8_3 TimestampsTT8(b)];
%                 else if CellNumbersTT8(b)==4
%                         Cell_8_4=[Cell_8_4 TimestampsTT8(b)];
%                     else if CellNumbersTT8(b)==5
%                             Cell_8_5=[Cell_8_5 TimestampsTT8(b)];
%                         else if CellNumbersTT8(b)==6
%                                 Cell_8_6=[Cell_8_6 TimestampsTT8(b)];
%                             else if CellNumbersTT8(b)==7
%                                     Cell_8_7=[Cell_8_7 TimestampsTT8(b)];
%                                 else if CellNumbersTT8(b)==8
%                                         Cell_8_8=[Cell_8_8 TimestampsTT8(b)];
%                                     else if CellNumbersTT8(b)==9
%                                             Cell_8_9=[Cell_8_9 TimestampsTT8(b)];
%                                         else if CellNumbersTT8(b)==10
%                                                 Cell_8_10=[Cell_8_10 TimestampsTT8(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_8_0=Cell_8_0 / 1000000;
% Cell_8_1=Cell_8_1 / 1000000;
% Cell_8_2=Cell_8_2 / 1000000;
% Cell_8_3=Cell_8_3 / 1000000;
% Cell_8_4=Cell_8_4 / 1000000;
% Cell_8_5=Cell_8_5 / 1000000;
% Cell_8_6=Cell_8_6 / 1000000;
% Cell_8_7=Cell_8_7 / 1000000;
% Cell_8_8=Cell_8_8 / 1000000;
% Cell_8_9=Cell_8_9 / 1000000;
% Cell_8_10=Cell_8_10 / 1000000;

%% second data set
%TT7
% a=size(CellNumbersTT52);
% 
% Cell_5_02=[];
% Cell_5_12=[];
% Cell_5_22=[];
% Cell_5_32=[];
% Cell_5_42=[];
% Cell_5_52=[];
% Cell_5_62=[];
% Cell_5_72=[];
% Cell_5_82=[];
% Cell_5_92=[];
% Cell_5_102=[];
% 
% for b=1:a(2)
%     if CellNumbersTT52(b)==0
%         Cell_5_02=[Cell_5_02 TimestampsTT52(b)];
%     else if CellNumbersTT52(b)==1
%             Cell_5_12=[Cell_5_12 TimestampsTT52(b)];
%         else if CellNumbersTT52(b)==2
%                 Cell_5_22=[Cell_5_22 TimestampsTT52(b)];
%             else if CellNumbersTT52(b)==3
%                     Cell_5_32=[Cell_5_32 TimestampsTT52(b)];
%                 else if CellNumbersTT52(b)==4
%                         Cell_5_42=[Cell_5_42 TimestampsTT52(b)];
%                     else if CellNumbersTT52(b)==5
%                             Cell_5_52=[Cell_5_52 TimestampsTT52(b)];
%                         else if CellNumbersTT52(b)==6
%                                 Cell_5_62=[Cell_5_62 TimestampsTT52(b)];
%                             else if CellNumbersTT52(b)==7
%                                     Cell_5_72=[Cell_5_72 TimestampsTT52(b)];
%                                 else if CellNumbersTT52(b)==8
%                                         Cell_5_82=[Cell_5_82 TimestampsTT52(b)];
%                                     else if CellNumbersTT52(b)==9
%                                             Cell_5_92=[Cell_5_92 TimestampsTT52(b)];
%                                         else if CellNumbersTT52(b)==10
%                                                 Cell_5_102=[Cell_5_102 TimestampsTT52(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_5_02=Cell_5_02 / 1000000;
% Cell_5_12=Cell_5_12 / 1000000;
% Cell_5_22=Cell_5_22 / 1000000;
% Cell_5_32=Cell_5_32 / 1000000;
% Cell_5_42=Cell_5_42 / 1000000;
% Cell_5_52=Cell_5_52 / 1000000;
% Cell_5_62=Cell_5_62 / 1000000;
% Cell_5_72=Cell_5_72 / 1000000;
% Cell_5_82=Cell_5_82 / 1000000;
% Cell_5_92=Cell_5_92 / 1000000;
% Cell_5_102=Cell_5_102 / 1000000;
% 
% %% third data set
% % TT2
% a=size(CellNumbersTT23);
% 
% Cell_2_03=[];
% Cell_2_13=[];
% Cell_2_23=[];
% Cell_2_33=[];
% Cell_2_43=[];
% Cell_2_53=[];
% Cell_2_63=[];
% Cell_2_73=[];
% Cell_2_83=[];
% Cell_2_93=[];
% Cell_2_103=[];
% 
% for b=1:a(2)
%     if CellNumbersTT23(b)==0
%         Cell_2_03=[Cell_2_03 TimestampsTT23(b)];
%     else if CellNumbersTT23(b)==1
%             Cell_2_13=[Cell_2_13 TimestampsTT23(b)];
%         else if CellNumbersTT23(b)==2
%                 Cell_2_23=[Cell_2_23 TimestampsTT23(b)];
%             else if CellNumbersTT23(b)==3
%                     Cell_2_33=[Cell_2_33 TimestampsTT23(b)];
%                 else if CellNumbersTT23(b)==4
%                         Cell_2_43=[Cell_2_43 TimestampsTT23(b)];
%                     else if CellNumbersTT23(b)==5
%                             Cell_2_53=[Cell_2_53 TimestampsTT23(b)];
%                         else if CellNumbersTT23(b)==6
%                                 Cell_2_63=[Cell_2_63 TimestampsTT23(b)];
%                             else if CellNumbersTT23(b)==7
%                                     Cell_2_73=[Cell_2_73 TimestampsTT23(b)];
%                                 else if CellNumbersTT23(b)==8
%                                         Cell_2_83=[Cell_2_83 TimestampsTT23(b)];
%                                     else if CellNumbersTT23(b)==9
%                                             Cell_2_93=[Cell_2_93 TimestampsTT23(b)];
%                                         else if CellNumbersTT23(b)==10
%                                                 Cell_2_103=[Cell_2_103 TimestampsTT23(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_2_03=Cell_2_03 / 1000000;
% Cell_2_13=Cell_2_13 / 1000000;
% Cell_2_23=Cell_2_23 / 1000000;
% Cell_2_33=Cell_2_33 / 1000000;
% Cell_2_43=Cell_2_43 / 1000000;
% Cell_2_53=Cell_2_53 / 1000000;
% Cell_2_63=Cell_2_63 / 1000000;
% Cell_2_73=Cell_2_73 / 1000000;
% Cell_2_83=Cell_2_83 / 1000000;
% Cell_2_93=Cell_2_93 / 1000000;
% Cell_2_103=Cell_2_103 / 1000000;
% 
% % TT3
% a=size(CellNumbersTT33);
% 
% Cell_3_03=[];
% Cell_3_13=[];
% Cell_3_23=[];
% Cell_3_33=[];
% Cell_3_43=[];
% Cell_3_53=[];
% Cell_3_63=[];
% Cell_3_73=[];
% Cell_3_83=[];
% Cell_3_93=[];
% Cell_3_103=[];
% 
% for b=1:a(2)
%     if CellNumbersTT33(b)==0
%         Cell_3_03=[Cell_3_03 TimestampsTT33(b)];
%     else if CellNumbersTT33(b)==1
%             Cell_3_13=[Cell_3_13 TimestampsTT33(b)];
%         else if CellNumbersTT33(b)==2
%                 Cell_3_23=[Cell_3_23 TimestampsTT33(b)];
%             else if CellNumbersTT33(b)==3
%                     Cell_3_33=[Cell_3_33 TimestampsTT33(b)];
%                 else if CellNumbersTT33(b)==4
%                         Cell_3_43=[Cell_3_43 TimestampsTT33(b)];
%                     else if CellNumbersTT33(b)==5
%                             Cell_3_53=[Cell_3_53 TimestampsTT33(b)];
%                         else if CellNumbersTT3(b)==6
%                                 Cell_3_63=[Cell_3_63 TimestampsTT33(b)];
%                             else if CellNumbersTT33(b)==7
%                                     Cell_3_73=[Cell_3_73 TimestampsTT33(b)];
%                                 else if CellNumbersTT33(b)==8
%                                         Cell_3_83=[Cell_3_83 TimestampsTT33(b)];
%                                     else if CellNumbersTT33(b)==9
%                                             Cell_3_93=[Cell_3_93 TimestampsTT33(b)];
%                                         else if CellNumbersTT33(b)==10
%                                                 Cell_3_103=[Cell_3_103 TimestampsTT33(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_3_03=Cell_3_03 / 1000000;
% Cell_3_13=Cell_3_13 / 1000000;
% Cell_3_23=Cell_3_23 / 1000000;
% Cell_3_33=Cell_3_33 / 1000000;
% Cell_3_43=Cell_3_43 / 1000000;
% Cell_3_53=Cell_3_53 / 1000000;
% Cell_3_63=Cell_3_63 / 1000000;
% Cell_3_73=Cell_3_73 / 1000000;
% Cell_3_83=Cell_3_83 / 1000000;
% Cell_3_93=Cell_3_93 / 1000000;
% Cell_3_103=Cell_3_103 / 1000000;
% 
% %%fourth data set
% % TT5
% a=size(CellNumbersTT54);
% 
% Cell_5_04=[];
% Cell_5_14=[];
% Cell_5_24=[];
% Cell_5_34=[];
% Cell_5_44=[];
% Cell_5_54=[];
% Cell_5_64=[];
% Cell_5_74=[];
% Cell_5_84=[];
% Cell_5_94=[];
% Cell_5_104=[];
% 
% for b=1:a(2)
%     if CellNumbersTT54(b)==0
%         Cell_5_04=[Cell_5_04 TimestampsTT54(b)];
%     else if CellNumbersTT54(b)==1
%             Cell_5_14=[Cell_5_14 TimestampsTT54(b)];
%         else if CellNumbersTT54(b)==2
%                 Cell_5_24=[Cell_5_24 TimestampsTT54(b)];
%             else if CellNumbersTT54(b)==3
%                     Cell_5_34=[Cell_5_34 TimestampsTT54(b)];
%                 else if CellNumbersTT54(b)==4
%                         Cell_5_44=[Cell_5_44 TimestampsTT54(b)];
%                     else if CellNumbersTT54(b)==5
%                             Cell_5_54=[Cell_5_54 TimestampsTT54(b)];
%                         else if CellNumbersTT54(b)==6
%                                 Cell_5_64=[Cell_5_64 TimestampsTT54(b)];
%                             else if CellNumbersTT54(b)==7
%                                     Cell_5_74=[Cell_5_74 TimestampsTT54(b)];
%                                 else if CellNumbersTT54(b)==8
%                                         Cell_5_84=[Cell_5_84 TimestampsTT54(b)];
%                                     else if CellNumbersTT54(b)==9
%                                             Cell_5_94=[Cell_5_94 TimestampsTT54(b)];
%                                         else if CellNumbersTT54(b)==10
%                                                 Cell_5_104=[Cell_5_104 TimestampsTT54(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_5_04=Cell_5_04 / 1000000;
% Cell_5_14=Cell_5_14 / 1000000;
% Cell_5_24=Cell_5_24 / 1000000;
% Cell_5_34=Cell_5_34 / 1000000;
% Cell_5_44=Cell_5_44 / 1000000;
% Cell_5_54=Cell_5_54 / 1000000;
% Cell_5_64=Cell_5_64 / 1000000;
% Cell_5_74=Cell_5_74 / 1000000;
% Cell_5_84=Cell_5_84 / 1000000;
% Cell_5_94=Cell_5_94 / 1000000;
% Cell_5_104=Cell_5_104 / 1000000;
% 
% %% fifth data set
% % TT8
% a=size(CellNumbersTT85);
% 
% Cell_8_05=[];
% Cell_8_15=[];
% Cell_8_25=[];
% Cell_8_35=[];
% Cell_8_45=[];
% Cell_8_55=[];
% Cell_8_65=[];
% Cell_8_75=[];
% Cell_8_85=[];
% Cell_8_95=[];
% Cell_8_105=[];
% 
% for b=1:a(2)
%     if CellNumbersTT85(b)==0
%         Cell_8_05=[Cell_8_05 TimestampsTT85(b)];
%     else if CellNumbersTT85(b)==1
%             Cell_8_15=[Cell_8_15 TimestampsTT85(b)];
%         else if CellNumbersTT85(b)==2
%                 Cell_8_25=[Cell_8_25 TimestampsTT85(b)];
%             else if CellNumbersTT85(b)==3
%                     Cell_8_35=[Cell_8_35 TimestampsTT85(b)];
%                 else if CellNumbersTT85(b)==4
%                         Cell_8_45=[Cell_8_45 TimestampsTT85(b)];
%                     else if CellNumbersTT85(b)==5
%                             Cell_8_55=[Cell_8_55 TimestampsTT85(b)];
%                         else if CellNumbersTT85(b)==6
%                                 Cell_8_65=[Cell_8_65 TimestampsTT85(b)];
%                             else if CellNumbersTT85(b)==7
%                                     Cell_8_75=[Cell_8_75 TimestampsTT85(b)];
%                                 else if CellNumbersTT85(b)==8
%                                         Cell_8_85=[Cell_8_85 TimestampsTT85(b)];
%                                     else if CellNumbersTT85(b)==9
%                                             Cell_8_95=[Cell_8_95 TimestampsTT85(b)];
%                                         else if CellNumbersTT85(b)==10
%                                                 Cell_8_105=[Cell_8_105 TimestampsTT85(b)];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% Cell_8_05=Cell_8_05 / 1000000;
% Cell_8_15=Cell_8_15 / 1000000;
% Cell_8_25=Cell_8_25 / 1000000;
% Cell_8_35=Cell_8_35 / 1000000;
% Cell_8_45=Cell_8_45 / 1000000;
% Cell_8_55=Cell_8_55 / 1000000;
% Cell_8_65=Cell_8_65 / 1000000;
% Cell_8_75=Cell_8_75 / 1000000;
% Cell_8_85=Cell_8_85 / 1000000;
% Cell_8_95=Cell_8_95 / 1000000;
% Cell_8_105=Cell_8_105 / 1000000;