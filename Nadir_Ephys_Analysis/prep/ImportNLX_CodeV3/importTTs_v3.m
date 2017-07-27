%  import TTs

function [Timestamp, CellNumbers, DataPoints] = importTTs_v3 (FileName)

        FieldSelection(1) = 1;
        FieldSelection(2) = 0; %ScNumbers
        FieldSelection(3) = 1;
        FieldSelection(4) = 0; %Params 
        FieldSelection(5) = 1;

        ExtractHeader = 0;
         
        ExtractMode = 1;
 
[Timestamp, CellNumbers, DataPoints] = Nlx2MatSpike_v3(FileName, FieldSelection, ExtractHeader,ExtractMode);
         
%end
