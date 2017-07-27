%  import EVs

function [Timestamps, EventIDs, TTLs, EventStrings] = importEV_v3

        FieldSelection(1) = 1;
        FieldSelection(2) = 1;
        FieldSelection(3) = 1;
        FieldSelection(4) = 0; %Extras
        FieldSelection(5) = 1;

        ExtractHeader = 0;
         
        ExtractMode = 1;
 
[Timestamps, EventIDs, TTLs, EventStrings] = Nlx2MatEV(eventfilename.name, FieldSelection, ExtractHeader,ExtractMode);
         
end
