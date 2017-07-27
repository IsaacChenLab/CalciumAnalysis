%  import EVs

function [Timestamps, EventIDs, TTLs, EventStrings] = importEV;

        FieldSelection(1) = 1;
        FieldSelection(2) = 1;
        FieldSelection(3) = 1;
        FieldSelection(4) = 0; %Extras
        FieldSelection(5) = 1;

        ExtractHeader = 0;
         
        ExtractMode = 1;
 
[Timestamps, EventIDs, TTLs, EventStrings] = Nlx2MatEV('Events.nev', FieldSelection, ExtractHeader,ExtractMode);
         
end
