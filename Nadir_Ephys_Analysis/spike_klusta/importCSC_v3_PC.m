%  extract CSC samples from CSC data
 
function [Timestamp, SampleFrequency, Samples] = importCSC_v3 (FileName)

        FieldSelection(1) = 1;
        FieldSelection(2) = 0;
        FieldSelection(3) = 1;
        FieldSelection(4) = 1;
        FieldSelection(5) = 1;

        ExtractHeader = 0;
         
        ExtractMode = 1;
 
[Timestamp, SampleFrequency, NumValSamples, Samples] = Nlx2MatCSC(FileName, FieldSelection, ExtractHeader,ExtractMode);

         
end
