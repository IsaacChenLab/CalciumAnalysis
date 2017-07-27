
%  extract CSC samples from CSC data
 
function [Timestamp, SampleFrequency, Samples] = importCSC(FileName)

        
        FieldSelection = [1,0,1,1,1];

        ExtractHeader = 0;
         
        ExtractMode = 1;
        
        test = FileName;

 
[Timestamp, SampleFrequency, ~, Samples] = Nlx2MatCSC(FileName, FieldSelection, ExtractHeader,ExtractMode,[]);
         
end
