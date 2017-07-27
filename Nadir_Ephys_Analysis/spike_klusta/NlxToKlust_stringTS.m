% # modification of NlxToKlust that uses strings at the end of the file names instead of
% # timestamps for reading in appended files from Cheetah 5.7.2 (i.e. CSC1_0002)
% # example call: KlustaData = NlxToKlust_stringTS ([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16], '0001', 'kwikdata1to16.dat', 0)

function done = NlxToKlust(desiredCSCs, endString, filename, platform);

    myid = fopen(filename,'w');
    chunkSize = 1024;       % 1 megabyte divided by 512, in bits, divided by 16
    
    if isempty(endString)
        tempfilename=dir(['CSC' num2str(desiredCSCs(1)) '.ncs']);
    else
        tempfilename=dir(['CSC' num2str(desiredCSCs(1)) '_' endString '*.ncs']);
    end
    
    if platform==1
        [tempTimestamp,tempSampleFrequency,tempSamples] = importCSC_v3_PC(tempfilename.name);
    else
        [tempTimestamp,tempSampleFrequency,tempSamples] = importCSC_v3(tempfilename.name);
    end
    
    numSamples = length(tempSamples);
    
    numChunks = floor(numSamples / chunkSize);
    lastChunkStart = (numChunks * chunkSize) + 1;
    disp('Records:');
    disp(numSamples);
    disp('Number of chunks:');
    disp(numChunks);
    disp('Index of start of last chunk:');
    disp(lastChunkStart);
    
    for i = 1:chunkSize:(numChunks*chunkSize)      
        clear tempBlock;
        tempBlock = zeros(length(desiredCSCs),512,chunkSize);
        disp('Percent complete:');
        disp((i/numSamples)*100);
        for j=1:length(desiredCSCs)
            if isempty(endString)
                tempfilename=dir(['CSC' num2str(desiredCSCs(j)) '.ncs']);
            else
                tempfilename=dir(['CSC' num2str(desiredCSCs(j)) '_' endString '*.ncs']);
            end
            if platform==1
                [tempBlock(j,:,:)] = Nlx2MatCSC(tempfilename.name, [0 0 0 0 1], 0, 2, [i (i+chunkSize-1)]);
            else
                [tempBlock(j,:,:)] = Nlx2MatCSC_v3(tempfilename.name, [0 0 0 0 1], 0, 2, [i (i+chunkSize-1)]);
            end
        end
        
        tempSamples = reshape(tempBlock,length(desiredCSCs),512*chunkSize);
        fwrite(myid,tempSamples,'int16');
    end
    
    clear tempBlock;
    if platform==1
        tempBlock = zeros(length(desiredCSCs),512,(numSamples - lastChunkStart + 1));
    else
        tempBlock = zeros(length(desiredCSCs),512,(numSamples - lastChunkStart));
    end
    
    %disp(size(tempBlock));
	for j=1:length(desiredCSCs)
        if isempty(endString)
            tempfilename=dir(['CSC' num2str(desiredCSCs(j)) '.ncs']);
        else
            tempfilename=dir(['CSC' num2str(desiredCSCs(j)) '_' endString '*.ncs']);
        end
        if platform==1
            [tempBlock(j,:,:)] = Nlx2MatCSC(tempfilename.name, [0 0 0 0 1], 0, 2, [lastChunkStart numSamples]);
            % disp(size(tempBlock(j,:,:)));
        else
            [tempBlock(j,:,:)] = Nlx2MatCSC_v3(tempfilename.name, [0 0 0 0 1], 0, 2, [lastChunkStart numSamples]);
           % [temporary] = Nlx2MatCSC_v3(tempfilename.name, [0 0 0 0 1], 0, 2, [lastChunkStart numSamples]);
           % disp(size(temporary));
           % tempBlock(j,:,:) = temporary;
        end
    end
    
    clear tempSamples;
    if platform==1
        tempSamples = reshape(tempBlock,length(desiredCSCs),512*(numSamples - lastChunkStart + 1));
    else
        tempSamples = reshape(tempBlock,length(desiredCSCs),512*(numSamples - lastChunkStart));
    end
    
	fwrite(myid,tempSamples,'int16');

    
    fclose(myid);
    
    done = 1;
end