function localLoopVector = PARGRID_getCurrentLoopSection( loopVector, processNr, totalNumberOfProcesses )
%PARGRID_GETCURRENTLOOPSECTION Summary of this function goes here
%   Detailed explanation goes here

    LPieceSize = floor(numel(loopVector) / totalNumberOfProcesses);
    LRemainder = numel(loopVector) - LPieceSize*totalNumberOfProcesses;

    index = (processNr-1)*LPieceSize+1 : processNr*LPieceSize;
    if (processNr <= LRemainder)
        index = [index, processNr*LPieceSize + 1];
    end
    
    index = index + (min(processNr, LRemainder+1)-1);
    localLoopVector = loopVector(index);

end
