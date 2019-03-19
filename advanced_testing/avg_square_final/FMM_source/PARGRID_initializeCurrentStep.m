function [PARGRID_loopID, PARGRID_DataStorageFileName, PARGRID_lastFileName, PARGRID_CurrentLoopVec, PARGRID_OriginalLoopVec] = ...
    PARGRID_initializeCurrentStep( PARGRID_loopID, PARGRID_processID, ...
    PARGRID_totalNumberOfProcesses, loopVector, storageFileName)
%PARGRID_INITIALIZECURRENTSTEP Summary of this function goes here
%   Detailed explanation goes here

    if (~iscell(loopVector))
        PARGRID_DataStorageFileName = storageFileName;
        PARGRID_lastFileName = '';
        PARGRID_loopID = PARGRID_loopID + 1;    
        PARGRID_CurrentLoopVec = PARGRID_getCurrentLoopSection( loopVector, PARGRID_processID, PARGRID_totalNumberOfProcesses );
        PARGRID_OriginalLoopVec = loopVector;
    else
        PARGRID_DataStorageFileName = storageFileName;
        PARGRID_lastFileName = '';
        PARGRID_loopID = PARGRID_loopID + 1;  
        PARGRID_CurrentLoopVec = cell(length(loopVector));
        for n = 1 : length(loopVector)
            PARGRID_CurrentLoopVec{n} = PARGRID_getCurrentLoopSection( loopVector{n}, PARGRID_processID, PARGRID_totalNumberOfProcesses );
        end
        PARGRID_OriginalLoopVec = loopVector{1};        
    end
end
