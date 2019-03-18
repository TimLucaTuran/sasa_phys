function lastFileName = PARGRID_writeCurrentResults(storageData, storageDataNames, ...
            loopID, processID, indexField, currentIndex, DataStorageFileName, lastFileName, PARGRID_OriginalLoopVec)
%PARFOR_WRITECURRENTRESULTS Summary of this function goes here
%   Detailed explanation goes here

    if (iscell(DataStorageFileName))
        error('DataStorageFileName has to be a string when using in the <PARFOR_writeCurrentResults> function!');
    end

    storageDataNames = reshape(storageDataNames, [length(storageDataNames), 1]);
    
    % generate the current storage file for the current loop ...
    PARGRID_filename = PARGRID_createStorageFileName( DataStorageFileName, loopID, processID );

    if (~isempty(lastFileName))
        fdata = load('-mat', lastFileName);

        for d = 1 : length(storageData)
            dummy = fdata.(storageDataNames{d});
            fdata.(storageDataNames{d}) = cat(1, dummy, {storageData{d}});
        end
        
        fdata.PARGRID_Index = [fdata.PARGRID_Index; {currentIndex}];
    else
        fdata.PARGRID_OriginalLoopVec = PARGRID_OriginalLoopVec;
        fdata.PARGRID_OriginalIndexField = indexField;
        fdata.PARGRID_StorageDataNames = [storageDataNames; {'PARGRID_OriginalIndexField'}; {'PARGRID_Index'}];
        fdata.PARGRID_Index = {currentIndex};
        for d = 1 : length(storageData)
            fdata.(storageDataNames{d}) = {storageData{d}};
        end
    end
    
    tempFileChar = '#';
    if (isOctave)
        % ---------------------------------------------------------------------
        % On this point all variables which are needed in the following are
        % only: PARGRID_fData, PARGRID_fdataNames, PARGRID_filename,
        % PARGRID_counter, PARGRID_TempDummy. All others can be overwritten. Be
        % sure that no storageVariable-Name has one of the three reserved names!!!
        PARGRID_fData = fdata;
        clear fdata;

        % new version, because of Octave compatibility ...
        PARGRID_fdataNames = fieldnames(PARGRID_fData);
        for PARGRID_counter = 1 : length(PARGRID_fdataNames)
            PARGRID_TempDummy = PARGRID_fData.(PARGRID_fdataNames{PARGRID_counter});
            eval([PARGRID_fdataNames{PARGRID_counter}, '=PARGRID_TempDummy;']);
        end

        clear PARGRID_fData PARGRID_TempDummy PARGRID_counter;
        feval(@save, '-v7', [PARGRID_filename, tempFileChar], PARGRID_fdataNames{:});
        % ---------------------------------------------------------------------
    else    
        save('-mat', '-v7.3', [PARGRID_filename, tempFileChar], '-struct', 'fdata');
    end
    
    [status, message] = movefile([PARGRID_filename, tempFileChar], PARGRID_filename, 'f');
    
    if (~status)
        disp(['Problems to rename storage file with message: ', message]);
        lastFileName = [PARGRID_filename, tempFileChar];
    else
        lastFileName = PARGRID_filename;
    end
    
% endfunction