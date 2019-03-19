function PARGRID_gatherData(pID_Min, pID_Max, loopID_Min, loopID_Max, filenameDataStorage, dataFieldsToExclude, indicesToExtractOnly)
%PARFOR_GATHERDATA Summary of this function goes here
%   Detailed explanation goes here

    if (nargin < 6)
        dataFieldsToExclude = {};
    end
    
    if (nargin < 7)
        indicesToExtractOnly = [];
    end
    
    if (~iscell(filenameDataStorage))
        filenameDataStorage = {filenameDataStorage};
    end
    
    for loopID = loopID_Min : loopID_Max
        for fn = 1 : length(filenameDataStorage)

            % generate the current storage file for the current loop ...
            filenameFinal = PARGRID_createStorageFileName( filenameDataStorage{fn}, loopID, 'FINAL');

            fDataFinal = struct();
            fDataFinalIsEmpty = true;

            for p = pID_Min : pID_Max

                disp(['GATHER: ', num2str(p), '/', num2str(pID_Max)]);
                filename = PARGRID_createStorageFileName( filenameDataStorage{fn}, loopID, p);

                if (exist(filename, 'file'))

                    % this delivers a cell array of the data vars stored in the
                    % mat-file ...
                    % varNamesInFile = who('-file', filename);
                    varNamesInFile = load('-mat', filename, 'PARGRID_StorageDataNames');
                    varNamesInFile = varNamesInFile.PARGRID_StorageDataNames;
                    PARGRID_OriginalLoopVec = load('-mat', filename, 'PARGRID_OriginalLoopVec');
                    PARGRID_OriginalLoopVec = PARGRID_OriginalLoopVec.PARGRID_OriginalLoopVec;
                    PARGRID_OriginalLoopVec = reshape(PARGRID_OriginalLoopVec, [length(PARGRID_OriginalLoopVec), 1]);

                    % remove the var-names that are included in the
                    % <dataFieldsToExclude> structure ...
                    for ex = 1 : length(dataFieldsToExclude)
                        cVar = 1;
                        while (cVar <= length(varNamesInFile))

                            % compare the strings ...
                            if (strcmp(varNamesInFile{cVar}, dataFieldsToExclude{ex}))
                                % delete element from list ...
                                varNamesInFile(cVar) = [];
                            else
                                cVar = cVar + 1;
                            end
                        end
                    end

                    %fdata = load('-mat', filename);
                    fdata = feval(@load, '-mat', filename, varNamesInFile{:});
                    fNames = fieldnames(fdata);

                    fDataFinalIsEmpty = false;

                    if (~isempty(indicesToExtractOnly))
                        indexToExtract = [];
                        for ex = 1 : length(indicesToExtractOnly)
                            indexToExtract = [indexToExtract, find(cell2mat(fdata.PARGRID_OriginalIndexField) == indicesToExtractOnly(ex))];
                        end
                    else
                        indexToExtract = 1:length(fdata.PARGRID_OriginalIndexField);
                    end

                    indexToSkip = 1:length(fdata.PARGRID_OriginalIndexField);
                    indexToSkip(indexToExtract) = [];

                    for f = 1 : length(fNames)

                        % if data entry is no cell array --> continue because
                        % there's nothing to merge ...
                        if (~iscell(fdata.(fNames{f})) && ~strcmp(fNames{f}, 'PARGRID_OriginalIndexField'))
                            continue;
                        end

                        fdata.(fNames{f})(indexToSkip) = [];

                        % get the current field of the fDataFinal structure ...
                        % if it is not exisiting then assume an empty array.
                        if (isfield(fDataFinal, fNames{f}))
                            cField = fDataFinal.(fNames{f});
                        else
                            if (iscell(fdata.(fNames{f})))
                                cField = {};
                            else
                                cField = [];
                            end
                        end

                        % merge data ...
                        if (~isempty(fdata.(fNames{f})))
                            cField = cat(1, cField, reshape(fdata.(fNames{f}), [length(fdata.(fNames{f})), 1]));
                        end                    

                        fDataFinal.(fNames{f}) = cField;
                    end

                    clear fdata;
                end
            end

            % only write the new file if the we have some data to write ...
            if (~fDataFinalIsEmpty)

                % now sort the fields inside PARFOR_Indexer ...
                [dummy, index] = sort(fDataFinal.PARGRID_OriginalIndexField);

                % resort all data in fDataFinal that it matches the new sorted
                % PARFOR_Indexer ...
                fNames = fieldnames(fDataFinal);
                for f = 1 : length(fNames)
                    cField = fDataFinal.(fNames{f});
                    fDataFinal.(fNames{f}) = cField(index);
                end

                % now check, if we have some multiple indices in PARFOR_Indexer.
                % This suggests that multiple processes worked on the same loop
                % index. But now we need only one result. So all multiple solutions
                % must be deleted ...
                i = 1;
                while (i < length(fDataFinal.PARGRID_OriginalIndexField))

                    % because PARFOR_Indexer is sorted it is sufficient to compare
                    % neighboured elements. If they are equal then we have a
                    % multiple data entry which must be deleted ...
                    if (fDataFinal.PARGRID_OriginalIndexField(i) == fDataFinal.PARGRID_OriginalIndexField(i+1))

                        % delete the entry ...
                        fNames = fieldnames(fDataFinal);
                        for f = 1 : length(fNames)
                            cField = fDataFinal.(fNames{f});
                            cField(i) = [];
                            fDataFinal.(fNames{f}) = cField;
                        end                                
                    else
                        i = i + 1;
                    end
                end

                % remove the field PARFOR_Indexer from the output data and save the
                % file ...
                fDataFinal = rmfield(fDataFinal, 'PARGRID_Index');

                % check if all "data points" were calculated, if some
                % process crashed, then it may by possible that not the
                % complete dataset is existing ...
                PARGRID_OriginalLoopVec = sort(PARGRID_OriginalLoopVec);
                allDataIsCalculated = isequal(fDataFinal.PARGRID_OriginalIndexField, PARGRID_OriginalLoopVec);
                if (~allDataIsCalculated)
                    disp(' ');
                    disp('WARNING: ');
                    disp('------------------------------------------------------------');
                    disp('The calculated data does not match the original data size!');
                    disp('Some data is missing! This might be caused by crashed');
                    disp('processes, which finished abnormally!');
                    disp(' ');
                end
                
                if (isOctave)
                    % ---------------------------------------------------------------------
                    % On this point all variables which are needed in the following are
                    % only: PARGRID_fData, PARGRID_fdataNames, PARGRID_filename,
                    % PARGRID_counter, PARGRID_TempDummy. All others can be overwritten. Be
                    % sure that no storageVariable-Name has one of the three reserved names!!!
                    PARGRID_fData = fDataFinal;
                    clear fDataFinal;

                    % new version, because of Octave compatibility ...
                    PARGRID_fdataNames = fieldnames(PARGRID_fData);
                    for PARGRID_counter = 1 : length(PARGRID_fdataNames)
                        PARGRID_TempDummy = PARGRID_fData.(PARGRID_fdataNames{PARGRID_counter});
                        eval([PARGRID_fdataNames{PARGRID_counter}, '=PARGRID_TempDummy;']);
                    end

                    clear PARGRID_fData PARGRID_TempDummy PARGRID_counter;
                    feval(@save, '-v7', filenameFinal, PARGRID_fdataNames{:});
                    % ---------------------------------------------------------------------
                else    
                    save('-mat', '-v7.3', filenameFinal, '-struct', 'fDataFinal');
                end
                
                                

                % delete all seperate processes storage files ...
                % but only if dataFieldsToExclude is empty -> otherwise we
                % would delete storaged data ...
                if (allDataIsCalculated && isempty(dataFieldsToExclude) && isempty(indicesToExtractOnly))
                    for p = pID_Min : pID_Max
                        filename = PARGRID_createStorageFileName( filenameDataStorage{fn}, loopID, p);

                        if (exist(filename, 'file'))
                            try
                                delete(filename);
                            catch
                                warning('Error occured while trying to delete partial storage files!');
                            end
                        end
                    end
                end
            end
        end
    end
    
% endfunction
