function sFilename = PARGRID_createStorageFileName( baseName, loopID, pID )
%PARFOR_CREATECONTROLFILENAME Summary of this function goes here
%   Detailed explanation goes here

    if (ischar(pID) && strcmp(upper(pID), 'FINAL'))
        [pathstr, name, ext] = fileparts(baseName);
        if (isempty(pathstr)); pathstr = '.'; end;
        sFilename = [pathstr, filesep, name, '___LoopID#', num2str(loopID), '___FINAL', ext];
    else
        [pathstr, name, ext] = fileparts(baseName);
        if (isempty(pathstr)); pathstr = '.'; end;
        sFilename = [pathstr, filesep, name, '___LoopID#', num2str(loopID), '___pID#', num2str(pID), ext];
    end
    
% endfunction
    