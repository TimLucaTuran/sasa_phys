function generateStartScriptsOCT(additionalArgs, codeToRunBefore)
% function generateStartScriptsOCT(additionalArgs, codeToRunBefore)
%
% THIS FUNCTION/SCRIPT WILL GENERATE ALL NECESSARY M-FILES AND
% LINUX-SCRIPTS WHICH ARE NEEDED TO RUN A MATLAB-SCRIPT IN THE
% <PARGRID>-ENVIRONMENT.
%
% The function for which the scripts will be generated is defined below
% (see <functionName>). This function must take at least three input
% arguments whereby the considered first three arguments are "reserved" for
% the action of the PARGRID environment. These first three parameters are:
%
%   PARGRID_processID ...
%       the number of the current process
%   PARGRID_totalNumberOfProcesses ...
%       total number of processes in which a special loop is divided in
%   PARGRID_storageID ...
%       the number of the storageFile which attributed to the current
%       process. Normally this number should be identical to
%       <PARGRID_processID>.
%
% If the function to be run needs additional arguments they can be appended
% after <PARGRID_storageID>.
%
% Usually a function/process is considered to have no addtional arguments,
% so that the <generateStartScriptsOCT>-function can be run without any
% argument. The values for the first three reserved parameters are
% automatically generated.
% If some process needs additonal parameters, then these parameters can be
% given in a cell array named <additionalArgs> (see argument list).
% For detailed info see below. If one needs to run any additonal code prior
% to the actual function/process call this code can be given in the
% parameter <codeToRunBefore>.
%
% additionalArgs ...
%   List of additonal parameters in form of a cell array. But the values
%   have to be given in form of a string so that they can be included into
%   the function call directly. What does it mean? Just to give an example:
%   If you want to give a number, let's say <7>, as the additonal argument
%   then <additionalArgs> is {'7'} and NOT {7}. If you have a more complex
%   argument which cannot be cast into a string that easy, then you have
%   the possibilty to calculate the value in code which is run before the
%   actual funcxtion call. See <codeToRunBefore> for detailed info.
%
% codeToRunBefore ...
%   List of additonal code lines in form of a cell array. The values
%   have to be given in form of a string so that they can be included into
%   the script directly. What does it mean? Just to give an example:
%   You want to give two random numbers as additonal arguments. Then you
%   can do something like: set <codeToRunBefore> to
%   {'arg1 = rand(100);', 'arg2 = rand(200)'} and set <additionalArgs> to:
%   {'arg1', 'arg2'}.


%%%%%%%%%%%%%%%%%%   TO BE SET BY USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
functionName = 'avg_wire_final';
storageName = functionName;
totalNumberOfProcesses = 128; 
diff_storageID_pID = 0;
        
% ---------------------------------------------------------------------
CALCULATIONPATH = '/lustre/vol/sperrhake/programming/calc/';
QUEUE_TO_SUBMIT_PARAS = '-d /lustre/vol/sperrhake/programming/calc/ -q fnadd5d -l ncpus=2,nodes=1';
PROCESSES_TO_RUN = 1:totalNumberOfProcesses;
UNZIP_CODEFILES = false;
REMOVE_STARTSCRIPTS = false;

% ---------------------------------------------------------------------
% define paths
lustre_source = [CALCULATIONPATH,'/FMM_source/'];
lustre_data = '/lustre/vol/sperrhake/programming/data/';
function_dest_dir = [lustre_data,functionName,'/'];

mkdir(lustre_data,functionName);
mkdir([lustre_data,functionName,'/FMM_source']);
mkdir([lustre_data,functionName,'/materials']);
mkdir([lustre_data,functionName,'/draw_files']);
% mkdir([lustre_data,functionName,'/build_files']);

copyfile([functionName,'.m'],[function_dest_dir,functionName,'.m']);
%copyfile('build_struct_from_particles.m',[function_dest_dir,'build_struct_from_particles.m'])
copyfile(lustre_source,[function_dest_dir,'FMM_source/']);
copyfile([CALCULATIONPATH,'materials/'],[function_dest_dir,'materials/']);
copyfile([CALCULATIONPATH,'draw_files/'],[function_dest_dir,'draw_files/']);
copyfile([mfilename,'.m'],function_dest_dir)

OCTAVE_PATH = '/cluster/bin/octave';
% OCTAVE_PATH = '/lustre/vol/kroll/install/bin/octave';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USE_NODEBIND = false;
NODEBIND_LIST = {'node13', 'node14', 'node15', 'node16'};
NODEBIND_MAXPROC = 8;

GATHER_LOOPMIN = 1;
GATHER_LOOPMAX = 1;
GATHER_FILE = [CALCULATIONPATH, storageName, '.mat'];
% GATHER_FILE = ['/home/sperrhake/programming/data', storageName, '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the default values for the input arguments if needed ...
if (nargin < 2)
    codeToRunBefore = {};
end
if (nargin < 1)
    additionalArgs = {};
end

% if <additionalArgs> or <codeToRunBefore> are strings then --> don't
% send an error but bracket them into a cell array.
if (ischar(codeToRunBefore))
    codeToRunBefore = {codeToRunBefore};
end
if (ischar(additionalArgs))
    additionalArgs = {additionalArgs};
end

% check if <additionalArgs> or <codeToRunBefore> are cell-arrays ...
if (~iscell(codeToRunBefore) || ~iscell(additionalArgs))
    error('<codeToRunBefore> and <additionalArgs> have to be cell arrays!');
end

% check if the contents of <additionalArgs> or <codeToRunBefore> are
% strings. If not --> send an error.
for aArgs = 1 : length(additionalArgs)
    if (~ischar(additionalArgs{aArgs}))
        error('All elements of <additionalArgs> have to be strings!');
    end
end
for cLines = 1 : length(codeToRunBefore)
    if (~ischar(codeToRunBefore{cLines}))
        error('All elements of <codeToRunBefore> have to be strings!');
    end
end

% ---------------------------------------------------------------------
% this will become the linux shell script which runs the qsubmission
octName = ['OCTSTART___', functionName, '.start'];

% include the unzip command ...
if (UNZIP_CODEFILES)
    linePATH = ['unzip -jo -d ', CALCULATIONPATH, ' GRID_SourceCode_DSRRNormal.zip'];
    dlmwrite(octName, linePATH, '');
else
    dlmwrite(octName, ' ', '');
end

% include the cd command ...
linePATH = ['cd ', CALCULATIONPATH];
dlmwrite(octName, linePATH, '-append', 'delimiter', '');
linePATH = ' ';
dlmwrite(octName, linePATH, '-append', 'delimiter', '');

% loop over all processes ...
NODEBIND_CURRENTNODE = 1;
cpu = 1;
for n = 1 : totalNumberOfProcesses
    
    % the GRIDSTART_N___ - script ...
    scriptName = ['GRIDSTART_', num2str(n) ,'___', functionName, '.m'];
    dlmwrite(scriptName, '', '');
    
    % include the <codeToRunBefore> ...
    for cLines = 1 : length(codeToRunBefore)
        dlmwrite(scriptName, codeToRunBefore{cLines}, '-append', 'delimiter', '');
    end
    
    % this is the actual function call to <functionName> ...
    lineText = [functionName, '(', num2str(n), ',', num2str(totalNumberOfProcesses), ...
        ',', num2str(n + diff_storageID_pID)];
    
    % include the <additionalArgs> ...
    for aArgs = 1 : length(additionalArgs)
        lineText = [lineText, ',', additionalArgs{aArgs}];
    end
    lineText = [lineText, ');'];
    dlmwrite(scriptName, lineText, '-append', 'delimiter', '');
    
    % now, check what processes have to be run ...
    if (~isempty(find(PROCESSES_TO_RUN == n, 1)))
        
        % this will become the startscript which is send to qsub ...
        tempOCTScript = ['OS___', functionName, '_#', num2str(n), '___', num2str(datenum(clock)), '.start'];
        startLine = ['echo "cd ', CALCULATIONPATH, '" > ', tempOCTScript];
        dlmwrite(octName, startLine, '-append', 'delimiter', '');
        startLine = ['echo "',OCTAVE_PATH,' < ', scriptName, '" >> ', tempOCTScript];
        dlmwrite(octName, startLine, '-append', 'delimiter', '');
        
        % append the qsub call to the start-script ...
        if (USE_NODEBIND)
            startLine = ['qsub ', appendNodeBind(QUEUE_TO_SUBMIT_PARAS, ['nodes=', NODEBIND_LIST{NODEBIND_CURRENTNODE} ]), ' ', tempOCTScript];
        else
            startLine = ['qsub ', QUEUE_TO_SUBMIT_PARAS, ' ', tempOCTScript];
        end
        dlmwrite(octName, startLine, '-append', 'delimiter', '');
        
        if (REMOVE_STARTSCRIPTS)
            % remove the temporary start script files ...
            startLine = ['rm ', tempOCTScript];
            dlmwrite(octName, startLine, '-append', 'delimiter', '');
        end
        
        startLine = ' ';
        dlmwrite(octName, startLine, '-append', 'delimiter', '');
        
        if (cpu >= NODEBIND_MAXPROC)
            cpu = 1;
            NODEBIND_CURRENTNODE = NODEBIND_CURRENTNODE + 1;
            if (NODEBIND_CURRENTNODE > length(NODEBIND_LIST))
                warning('<NODEBIND_CURRENTNODE> exceeds limits of <NODEBIND_LIST>!');
            end
            NODEBIND_CURRENTNODE = min(NODEBIND_CURRENTNODE, length(NODEBIND_LIST));
        else
            cpu = cpu+1;
        end
    end
end

% write also the GRIDGATHER script ...
scriptName = ['GRIDGATHER___', functionName, '.m'];
startLine = ['PARGRID_gatherData(', num2str(1), ',' num2str(totalNumberOfProcesses), ',', ...
    num2str(GATHER_LOOPMIN), ',', num2str(GATHER_LOOPMAX), ',''', GATHER_FILE ''')'];
dlmwrite(scriptName, startLine, '');

% write the linux shell script to run the GATHER process ...
octName = ['OCTSTART___', functionName, '.gather'];
startLine = ['cd ', CALCULATIONPATH];
dlmwrite(octName, startLine, '');
startLine = [OCTAVE_PATH,' < ', scriptName];
dlmwrite(octName, startLine, '-append', 'delimiter', '');

% endfunction


function QUEUE_TO_SUBMIT_PARAS = appendNodeBind(QUEUE_TO_SUBMIT_PARAS, strToAppend)

loc = strfind(QUEUE_TO_SUBMIT_PARAS, '-l');
if (~isempty(loc))
    loc = loc(1);
    token = QUEUE_TO_SUBMIT_PARAS(1 : loc+1);
    remain = QUEUE_TO_SUBMIT_PARAS(loc+2 : end);
    
    [token2, remain2] = strtok(remain, '-');
    token2 = [' ', strtrim(token2), ',', strtrim(strToAppend), ' '];
    
    QUEUE_TO_SUBMIT_PARAS = [token, token2, remain2];
else
    QUEUE_TO_SUBMIT_PARAS = [strtrim(QUEUE_TO_SUBMIT_PARAS), ' -l ', strToAppend];
end

% endfunction
