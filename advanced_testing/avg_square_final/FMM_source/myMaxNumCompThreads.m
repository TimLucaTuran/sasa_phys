function N = myMaxNumCompThreads( varargin )
%MYMAXNUMCOMPTHREADS is identical in its behaviour like the MATLAB internal
%function <maxNumCompThreads>. The current version only provides the
%extension to Octave. Difference: If set to 'automatic', then in Octave
%the number of threads is always set to 1, whereby in MATLAB the number is
%chosen to be equal to the core number of the machine the process is
%running on. 

    if (nargin > 1)
        error('myMaxNumCompThreads:nargin', 'The number of input arguments is limited to 1!');
    end

    if (isOctave)
        N = str2num(getenv('OMP_NUM_THREADS'));
        if (nargin == 1)
            if (ischar(varargin{1}))
                if (strcmpi(varargin{1}, 'automatic'))
                    setenv('OMP_NUM_THREADS', '1');
                else
                    error('myMaxNumCompThreads:automatic', 'The type of input argument should be integer or ''automatic''!');
                end
            elseif (isreal(varargin{1}) && varargin{1} >= 1)
                num = round(varargin{1});
                setenv('OMP_NUM_THREADS', num2str(num));
            else
                error('myMaxNumCompThreads:argtype', 'The type of input argument should be postive integer or ''automatic''!');
            end
        end
        
    else
        if (nargin == 0)
            N = maxNumCompThreads;
        else
            N = maxNumCompThreads(varargin{1});
        end
    end
end
