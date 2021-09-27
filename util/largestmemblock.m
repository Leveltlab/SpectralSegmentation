function largestmemblock(varargin)
% LARGESETMEMBLOCK(VERBOSITY,PRECISION) 
% Finds (roughly) the largest memory block available by
% trying to allocate successively larger contiguous blocks of memory.
% Precise roughly to about 90% of the largest memory block available.
% 
% Arguments:
%    VERBOSITY: '-v' flag. Specify all iterations tried
%    PRECISION: Fraction between [0.7, 0.9999]. 
%      Value closer to 0.7 return a quicker result, but up to 70% precision
%      Value closer to 1 returns a more precise result, but takes much
%      longer to converge
%      Default is 0.95
% 
%  Not as accurate as feature('memstats'), but useful as an approximation
%  on 32-bit and 64-bit Mac and Linux distributions
% 
%  % No verbosity
%  >> largestmemblock   
% 
%  % Verbose output
%  >> largestmemblock('-v')
% 
%  % More accurate estimate of largest memory block available, no verbosity
%  >> largestmemblock([],0.999)
% 
%  % More accurate value, longer execution time
% 
%  v1.0 - November 2009, Manu Raghavan

if nargin>=1 && ~isempty(varargin{1})
    verbosity = true;
else
    verbosity = false;  
end

precision = 0.95;
if nargin==2
    precision = varargin{2};
    if(precision>0.9999 || precision<0.7)
        error('Optional second argument PRECISION must range between 0.7 and 0.9999');
    end
end

mexExtension = mexext;

if (strcmp(mexExtension,'mexglx') || ... 
    strcmp(mexExtension,'mexmaci') || ...
    strcmp(mexExtension,'mexw32'))    
    % Initial size
    B = intmax('int32'); %2 Gigabytes    
end

if (strcmp(mexExtension,'mexa64') || ... 
    strcmp(mexExtension,'mexmaci64') || ...
    strcmp(mexExtension,'mexw64'))    
    % Initial size
    %B = 8*2^40; %8 Terabytes
    B = 60*2^30; %8 Gigabytes
end

% Iterative memory allocation attempts
while true
    try   
        if(verbosity)
            fprintf(1,'Trying to allocate %d KB.', B/1024);
        end
        A = zeros(round(B), 1, 'uint8');        
        break;
    catch ME
        if(verbosity)
            fprintf(2,'\n   Unsuccessful. Reducing block size and trying again...\n');
        end
        % Rethrow any non-out-of-memory errors
        if ~strcmp(ME.identifier, 'MATLAB:nomem');
            rethrow(ME);
        end
    end
    % Reduce by a factor of 10% and try again
    B = B*precision;
end

% Report the results
fprintf(...
    '\nMaximum contiguous block is approximately:\n (Accurate to %2.2f %% of real memory value)\n    %15.0f KB\n    %15.0f MB\n', ...
    precision*100, B/1024, B/(1024^2));

% Clear the big variable
clear A
  