function OUTMAT = sweep_time_freq_mat(f1, f2, fres, nsamps, amp, rowsOUT)

% This function is for constructing time-frequency representations of
% linearly swept sinusoids. It arranges 'amp' onto a matrix 'OUTMAT' 
% according to the frequency coordinates defined by f1, f2, fres, and 
% across time duration set by nsamps. All other values in 'OUTMAT' are 
% zero. f is the vector of frequency sampling points along the y axis 
% of outTMAT. szOUT is the size of the final matrix outMAT. If this is 
% larger than MAT, zeros will be padded to achieve szOUT.

% Inputs:
% f1 =      Sweep starting frequency (Hz, >0)
% f2 =      Sweep ending frequency (Hz, >0)
% fres =    Frequency resolution (Hz)
% nsamps =  Duration for sweep to traverse from f1 to f2 (# of samples).
% amp =     Vector of amplitude values to enter to outMAT
% rowsOUT = N of rows of the final output array. MAT will be padded with 
%           zeros to meet this size. int. > nf

% Outputs:
% OUTMAT =  Matrix of time-frequency-intensity data.
%           Dimensions are [nf x length(amp) x (0:1)]

assert(f1>0 || f2>0, 'ERROR: VARIABLES "f1" & "f2" MUST BE POSITIVE')

% Create Vector of sample indices
t = 1:1:nsamps;         

% Create Vector of frequency indices
if f1 == f2 % if f1 is equal to f2, there is only one freq index.
    nf = 1;                             
    f_indices = ones(1, nf) .* f1;
            fmin = f1;
elseif f1 > f2 % case where the sweep decends in freq.
            nf = (f1 - f2) / fres;     
            fmin = f2;
            f_indices = linspace(1, nf, nf);
elseif f1 < f2 % case where the sweep ascends in freq.
            nf = (f2 - f1) / fres;             
            fmin = f1;
            f_indices = linspace(nf, 1, nf);    
end

assert(rowsOUT >= nf, 'ERROR: VARIABLES "rowsOUT" MUST BE > (range([f1, f2]) / fres)')

% ratio of 't' to 'f'
r = floor(nsamps / nf);

% Hold each frequency coordinate for 'r' time samples
f_reshaped = repelem(f_indices, r);

if length(f_reshaped) < nsamps
    diff = nsamps-length(f_reshaped);
    xtra = linspace(f_reshaped(end)+1, f_reshaped(end)+floor((diff/r)), floor(diff/r));
    xtra = repelem(xtra, r);
    f_reshaped = [f_reshaped, xtra];
    if length(f_reshaped) < nsamps
        diff = nsamps-length(f_reshaped);
        f_reshaped = [f_reshaped, f_reshaped(end) * ones(1, diff)];
    end
elseif length(f_reshaped) > nsamps
    killsamps = length(f_reshaped) - nsamps;
    rkill = ceil(killsamps / r);
    f_reshaped(rkill:rkill:end,:) = [];
    f_reshaped = f_reshaped(1: end-killsamps);
else
end

% Build Matrix
if f1 ~= f2
    % Accumulate 'amp' to matrix according to coordinates in t and f
    MAT = accumarray([f_reshaped', t'], amp, [max(f_reshaped), nsamps]); % the max(f_reshaped) was 'nf' before but changed to handle the new "if f_reshaped is short" stuff
elseif f1 == f2
    MAT = amp .* (ones(1, nsamps));
end

% Add empty rows to complete frequency range as per rowsOUT
n_bottom_rows = fmin / fres;
bottom_rows = zeros(n_bottom_rows, nsamps);
n_top_rows = rowsOUT - (size(MAT, 1) + n_bottom_rows);
top_rows = zeros(n_top_rows, nsamps);
OUTMAT = [top_rows; MAT; bottom_rows];

end