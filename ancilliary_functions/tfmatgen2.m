function [n_sweep_samps, n_silence_samps, matOUT] = tfmatgen2(fc1, fc2, f_vec, t_vec, ...
    duration_sweep, duration_silence)
% This function generates a time-frequency representaiton of a test signal 
% with known chracteristics, with size matching the sizes of arguments 
% t_vec and f_vec. No transforms or signal analyses are used. 
% 
% The test signal must be a linearly swept sinusoid with known start and 
% end frequencies, and duration. A period of silence with known duration 
% is expected before and after the sinusoid. The frequencies and durations 
% must be within the bounds of t_vec and f_vec.
%
% Inputs:
% fc1               = Start Frequency of carrier signal sweep (Hz)
% fc2               = End frequency of carrier signal sweep (Hz)
% f_vec             = Frequency vector corresponding to rows of matrix (Hz)
% t_vec             = Time vector corresponding to columns of matric (s)
% duration_sweep    = Duration of carrier and AM sweep signal (s)
% duration_silence  = Duration of silence before and after test signal (s)\
% freqscale         = Selects linear ('lin') or logarithmic ('log') 
%                       scaling for the frequency axis. Must match the
%                       scaling of f_vec. (string)
%
% Outputs:
% matOUT            = The time-frequency representation of the prescribed 
%                       signal. The units are magnitude squared (W).
% n_sweep_samps     = The number of samples (columns) of matOUT that
%                       contain the sweep
% n_silence_samps   = The number of samples (columns) of matOUT that
%                       contain the silence before the sweep
%
% Ben Jancovich, 2022
% Centre for Marine Science and Innovation
% School of Biological, Earth and Environmental Sciences
% University of New South Wales, Sydney, Australia
%
%% Prepare frequency and vector based on freqscale and freqdir
%
% Extend f_vec below fmin to handle frequencies below fmin. This prevents
% indexing errors when using negative frequencies.
% Anything generated in this range will be discarded and not included in
% analysis, so the <fmin part of this vector can be non-uniformly spaced 
% and differently scaled.
%
% Determine if f_vec ascends or descends:
if f_vec(1) < f_vec(end)
    freqdir = 'normal';
elseif f_vec(1) > f_vec(end)
    freqdir = 'reverse';
else
    error('ERROR - F_VEC DOES NOT ASCEND OR DESCEND');
end

% Make some extra rows to insert below fmin:
switch freqdir
    case 'normal'
        extrarows = linspace(-max(f_vec), min(f_vec)-1, length(f_vec))';
        f_vec = [extrarows; f_vec];
    case 'reverse'
        extrarows = linspace(min(f_vec)-1, -max(f_vec), length(f_vec))';
        f_vec = [f_vec; extrarows];
end

%% Construct TFR as table

% Determine size of final output matrix
rows_out = length(f_vec);
cols_out = length(t_vec);

% Convert freq and time vectors from numbers to strings for table indexing
f_labels = arrayfun(@(z) num2str(z, 15), f_vec, 'UniformOutput', 0);
t_labels = arrayfun(@(z) num2str(z, 15), t_vec, 'UniformOutput', 0);

% Create empty matrix with the desired dimensions
emptymat = zeros(rows_out, cols_out);

% Convert matrix to a table so we can index using non-integer freqs and times
table = array2table(emptymat, VariableNames=t_labels, RowNames=f_labels);

%% Generate sweep time index

% Discrete sampling means sweep does not start at precisely the perscribed time.
% Find time index closest to prescribed start time
[~,sweep_tidx_start] = min(abs(t_vec - duration_silence));

% Discrete sampling means sweep does not end at precisely the perscribed time.
% Find time index closest to prescribed end time
[~,sweep_tidx_end] = min(abs(t_vec - (duration_silence+duration_sweep)));

% Construct vector of sweep time indices
t_indices_sweep = t_vec(sweep_tidx_start:sweep_tidx_end);

%% Construct sweep frequency index

% Compressed frequency resolution means sweep may not start at precisely 
% the prescribed frequency. Find freq index closest to prescribed start
% frequency. 
[~,sweep_fidx_start] = min(abs(f_vec - fc1));
fc1 = f_vec(sweep_fidx_start);

% Compressed frequency resolution means sweep may not end at precisely 
% the perscribed frequency. Find freq index closest to prescribed end
% frequency. 
[~,sweep_fidx_end] = min(abs(f_vec - fc2));
fc2 = f_vec(sweep_fidx_end);

% Number of columns containing the sweep
n_cols_sweep = length(t_indices_sweep);

% Construct Vector of sweep frequency indices
f_indices_sweep_ideal = linspace(fc1, fc2, n_cols_sweep);

% Compressed frequency resolution means sweep f_indices may not align
% precisely with those of f_vec. Find closest freq indices in f_vec.
for i = 1:length(f_indices_sweep_ideal)
    [~,sweep_fidx_all(i)] = min(abs(f_vec - f_indices_sweep_ideal(i)));
end

% Extract closest matching indices from f_vec
f_indices_sweep = f_vec(sweep_fidx_all);

%% Write sweep data into table

% Convert sweep indices from numbers to strings for table indexing
f_sweep_labels = arrayfun(@(z) num2str(z, 15), f_indices_sweep, 'UniformOutput', 0);
t_sweep_labels = arrayfun(@(z) num2str(z, 15), t_indices_sweep, 'UniformOutput', 0);

% Write ones to the cells of the table at freq and time indices in
% f_indices_sweep and t_indices_sweep
for i=1:length(f_indices_sweep)
    f_index = f_sweep_labels(i);
    t_index = t_sweep_labels(i);
    table{f_index, t_index} = 1;
end

% Trim off the rows below fmin
extrarow_labels = arrayfun(@(z) num2str(z, 15), extrarows, 'UniformOutput', 0);
table(extrarow_labels,:) = [];

% Convert table back to a regular matrix
matOUT = table2array(table);

% Number of samples in sweep
n_sweep_samps = length(t_indices_sweep);

% Number of zeros to pad as silence
n_silence_samps = (cols_out - n_sweep_samps)/2;