function out = easySRC(in, in_fs, out_fs, bw)
% This function uses a multirate converter object to convert arbitrary
% sample rates of time series data. Output rate tolerance is not exposed 
% but can be changed in line 22 below. Recommended to leave it at 0.0001.
%
% in        = Input signal, time series
% in_fs     = Sampling frequency of input signal
% out_fs    = Desired sampling frequency of output signal
% bw        = The two sided frequency range (DC centered) of the
%             information carrying portion of the output signal. Eg. for a
%             44.1kHz output signal, with information retention required 
%             up to 20kHz, enter 40 for "bw".
%
% Ben Jancovich, 2022
% University of New South Wales
% This work is licensed under the Creative Commons Attribution 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%
% Instantiate Sample Rate Converter (SRC) object
tol = 0.0001;
src = dsp.SampleRateConverter("InputSampleRate", in_fs,...
    "OutputSampleRate", out_fs, "Bandwidth", bw, "OutputRateTolerance",...
    tol);

% Compute decimation rate based on SRC settings
[~, M] = getRateChangeFactors(src);

% If length of signal "in" is NOT a multiple of decimation rate "M", 
% pad the end of signal "in" with zeros until it is.
in_len = length(in);
if mod(in_len, M) ~= 0
    remainder = mod(in_len, M);
    extrasamps = (M - remainder);
    in = [in; zeros(extrasamps, 1)];
    % Run SRC
    out = src(in);
    % Trim extra zeros
    trimsamps = ceil(extrasamps/M);
    out = out(1 : end-trimsamps);
else
    out = src(in);
end
% Kill SRC Object
release(src);
end