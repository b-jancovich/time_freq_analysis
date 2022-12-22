function matOUT = TFRunitconvertNorm(matIN, units)
% This function ingests a matrix time-frequency representation and converts
% the intensity units according to the argument "units". The result is then
% normalized to max(matIN) = 1.

% Input:
% matIN = [NxM] matrix of time-frequency data (expected format is real magnitude)
% units = string, sets the units of the output, matOUT.
%           Options: 
%               'magnitude'   - Units remain unchanged, function only applies normalization.
%               'magdb' - Magnitude, expressed on a log scale, ie Decibels.
%               'power' - Power expressed on a linear scale, ie watts.
%               'powdb' - Power expressed on a logarithmic scale, ie dBW.

% Check input is real.
assert(isreal(matIN), 'Error, input matrix must be real magnitude, not complex values.')

% Normalize to max = 1
mat = matIN ./ max((matIN), [], 'all');

% Convert units
switch units
    case 'magnitude'
        matOUT = mat;
    case 'magdb'
        matOUT = 20*log10(mat); 
    case 'power'
        matOUT = mat .^2;  
    case 'powdb'
        matOUT = 10*log10(mat .^2);   
end
end