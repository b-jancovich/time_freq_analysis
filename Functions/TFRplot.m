function TFRplot(t_vec, f_vec, TFR, freqlim, timelim, rotated)
% This function plots time-frequency representations as heat mapsusing the
% built-in MATLAB function "surf". Nothing interesting or novel, 
% just a shortcut for formatting nice plots.
%
% t_vec     = Time vector, must be same length as 'n' columns in TFR
% f_vec     = Frequency vector, must be same length as 'm' rows in TFR
% TFR       = Matrix of time frequency intensity [m x n]
% freqlim   = plot limits on frequency axis (same units as f_vec)
% timelim   = plot limits on time axis (same units as t_vec)
% rotated   = Optional. Some TFR algorithms return data with a different 
%               orientation. If data us upsidedown or otherwise not 
%               oriented as expected, try setting this to 1. Default is 0.
%
% Ben Jancovich, 2022
% University of New South Wales
% This work is licensed under the Creative Commons Attribution 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%
if nargin < 6
    rotated = 0;
end

switch rotated
    case 1
        rotate = [-90 -90];
    case 0
        rotate = [90 90];
end

surf(f_vec, t_vec, TFR', EdgeColor = 'none', FaceColor='texturemap');
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir='reverse', View=rotate, fontsize=12, FontName='Calibri');
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;