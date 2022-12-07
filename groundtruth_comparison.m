%% Comparison of Time-Frequency Analysis Algorithms with Ground Truth.

% This script generates a carrier signal (linearly frequency swept sine)
% and an amplitude modulation signal (linearly frequency swept square)
% and multiplies them together, creating an amplitude modulated test signal.
% It then constructs a matrix that can be thought of as a synthetic
% spectrogram representing the "Ground Truth" of the test signal's
% power as a function of both time & frequency. This spectrogram is
% constructed analytically using knowledge of the test signal's TF
% charactertistics, rather than by transforming or otherwise analysing
% the signal. This ground truth matrix serves as a theoretically perfect
% spectrogram of the test signal. Several transform based TF analyses are
% then performed and compared against the ground truth.

% Note: if "carrier_freq1" x "mod_freq1", "carrier_freq2" x "mod_freq2",
% or any of the products of their first 4 subharmonics reach a negative
% value, this script will break.

% Methods compared:
% Short Time Fourier Transform (short window)
% Short Time Fourier Transform (long window)
% Continuous Wavelet Transform
% Fractional Adaptive Superresolution Wavelet Transform

% Ben Jancovich, 2022
% Centre for Marine Science and Innovation
% School of Biological, Earth and Environmental Sciences
% University of New South Wales, Sydney, Australia

clear
close
clc

%% User Variables

% STFT FFT size is locked at n_freqs_total*2 = 1250. This results in equal
% number of frequency points between fmin and fmax for all methods. This is set
% using fmin, fmax and fres for other methods.

% Test Signal Parameters:
carrier_freq1 = 50;         % Sine Sweep start frequency (Hz)
carrier_freq2 = 30;         % Sine Sweep end frequency (Hz)
mod_freq1 = 2;              % Amplitude Modulation sweep start frequency (Hz)
mod_freq2 = 6;              % Amplitude Modulation sweep end frequency (Hz)
duration_sweep = 5;         % Sweep Duration
silence = 1;                % Silence to pad start and end (seconds)
fs = 250;                   % Sampling Frequency (Hz)

% Ground Truth Parameters:
sigma = 1;                  % Standard deviation of gaussian filter

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest
fmax = fs/2;                % Highest frequency of interest
f_res = 0.2;                % Frequency resolution (Hz)
powerscaling = 'lin';       % Plot power as 'lin' (W) or 'log' (dBW)  

% STFT Window Sizes
n_fft = 2*((fs/2)/f_res);   % spectrogram FFT length (samples)
win1 = n_fft/5;             % Window length for long STFT (samples)
win2 = 50;                  % Window length for short STFT (samples)
overlap1 = 75;              % Window overlap % for long STFT
overlap2 = 75;              % Window Overlap % for short STFT

% CWT Parameters
tbp = 120;      % Time bandwidth product of Morse wavelet.
gamma = 3;      % Symmetry of the wavelet. 3 = symmetric.
vpo = 48;       % Voices per Octave. Must be in range 10 : 48

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 50];    % Interval of superresolution orders
mult = 1;           % Multiplicative ('1') or additive ('0') superresolution
dcfilt = 10;        % Cutoff of DC filter (Hz)

%% Generate time & Frequency Vectors

% Frequency & samp counts
n_samps_sweep = duration_sweep * fs;
n_samps_total = n_samps_sweep + (2*(silence*fs));
n_freqs_sweep = (carrier_freq1-carrier_freq2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vectors
t_vec_sweep = linspace(0, n_samps_sweep/fs, n_samps_sweep);
t_vec_total = linspace(0, n_samps_total/fs, n_samps_total);

% Frequency vectors
f_vec_sweep = linspace(carrier_freq1, carrier_freq2, n_freqs_sweep);
f_vec_total = linspace(fmin, fmax, n_freqs_total);

%% Generate signals

% Generate Sine Sweep (Carrier Signal)
signal = chirp(t_vec_sweep, carrier_freq1, t_vec_sweep(end), carrier_freq2, 'linear');

% Generate Square Sweep (Amplitude Modulator)
mod = rescale(sign(chirp(t_vec_sweep, mod_freq1, t_vec_sweep(end), mod_freq2, 'linear')));

% Pad signal and mod with zeros to insert silence at start & end.
signal_sil = [zeros(1, silence * fs), signal, zeros(1, silence * fs)];
mod_sil = [zeros(1, silence * fs), mod, zeros(1, silence * fs)];

% Amplitude Modulation
signal = signal_sil .* mod_sil;

%% Ground Truth - Time-Freq Matrix

% Calculate Start and End Frequencies for Upper Sideband Components
usb1_f1 = carrier_freq1 + mod_freq1;
usb1_f2 = carrier_freq2 + mod_freq2;
usb2_f1 = carrier_freq1 + (mod_freq1 * 3);
usb2_f2 = carrier_freq2 + (mod_freq2 * 3);
usb3_f1 = carrier_freq1 + (mod_freq1 * 5);
usb3_f2 = carrier_freq2 + (mod_freq2 * 5);

% Calculate Start and End Frequencies for Lower Sideband Components
lsb1_f1 = carrier_freq1 - mod_freq1;
lsb1_f2 = carrier_freq2 - mod_freq2;
lsb2_f1 = carrier_freq1 - (mod_freq1 * 3);
lsb2_f2 = carrier_freq2 - (mod_freq2 * 3);
lsb3_f1 = carrier_freq1 - (mod_freq1 * 5);
lsb3_f2 = carrier_freq2 - (mod_freq2 * 5);

% Error handling
sideband_freqs = [usb1_f1, usb1_f2, usb2_f1, usb2_f2, usb3_f1, usb3_f2,...
    lsb1_f1, lsb1_f2, lsb2_f1, lsb2_f2, lsb3_f1, lsb3_f2];
assert(all(sideband_freqs >= 0), ['ERROR: THE PRODUCT OF CARRIER AND ' ...
    'MODULATOR FREQUENCIES YOU HAVE ENTERED HAS RESULTED IN A FREQUENCY ' ...
    'COMPONENT WITH NEGATIVE FREQUENCY. THIS IS NOT SUPPORTED. PLEASE ' ...
    'INCREASE LOWEST CARRIER OR MODULATOR FREQUENCY'])

% Build matrices representing carrier and sidebands
% (f1, f2, f_res, nsamps, amp, rowsOUT)
amp = ones(1, n_samps_sweep);

carrier_MAT = tfmatgen(carrier_freq1, carrier_freq2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb1_MAT = tfmatgen(usb1_f1, usb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb2_MAT = tfmatgen(usb2_f1, usb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb3_MAT = tfmatgen(usb3_f1, usb3_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb1_MAT = tfmatgen(lsb1_f1, lsb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb2_MAT = tfmatgen(lsb2_f1, lsb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb3_MAT = tfmatgen(lsb3_f1, lsb3_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);

% Combine matrices for all components
all_MAT = carrier_MAT + (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + (usb3_MAT .* 0.125) +...
    (lsb1_MAT .* 0.5) + (lsb2_MAT .* 0.25) +...
    (lsb3_MAT .* 0.125);

% Due to fmin being constant=0 in tfmatgen.m, and a variable in this
% script, the ground truth matrix is shifted up in frequency by the value
% of fmin/f_res+1. Use circshift to shift it back down.
all_MAT = circshift(flip(all_MAT, 1), -(fmin/f_res)-1, 1);

% That circshift moves some data that is below fmin to the high frequencies.
% It shouldn't be there. Zero it out.
all_MAT(end-(fmin/f_res):end, :) = 0;

% Apply amplitude modulation
modMAT = repmat(mod,[n_freqs_total, 1]);
all_MAT = all_MAT .* modMAT;

% Add some gaussian smoothing to the ground truth matrix'
all_MAT = imgaussfilt(all_MAT, sigma);

% Gaussian has changed range of values. Rescale to 0-1.
all_MAT = rescale(all_MAT);

% Add silence to start and end
all_MAT = [zeros(n_freqs_total, silence*fs), all_MAT, zeros(n_freqs_total, silence*fs)];

% Convert magnitude to power (W)
all_MAT = all_MAT.^2;

%% Compute Transforms

% Compute short windowed STFT
[stft_shortwin, stftshort_freq, stftshort_time] = spectrogram(signal, ...
    win2, ceil(win2*(overlap2/100)), n_fft, fs, "yaxis");

% Compute long windowed STFT
[stft_longwin, stftlong_freq, stftlong_time] = spectrogram(signal, ...
    win1, ceil(win1*(overlap1/100)), n_fft, fs, "yaxis");

% Compute CWT
[cwlet, f_cwlet] = cwt(signal, fs, VoicesPerOctave=vpo, WaveletParameters = [gamma, tbp], FrequencyLimits=[fmin fmax]);

% Compute SLT
slt = nfaslt(signal, fs, [dcfilt, fmax], n_freqs_total, c1, order, mult);


%% Unit Convertsions & Normalizations

% Data to be plotted
% Convert algorithm outputs to power (W)
stft_longwin_plot = (abs(stft_longwin).^2);     % spectrogram() returns complex data.
stft_shortwin_plot = (abs(stft_shortwin).^2);   % spectrogram() returns complex data.
cwlet_plot = (abs(cwlet) .^2);                  % cwt() returns complex data.
slt_plot = (slt .^2);                           % nfaslt returns magnitude.

% If powerscaling is set to plot power in dB, do log conversions:
switch powerscaling
    case 'log'
        % Convert linear power (W) to log power (dBW)
        stft_longwin_plot = 10*log10(stft_longwin_plot / 1);
        stft_shortwin_plot = 10*log10(stft_shortwin_plot / 1);
        cwlet_plot = 10*log10(cwlet_plot / 1);
        slt_plot = 10*log10(slt_plot / 1);

        % Log conversions will result in some -inf values.
        % Replace -inf with the minimum value in each matrix.
        stft_longwin_plot(isinf(stft_longwin_plot)) = min(stft_longwin_plot(~isinf(stft_longwin_plot)), [], 'all');
        stft_shortwin_plot(isinf(stft_shortwin_plot)) = min(stft_shortwin_plot(~isinf(stft_shortwin_plot)), [], 'all');
        cwlet_plot(isinf(cwlet_plot)) = min(cwlet_plot(~isinf(cwlet_plot)), [], 'all');
        slt_plot(isinf(slt_plot)) = min(slt_plot(~isinf(slt_plot)), [], 'all');
    case 'lin'
end

% % Normalize to max=1
stft_shortwin_plot = stft_shortwin_plot ./ max((stft_shortwin_plot), [], 'all');
stft_longwin_plot = stft_longwin_plot ./ max((stft_longwin_plot), [], 'all');
cwlet_plot = cwlet_plot ./ max((cwlet_plot), [], 'all');
slt_plot = slt_plot ./ max((slt_plot), [], 'all');

% Data for use in measures of error
% Convert algorithm outputs to power (W)
stft_longwin = abs(stft_longwin).^2;    % spectrogram returns complex data.
stft_shortwin = abs(stft_shortwin).^2;  % spectrogram returns complex data.
cwlet = abs(cwlet) .^2;                 % cwt returns complex data.
slt = slt .^2;              % nfaslt returns magnitude.

% Normalize to max = 1
stft_shortwin = stft_shortwin ./ max((stft_shortwin), [], 'all');
stft_longwin = stft_longwin ./ max((stft_longwin), [], 'all');
cwlet = cwlet ./ max((cwlet), [], 'all');
slt = slt ./ max((slt), [], 'all');


%% Compute time & frequency resolutions

% STFT Resolutions
short_STFT_fres = (fmax-fmin) / win2;
long_STFT_fres = (fmax-fmin) / win1;
short_STFT_tres = (win2-(ceil(win2*(overlap2/100))))*(1/fs);
long_STFT_tres = (win1-(ceil(win1*(overlap1/100))))*(1/fs);

% CWT Resolutions
idx1 = 1;
idx2 = 2;
cwletfrqs = flip(f_cwlet);
for i = 1:length(cwletfrqs)
    cwlet_fres_vec(i) = cwletfrqs(idx2) - cwletfrqs(idx1);
    if idx2 < length(cwletfrqs)
        idx1 = idx1 +1;
        idx2 = idx2 +1;
    end
end

%% Compute Error

% All matrices must be the same size as groundtruth.
% Scale all to (n_freqs_total x n_samps_total).
stft_shortwin_resz = imresize(stft_shortwin, [n_freqs_total, n_samps_total]);
stft_longwin_resz = imresize(stft_longwin, [n_freqs_total, n_samps_total]);
cwlet_resz= imresize(cwlet, [n_freqs_total, n_samps_total]);
slt_resz = imresize(slt, [n_freqs_total, n_samps_total]);
% The block above is incorrect. imresize() warps the image.
% Have experimented with resize methods, but because cwt has inherently
% non-linear frequency scaling, no linear resize method will work. Need to
% generate a ground truth for each method. Branching in Git for this
% change.

% Compare similarity to ground truth via Root Mean Square Error:
stft_shortwin_freqerror = mean(rmse(stft_shortwin_resz, all_MAT, 1));
stft_shortwin_timeerror = mean(rmse(stft_shortwin_resz, all_MAT, 2));
stft_shortwin_totalerror = rmse(stft_shortwin_resz, all_MAT, 'all');

stft_longwin_freqerror = mean(rmse(stft_longwin_resz, all_MAT, 2));
stft_longwin_timeerror = mean(rmse(stft_longwin_resz, all_MAT, 1));
stft_longwin_totalerror = rmse(stft_longwin_resz, all_MAT, 'all');

cwlet_freqerror = mean(rmse(cwlet_resz, all_MAT, 2));
cwlet_timeerror = mean(rmse(cwlet_resz, all_MAT, 1));
cwlet_totalerror = mean(rmse(cwlet_resz, all_MAT, 'all'));

slt_freqerror = mean(rmse(slt_resz, all_MAT, 2));
slt_timeerror = mean(rmse(slt_resz, all_MAT, 1));
slt_totalerror = rmse(slt_resz, all_MAT, 'all');

% Compare similarity to ground truth via Structural Similarity Index:
stft_shortwin_SSIM = ssim(stft_shortwin_resz, all_MAT);
stft_longwin_SSIM = ssim(stft_longwin_resz, all_MAT);
cwlet_SSIM = ssim(cwlet_resz, all_MAT);
slt_SSIM = ssim(slt_resz, all_MAT);

% Dynamic STFT Names
stftshort_name = ['STFT, ', num2str(win2), 'pt. Window']; %, num2str(overlap2), ' % Overlap'
stftlong_name = ['STFT, ', num2str(win1), 'pt. Window']; % , num2str(overlap1), ' % Overlap'



%% Plot Figure 1 - Time Domain Signal

% Test signal
figure(1)
t1 = tiledlayout(3, 1);
nexttile
plot(t_vec_total, signal_sil)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])
nexttile
plot(t_vec_total, mod_sil)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-.05 1.5])
nexttile
plot(t_vec_total, signal)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('c');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])

set(gcf, 'Position', [50 100 1000 500])
saveas(gcf,'Timedomain_Test_signal','svg')

%% Plot Figure 2 & 3 - Error

% Collate Error data for plotting
xlabels1 = categorical({'Freq Axis', 'Time Axis', 'Total RMSE'});
xlabels1 = reordercats(xlabels1,{'Freq Axis', 'Time Axis', 'Total RMSE'});
ydata1 = [stft_shortwin_freqerror, stft_longwin_freqerror, cwlet_freqerror, slt_freqerror;
    stft_shortwin_timeerror, stft_longwin_timeerror, cwlet_timeerror, slt_timeerror;
    stft_shortwin_totalerror, stft_longwin_totalerror, cwlet_totalerror, slt_totalerror];
xlabels2 = categorical({[num2str(win2), 'pt ', 'STFT'], [num2str(win1), 'pt ', 'STFT'], 'CWT', 'SWT'});
xlabels2 = reordercats(xlabels2,{[num2str(win2), 'pt ', 'STFT'], [num2str(win1), 'pt ', 'STFT'], 'CWT', 'SWT'});
ydata2 = [stft_shortwin_SSIM, stft_longwin_SSIM, cwlet_SSIM, slt_SSIM];

% Plot label rounding
np1 = 3;
np2 = 3;
np3 = 2;
np4 = 4;
np5 = 3;
np6 = 3;
np7 = 3;

% Plot RMSE
figure(2)
b1 = bar(xlabels1, ydata1);
xtips1 = b1(1).XEndPoints;
ytips1 = b1(1).YEndPoints;
labels1 = string(round(b1(1).YData, np1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips2 = b1(2).XEndPoints;
ytips2 = b1(2).YEndPoints;
labels2 = string(round(b1(2).YData, np1));
text(xtips2,ytips2,labels2,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips3 = b1(3).XEndPoints;
ytips3 = b1(3).YEndPoints;
labels3 = string(round(b1(3).YData, np1));
text(xtips3,ytips3,labels3,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips4 = b1(4).XEndPoints;
ytips4 = b1(4).YEndPoints;
labels4 = string(round(b1(4).YData, np1));
text(xtips4,ytips4,labels4,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
ylim([0 0.14])
lg = legend(stftshort_name, stftlong_name, 'CWT', 'SWT');
lg.Location = 'Northwest';
ylabel 'RMSE re. Ground Truth'
grid on
set(gca, FontSize=12, FontName='Calibri')
set(gcf, 'Position', [50 50 1000 300])
saveas(gcf,'Final_methods_analytical_RMSE_TF','svg')

% Plot SSI
figure(3)
b2 = bar(xlabels2, ydata2);
xtips1_2 = b2(1).XEndPoints;
ytips1_2 = b2(1).YEndPoints;
labels1_2 = string(round(b2(1).YData, np2));
text(xtips1_2,ytips1_2,labels1_2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel 'SSI re. Ground Truth'
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylim([0 (round(max(ydata2)*1.3, np2-2))])
grid on
set(gca, FontSize=12, FontName='Calibri')
xtickangle(90)

set(gcf, 'Position', [50 100 350 300])
saveas(gcf,'Final_methods_ERROR_analytical_SSI','svg')


%% Plot More Figures - Time-Freq Representations

% Common Axis limits
freqlim = [10 70];
timelim = [0 (n_samps_total/fs)];

% Init figure
figure (6)
t1 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], FontSize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot STFT with Short Window
nexttile
surf(stftshort_freq, stftshort_time, stft_shortwin_plot', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], FontSize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t1.TileSpacing = 'compact';
t1.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_shortSTFT','svg')

figure (7)
t2 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90],  FontSize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot STFT with Long Window
nexttile
surf(stftlong_freq, stftlong_time, stft_longwin_plot', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t2.TileSpacing = 'compact';
t2.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_longSTFT','svg')

figure (8)
t3 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot CWT
nexttile
surf(f_cwlet, t_vec_total, cwlet_plot', EdgeColor="none", FaceColor="texturemap")
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t3.TileSpacing = 'compact';
t3.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_CWT','svg')


figure (9)
t4 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot Superlets
nexttile
surf(f_vec_total, t_vec_total, slt_plot', EdgeColor="none", FaceColor="texturemap")
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t4.TileSpacing = 'compact';
t4.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_SWT','svg')


t2.TileSpacing = 'compact';
t2.Padding = 'loose';