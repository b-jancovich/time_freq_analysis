%% Comparison of Time-Frequency Analysis Algorithms with Ground Truth.
%
% This script generates a carrier signal (linearly frequency swept sine)
% and an amplitude modulation signal (linearly frequency swept square)
% and multiplies them together, creating an amplitude modulated test signal.
% It then constructs a matrix that can be thought of as a synthetic
% spectrogram representing the "Ground Truth" of the test signal's
% power, as a function of both time & frequency. This spectrogram is
% constructed analytically using knowledge of the test signal's TF
% charactertistics, rather than by transforming or otherwise analysing
% the signal, so is free from the artifacts, inaccuracies and loss of 
% detail associated with TF analysis. This ground truth matrix therfore 
% serves as a theoretically perfect spectrogram of the test signal. 
% Several transform based TF analyses are then performed and compared 
% against the ground truth using statistical tests of error.
%
% Methods compared:
% Short Time Fourier Transform (short window)
% Short Time Fourier Transform (long window)
% Continuous Wavelet Transform
% Fractional Adaptive Superresolution Wavelet Transform
%
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
% NOTE: fc1 and fc2 must NOT be the same frequency.

% Test Signal Parameters:
fc1 = 50;                   % Sine Sweep start frequency (Hz) Must be ~= fc2
fc2 = 30;                   % Sine Sweep end frequency (Hz) Must be ~= fc1
fam1 = 2;                   % Amplitude Modulation sweep start frequency (Hz)
fam2 = 7;                  % Amplitude Modulation sweep end frequency (Hz)
duration_sweep = 5;         % Sweep Duration
duration_silence = 1;       % Silence to pad start and end (seconds)
fs = 250;                   % Sampling Frequency (Hz)
phi_carrier = 0;            % Initial phase of carrier waveform (degrees)
phi_am = -90;               % Initial phase of AM waveform (degrees)
                            % chirp() returns cosine, so -90 ensures sweep
                            % begins at amplitude = 1 and holds for one 
                            % half cycle.

% Ground Truth Parameters:
sigma = 0.3;                % Standard deviation of gaussian filter

% Plotting Parameters
intensity_units = 'mag';    % Color axis units for TFR plotting. 
%                           Options are:
%                           Linear scaled magnitide = 'mag' 
%                           Log scaled magnitude (dB) = 'magdb' 
%                           Linear scaled power = 'pow' 
%                           Log scaled power (dBW) = 'powdb'

% Image Comparison Parameters:
resize_method = 'nearest';  % Interpolation method - Nearest minimses artifacts

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest
fmax = fs/2;                % Highest frequency of interest
f_res = 0.2;                % Frequency resolution (Hz)

% STFT Window Sizes
n_fft = 2*((fs/2)/f_res);   % spectrogram FFT length (samples)
win_long = 250;             % Window length for long STFT (samples)
win_short = 50;             % Window length for short STFT (samples)
overlap_long = 75;          % Window overlap % for long STFT
overlap_short = 75;         % Window Overlap % for short STFT

% CWT Parameters
time_bandwidth = 120;        % Time bandwidth product of Morse wavelet.
vpo = 48;                   % Voices per Octave. Must be in range 10 : 48

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 40];    % Interval of superresolution orders
mult = 1;           % Multiplicative ('1') or additive ('0') superresolution

%% Generate time & Frequency Vectors

% Frequency & samp counts
n_samps_sweep = duration_sweep * fs;
n_samps_total = n_samps_sweep + (2*(duration_silence*fs));
n_freqs_sweep = (fc1-fc2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vectors
% t_vec_sweep = linspace(0, n_samps_sweep/fs, n_samps_sweep);
t_vec_sweep = linspace(0, duration_sweep, n_samps_sweep);
t_vec_total = linspace(0, n_samps_total/fs, n_samps_total);

% Frequency vectors
f_vec_sweep = linspace(fc1, fc2, n_freqs_sweep)';
f_vec_total = linspace(fmin, fmax, n_freqs_total)';

%% Generate signals

% Generate Sine Sweep (Carrier Signal)
sig_c = chirp(t_vec_sweep, fc1, t_vec_sweep(end), fc2, 'linear', phi_carrier); %, phi_carrier

% Generate Square Sweep (Amplitude Modulator)
sig_am = rescale(sign(chirp(t_vec_sweep, fam1, t_vec_sweep(end), fam2, 'linear', phi_am)));

% Pad signal and mod with zeros to insert silence at start & end.
sig_c_sil = [zeros(1, duration_silence * fs), sig_c, zeros(1, duration_silence * fs)];
sig_am_sil = [zeros(1, duration_silence * fs), sig_am, zeros(1, duration_silence * fs)];

% Amplitude Modulation
signal = sig_c_sil .* sig_am_sil;

%% Compute Transforms

% Compute short windowed STFT
[stft_shortwin, stft_shortwin_f, stft_shortwin_t] = spectrogram(signal, ...
    win_short, ceil(win_short*(overlap_short/100)), n_fft, fs, "yaxis");

% Compute long windowed STFT
[stft_longwin, stft_longwin_f, stft_longwin_t] = spectrogram(signal, ...
    win_long, ceil(win_long*(overlap_long/100)), n_fft, fs, "yaxis");

% Compute CWT
[cwlet, cwlet_f] = cwt(signal, fs, VoicesPerOctave=vpo, timeBandWidth=time_bandwidth, FrequencyLimits=[fmin fmax]);

% Compute SLT
slt = nfaslt(signal, fs, [fmin, fmax], n_freqs_total, c1, order, mult);

%% Construct Ground Truth Matrices
% Because each analysis algorithm will return a matrix with a different
% number of rows (frequencies) and columns (times), and free-scale resizing them 
% will corrupt the results, each analysis algorithm must be compared with 
% its own groundtruth matrix having the same aspect ratio.

% Generate groundtruth TFR that is used for plotting only
[groundtruth_t, groundtruth_f, groundtruth] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    f_vec_total, t_vec_total, sigma, duration_sweep,...
    duration_silence, phi_am, fs, 0.5);

% Generate stft_shortwin groundtruth & Corresponding time and freq vectors:
[stft_shortwin_GT_t, stft_shortwin_GT_f, stft_shortwin_GT] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    stft_shortwin_f, stft_shortwin_t, sigma, duration_sweep,...
    duration_silence, phi_am, fs, f_res);

% Generate stft_longwin groundtruth & Corresponding time and freq vectors:
[stft_longwin_GT_t, stft_longwin_GT_f, stft_longwin_GT] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    stft_longwin_f, stft_longwin_t, sigma, duration_sweep,...
    duration_silence, phi_am, fs, f_res);

% Generate CWT groundtruth & Corresponding time and freq vectors:
[cwlet_GT_t, cwlet_GT_f, cwlet_GT] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    cwlet_f, t_vec_total, sigma, duration_sweep,...
    duration_silence, phi_am, fs, f_res, 'log', 'reverse');
% Note: the additional (optional) input arguments here are for logarithmic 
% frequency scaling and a high to low (reversed) frequency axis for CWT. 
% These are to match the GT to the output of the CWT algorithm.

% % Generate SLT groundtruth & Corresponding time and freq vectors:
[slt_GT_t, slt_GT_f, slt_GT] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    f_vec_total, t_vec_total, sigma, duration_sweep,...
    duration_silence, phi_am, fs, f_res);

%% Intensity Unit Conversions & Normalization

% Convert algorithm outputs real magnitude
stft_shortwin = abs(stft_shortwin);  % spectrogram returns complex data. Take absolute value to get magnitude.
stft_longwin = abs(stft_longwin);    % spectrogram returns complex data. Take absolute value to get magnitude.
cwlet = abs(cwlet);                  % cwt returns complex data. Take absolute value to get magnitude.
slt = sqrt(slt);                     % nfaslt returns squared magnitude (power).  Take square root to get magnitude.

switch intensity_units
    case 'mag'
        % do nothing
    case 'magdb'
        stft_shortwin = 20*log10(stft_shortwin); 
        stft_shortwin_GT = 20*log10(stft_shortwin_GT); 
        stft_longwin = 20*log10(stft_longwin);    
        stft_longwin_GT = 20*log10(stft_longwin_GT);    
        cwlet = 20*log10(cwlet);                  
        cwlet_GT = 20*log10(cwlet_GT);                  
        slt = 20*log10(slt);
        slt_GT = 20*log10(slt_GT);
    case 'pow'
        stft_shortwin = stft_shortwin .^2; 
        stft_shortwin_GT = stft_shortwin_GT .^2; 
        stft_longwin = stft_longwin .^2;    
        stft_longwin_GT = stft_longwin_GT .^2;    
        cwlet = cwlet .^2;               
        cwlet_GT = cwlet_GT .^2;               
        slt = slt .^2;
        slt_GT = slt_GT .^2;
    case 'powdb'
        stft_shortwin = 10*log10(stft_shortwin .^2); 
        stft_shortwin_GT = 10*log10(stft_shortwin_GT .^2); 
        stft_longwin = 10*log10(stft_longwin .^2);    
        stft_longwin_GT = 10*log10(stft_longwin_GT .^2);    
        cwlet = 10*log10(cwlet .^2);                  
        cwlet_GT = 10*log10(cwlet_GT .^2);                  
        slt = 10*log10(slt .^2);
        slt_GT = 10*log10(slt_GT .^2);
end

% Normalize to max = 1
stft_shortwin = stft_shortwin ./ max((stft_shortwin), [], 'all');
stft_shortwin_GT = stft_shortwin_GT ./ max((stft_shortwin_GT), [], 'all');
stft_longwin = stft_longwin ./ max((stft_longwin), [], 'all');
stft_longwin_GT = stft_longwin_GT ./ max((stft_longwin_GT), [], 'all');
cwlet = cwlet ./ max((cwlet), [], 'all');
cwlet_GT = cwlet_GT ./ max((cwlet_GT), [], 'all');
slt = slt ./ max((slt), [], 'all');
slt_GT = slt_GT ./ max((slt_GT), [], 'all');

%% Algorithmic TFR Scaling

% Resample the stft_shortwin TRF to match its groundtruth size
stft_shortwin_resz = imresize(stft_shortwin, size(stft_shortwin_GT), Method=resize_method);

% % Resample the stft_shortwin TRF to match its groundtruth size
stft_longwin_resz = imresize(stft_longwin, size(stft_longwin_GT), Method=resize_method);
% 
% % Resample the CWT TRF to match its groundtruth size
cwlet_resz = imresize(cwlet, size(cwlet_GT), Method=resize_method);
% 
% % Resample the SLT TRF to match its groundtruth size
slt_resz = imresize(slt, size(slt_GT), Method=resize_method);

%% TESTING: Plot things to check for errors before continuing to Error Calculation

% Plot the Resized TFRs that will be used for error calculation,
% side-by-side with their corresponding grounstruth TFRs to make sure they
% are appropriately scaled, and there are no erronious frequency or time
% shifts:
figure(1)
tiledlayout(1,2)
nexttile
imagesc(abs(cwlet_resz))
title('cwt')
ylabel('Row Number')
nexttile
imagesc(cwlet_GT)
title('cwt GT')

tiledlayout(4,2)
nexttile
imagesc(abs(stft_shortwin_resz))
title('stft short')
ylabel('Row Number')
nexttile
imagesc(stft_shortwin_GT)
title('stft short GT')
nexttile
imagesc(abs(stft_longwin_resz))
title('stft long')
ylabel('Row Number')
nexttile
imagesc(stft_longwin_GT)
title('stft long GT')
nexttile
imagesc(abs(cwlet_resz))
title('cwt')
ylabel('Row Number')
nexttile
imagesc(cwlet_GT)
title('cwt GT')
nexttile
imagesc(slt_resz)
title('slt')
ylabel('Row Number')
xlabel('Column Number')
nexttile
imagesc(slt_GT)
title('slt GT')
xlabel('Column Number')
sgtitle(['Resized TFRs & Corresponding Groundtruths', newline,...
    'These Matrices Are Used to Compute Measures of Error'],...
    FontWeight='bold')

% Plot the resized TFR's used for error calculation, side-by-side with the 
% originals to ensure resizing hasn't caused any interpolation artifacts.
figure(10)
tiledlayout(1,2)
nexttile
imagesc(abs(stft_shortwin_resz))
title('stft shortr resized')
ylabel('Row Number')
xlabel('Column Number')
nexttile
imagesc(abs(stft_shortwin))
title('stft short')
xlabel('Column Number')

figure(11)
tiledlayout(1,2)
nexttile
imagesc(abs(stft_longwin_resz))
title('stft long resized')
ylabel('Row Number')
xlabel('Column Number')
nexttile
imagesc(abs(stft_longwin))
title('stft long')
xlabel('Column Number')

figure(12)
tiledlayout(1,2)
nexttile
imagesc(abs(cwlet_resz))
title('cwt resized')
ylabel('Row Number')
xlabel('Column Number')
nexttile
imagesc(abs(cwlet))
title('cwt')
xlabel('Column Number')

figure(13)
tiledlayout(1,2)
nexttile
imagesc(slt_resz)
title('slt resized')
ylabel('Row Number')
xlabel('Column Number')
nexttile
imagesc(slt)
title('slt')
xlabel('Column Number')

%% Compute Error

% Compare similarity to ground truth via Root Mean Square Error:
stft_shortwin_freqerror = mean(rmse(stft_shortwin_resz, stft_shortwin_GT, 1));
stft_shortwin_timeerror = mean(rmse(stft_shortwin_resz, stft_shortwin_GT, 2));
stft_shortwin_totalerror = rmse(stft_shortwin_resz, stft_shortwin_GT, 'all');

stft_longwin_freqerror = mean(rmse(stft_longwin_resz, stft_longwin_GT, 2));
stft_longwin_timeerror = mean(rmse(stft_longwin_resz, stft_longwin_GT, 1));
stft_longwin_totalerror = rmse(stft_longwin_resz, stft_longwin_GT, 'all');

cwlet_freqerror = mean(rmse(cwlet_resz, cwlet_GT, 2));
cwlet_timeerror = mean(rmse(cwlet_resz, cwlet_GT, 1));
cwlet_totalerror = mean(rmse(cwlet_resz, cwlet_GT, 'all'));

slt_freqerror = mean(rmse(slt_resz, slt_GT, 2));
slt_timeerror = mean(rmse(slt_resz, slt_GT, 1));
slt_totalerror = rmse(slt_resz, slt_GT, 'all');

% Gather RMSE scores
rmse_freq = [stft_shortwin_freqerror, stft_longwin_freqerror,...
    cwlet_freqerror, slt_freqerror];
rmse_time = [stft_shortwin_timeerror, stft_longwin_timeerror,...
    cwlet_timeerror, slt_timeerror];
rmse_total = [stft_shortwin_totalerror, stft_longwin_totalerror,...
    cwlet_totalerror, slt_totalerror];

% Sort scores by error, low to high
rmse_freq = sort(rmse_freq, 'ascend');
rmse_time = sort(rmse_time, 'ascend');
rmse_total = sort(rmse_total, 'ascend');

% Calculate percent difference between best and 2nd best performer
rmse_freq_pdiff = 100*(rmse_freq(2)-rmse_freq(1))./rmse_freq(1);
rmse_time_pdiff = 100*(rmse_time(2)-rmse_time(1))./rmse_time(1);
rmse_total_pdiff = 100*(rmse_total(2)-rmse_total(1))./rmse_total(1);

% Compare similarity to ground truth via Structural Similarity Index:
stft_shortwin_SSIM = ssim(stft_shortwin_resz, stft_shortwin_GT);
stft_longwin_SSIM = ssim(stft_longwin_resz, stft_longwin_GT);
cwlet_SSIM = ssim(cwlet_resz, cwlet_GT);
slt_SSIM = ssim(slt_resz, slt_GT);

% Gather SSI scores
SSI = [stft_shortwin_SSIM, stft_longwin_SSIM,...
    cwlet_SSIM, slt_SSIM];

% Sort scores by similarity, high to low
SSI = sort(SSI, 'descend');

% Determine percent difference between best and 2nd best performer
SSI_pdiff = 100*(SSI(1)-SSI(2))./SSI(2);

% Dynamic STFT Names
stftshort_name = ['STFT, ', num2str(win_short), 'pt. Window']; %, num2str(overlap2), ' % Overlap'
stftlong_name = ['STFT, ', num2str(win_long), 'pt. Window']; % , num2str(overlap1), ' % Overlap'

%% Plot Figure 1 - Time Domain Signal

% Test signal
figure(2)
t1 = tiledlayout(3, 1);
nexttile
plot(t_vec_total, sig_c_sil)
ylabel('Amplitude (arb.)', FontWeight='normal');
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])
nexttile
plot(t_vec_total, sig_am_sil)
ylabel('Amplitude (arb.)', FontWeight='normal');
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-.2 1.2])
nexttile
plot(t_vec_total, signal)
ylabel('Amplitude (arb.)', FontWeight='normal');
ttl = title('c');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])

t1.Padding = "loose"
set(gcf, 'Position', [50 100 1000 500])
saveas(gcf,'Timedomain_Test_signal','svg')

%% Plot Figure 2 & 3 - Error

% Collate Error data for plotting
xlabels1 = categorical({'Freq Axis', 'Time Axis', 'Total RMSE'});
xlabels1 = reordercats(xlabels1,{'Freq Axis', 'Time Axis', 'Total RMSE'});
ydata1 = [stft_shortwin_freqerror, stft_longwin_freqerror, cwlet_freqerror, slt_freqerror;
    stft_shortwin_timeerror, stft_longwin_timeerror, cwlet_timeerror, slt_timeerror;
    stft_shortwin_totalerror, stft_longwin_totalerror, cwlet_totalerror, slt_totalerror];
xlabels2 = categorical({[num2str(win_short), 'pt ', 'STFT'], [num2str(win_long), 'pt ', 'STFT'], 'CWT', 'SWT'});
xlabels2 = reordercats(xlabels2,{[num2str(win_short), 'pt ', 'STFT'], [num2str(win_long), 'pt ', 'STFT'], 'CWT', 'SWT'});
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
figure(3)
tiledlayout(1, 2)
nexttile
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
ylim([0 0.3])
lg = legend(stftshort_name, stftlong_name, 'CWT', 'SWT');
lg.Location = 'Northwest';
ylabel 'RMSE re. Ground Truth'
grid on
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
set(gca, FontSize=12, FontName='Calibri')
xtickangle(90)

% set(gcf, 'Position', [50 50 1000 300])
% saveas(gcf,'Final_methods_analytical_RMSE_TF','svg')

% Plot SSI
xlabels3 = categorical({[num2str(win_short), 'pt ', newline, 'STFT'], [num2str(win_long), 'pt ', newline, 'STFT'], 'CWT', 'SWT'});
% xlabels3 = reordercats(xlabels2,{[num2str(win_short), 'pt ', newline, 'STFT'], [num2str(win_long), 'pt ', newline, 'STFT'], 'CWT', 'SWT'});
nexttile
% figure(4)
b2 = bar(xlabels3, ydata2);
xtips1_2 = b2(1).XEndPoints;
ytips1_2 = b2(1).YEndPoints;
labels1_2 = string(round(b2(1).YData, np2));
text(xtips1_2,ytips1_2,labels1_2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel 'SSI re. Ground Truth'
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylim([0 (round(max(ydata2)*1.3, np2-2))])
grid on
set(gca, FontSize=12, FontName='Calibri')
xtickangle(90)

% set(gcf, 'Position', [50 100 350 300])
% saveas(gcf,'Final_methods_ERROR_analytical_SSI','svg')

set(gcf, 'Position', [50 50 1000 450])
saveas(gcf,'RSMSE & SSI','svg')

%% Plot More Figures - Time-Freq Representations

% Common Axis limits
freqlim = [10 70];
timelim = [0 7];
colorlabel = 'Magnitude (normalized)';

% Init figure
figure (5)
t1 = tiledlayout(3,2);

% Plot matrix "grund truth" time-frequency representation.
nexttile
TFRplot(groundtruth_t, groundtruth_f, groundtruth, freqlim, timelim)
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
c = colorbar;
c.Label.String = colorlabel;

% leave an empty tile
nexttile

% Plot STFT with Short Window
nexttile
TFRplot(stft_shortwin_t, stft_shortwin_f, stft_shortwin, freqlim, timelim)
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  
c = colorbar;
c.Label.String = colorlabel;

% Plot STFT with Long Window
nexttile
TFRplot(stft_longwin_t, stft_longwin_f, stft_longwin, freqlim, timelim)
ttl = title('c');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  
c = colorbar;
c.Label.String = colorlabel;

% Plot CWT
nexttile
TFRplot(t_vec_total, cwlet_f, cwlet, freqlim, timelim)
ttl = title('d');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
c = colorbar;
c.Label.String = colorlabel;

% Plot Superlets
nexttile
TFRplot(t_vec_total, f_vec_total, slt, freqlim, timelim)
ttl = title('e');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  
c = colorbar;
c.Label.String = colorlabel;

t4.TileSpacing = 'compact';
t4.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 6000])
saveas(gcf,'Groundtruth_vs_algoTFRs','svg')
