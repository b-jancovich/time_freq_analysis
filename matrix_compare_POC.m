clear

% Make a matrix representing the true signal
t = 1:1:10;
mtruth = flip(diag(rescale(square(t, 70))));

% Make a matrix representing a flawed time-frequency transformation of the
% true signal
m1 = flip([0 1 1 0 0 0 0 0 0 0
      0 0 1 1 0 0 0 0 0 0 
      0 0 0 1 1 0 0 0 0 0
      0 0 0 0 1 1 0 0 0 0 
      0 0 0 0 0 1 1 0 0 0 
      0 0 0 0 0 0 1 1 0 0 
      0 0 0 0 0 0 0 1 1 0
      0 0 0 0 0 0 0 0 1 1 
      0 0 0 0 0 0 0 0 0 0
      0 0 0 0 0 0 0 0 0 0]);

% Make another two, similarly flawed time-freq transformations
m2 = m1 + circshift(m1, [0, -1]);
m2 = m2 + circshift(m2, [-2, -1]);
m2(m2>1) = 1;
m3 = flip(eye(10, 10));

% Compare their similarity to ground truth via Root Mean Square Error:
m1_error = rmse(m1, mtruth, 'all');
m2_error = rmse(m2, mtruth, 'all');
m3_error = rmse(m3, mtruth, 'all');

% Plot
t = tiledlayout ('flow');
nexttile
surf(mtruth, EdgeColor="none")
set(gca, XDir="reverse", View=[90 90])
title('Ground Truth')
xlim([1 10])
ylim([1 10])
nexttile
surf(m1, EdgeColor="none")
set(gca, XDir="reverse", View=[90 90])
title(strcat('Time-Freq Algorithm 1', ', RMSE from Ground Truth= ', num2str(m1_error)))
xlim([1 10])
ylim([1 10])
nexttile
surf(m2, EdgeColor="none")
set(gca, XDir="reverse", View=[90 90])
title(strcat('Time-Freq Algorithm 2', ', RMSE from Ground Truth= ', num2str(m2_error)))
xlim([1 10])
ylim([1 10])
nexttile
surf(m3, EdgeColor="none")
set(gca, XDir="reverse", View=[90 90])
title(strcat('Time-Freq Algorithm 3', ', RMSE from Ground Truth= ', num2str(m3_error)))
xlim([1 10])
ylim([1 10])
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf, 'Position', [800 300 1000 600])
saveas(gcf,'time-freq_method_compare','svg')

