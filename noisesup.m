function OUT = noisesup(IN, threshold)

INmin = min(IN, [], 'all');
INmax = max(IN, [], 'all');
r = INmax - INmin;
thresh_scaled = r * threshold;
OUT = IN;
OUT(OUT < thresh_scaled) = INmin;

% % debug
% figure
% tiledlayout('flow')
% nexttile
% imagesc(IN)
% nexttile
% imagesc(OUT)
% pause

end