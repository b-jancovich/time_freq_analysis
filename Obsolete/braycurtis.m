function S = braycurtis(A, B)

D = abs(A - B) ./ (A + B);
D(A == 0 | B == 0) = nan;
S = 1 - D / max(D(:));
S(isnan(S)) = 0;
S = mean(S, 'all');

end