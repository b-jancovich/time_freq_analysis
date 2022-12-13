function score = simSURF(I, REF)

I_points = detectSURFFeatures(I);
REF_points = detectSURFFeatures(REF);
[I_features,~] = extractFeatures(I,I_points);
[REF_features,~] = extractFeatures(REF,REF_points);

[~, matchmetric] = matchFeatures(I_features, REF_features, ...
    Method="Approximate", MaxRatio=1, MatchThreshold=10);

score = mean(matchmetric);

end