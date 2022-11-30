function similarity = jacsim(IN, REF, tolerance)

IN(IN>tolerance) = 1;
REF(REF>tolerance) = 1;
IN(IN<=tolerance) = 0;
REF(REF<=tolerance) = 0;

similarity = jaccard(IN, REF);

% % Debug image 
% figure
% imshowpair(IN, REF)
% pause
% close

end
