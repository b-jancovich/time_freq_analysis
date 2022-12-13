function score = binaryerrdetect(I, REF, threshold)

I(I<threshold) = 0;
I(I>=threshold) = 1;
REF(REF<threshold) = 0;
REF(REF>=threshold) = 1;

score = (sum(I, 'all') - sum(REF, 'all')) ./ (size(I, 1) * size(I, 2));
end
