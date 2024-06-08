% This function takes a matrix as input and return the MAGNITUDE
% SPECTRA of each of its columns
function specMat = magSpec(mat)
    specMat = fft(mat);
    specMat = abs(specMat);
end