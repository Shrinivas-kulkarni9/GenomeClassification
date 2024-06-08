% This function takes as input the magnitude spectra and returns
% the correlation coefficient matrix
function mat = PCC(spectra)
    mat = corrcoef(spectra);
    mat = (1-mat)/2;
    mat = (mat+mat')/2;
end