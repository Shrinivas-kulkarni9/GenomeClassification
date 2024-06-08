function [] = newSeqClassify(testSeqAcNmb, mLen, disMat, lg, alabels, clusterNames)
    pLabel='';
    testS = getgenbank(testSeqAcNmb,'FileFormat', 'fasta');
    testSeq = testS.Sequence;
    testSeq = regexprep(testSeq,'[^A,^C, ^G, ^T]','','ignorecase');
    ns = numMappingPP(testSeq);
    nsLen = length(ns);
    I = mLen-nsLen;
    if(I>0)
        nsTemp = wextend('1','asym',ns,I);
        nsNew = nsTemp((I+1):length(nsTemp));
    elseif(I<0)
        nsNew=ns(1:mLen);
    else
        nsNew = ns;
    end
    fn = fft(nsNew);
    lgn = abs(fn); 
    testV = [];
    for d=1:length(lg)
        aa=(1-corrcoef(lgn,lg{d}))/2;
        testV = [testV aa(1,2)];
    end
 
    cn = unique(alabels);
        
    template = templateSVM(...
            'KernelFunction', 'polynomial', ...
            'PolynomialOrder', 2, ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
    cModel = fitcecoc(...
            disMat, ...
            alabels, ...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', cn); 
    
    clabel = predict(cModel,testV);    
    pLabel = clusterNames{clabel};
    fprintf('Perdicted Label: %s \n', pLabel);
end

