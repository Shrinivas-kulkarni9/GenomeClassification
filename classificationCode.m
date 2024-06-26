function [ accuracy, avg_acc, clNames ] = classificationCode( disMat,alabels, folds, totalSeq )

    cv = cvpartition(totalSeq,'kfold',folds);
    cMat1=cell(1,folds);
    cMat2=cell(1,folds);
    cMat3=cell(1,folds);
    cMat4=cell(1,folds);
    cMat5=cell(1,folds);
    cMat6=cell(1,folds);
    cn=unique(alabels);
    n = length(cn);
    ord = [];
    for xx=1:n
        ord=[ord xx];
    end

    parfor i = 1:folds
        AllTestInd = test(cv,i);
        testInd = find(AllTestInd);
        AllTrainInd = training(cv,i);
        trainInd = find(AllTrainInd);
        trainSet = disMat(trainInd,trainInd);
        testSet = disMat(testInd,trainInd);

        c1 = fitcdiscr(...
        trainSet, ...
        alabels(trainInd), ...
        'DiscrimType', 'linear', ...
        'Gamma', 0, ...
        'FillCoeffs', 'off', ...
        'ClassNames', cn);
        plabel1 = predict(c1,testSet);
        cMat1{i} = confusionmat(alabels(testInd),plabel1,'Order',ord);
        
        if(n==2)
            c2 = fitcsvm(...
            trainSet, ...
            alabels(trainInd), ...
            'KernelFunction', 'linear', ...
            'PolynomialOrder', [], ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true, ...
            'ClassNames', cn);
            plabel2 = predict(c2,testSet);
            cMat2{i} = confusionmat(alabels(testInd),plabel2,'Order',ord);
        else
            template = templateSVM(...
            'KernelFunction', 'linear', ...
            'PolynomialOrder', [], ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
            c2 = fitcecoc(...
            trainSet, ...
            alabels(trainInd), ...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', cn);
            plabel2 = predict(c2,testSet);
            cMat2{i} = confusionmat(alabels(testInd),plabel2,'Order',ord);
        end

        if(cn==2)
            c3 = fitcsvm(...
            trainSet, ...
            alabels(trainInd), ...
            'KernelFunction', 'polynomial', ...
            'PolynomialOrder', 2, ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true, ...
            'ClassNames', cn);
            plabel3 = predict(c3,testSet);
            cMat3{i} = confusionmat(alabels(testInd),plabel3,'Order',ord);
        else
             template = templateSVM(...
            'KernelFunction', 'polynomial', ...
            'PolynomialOrder', 2, ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
            c3 = fitcecoc(...
            trainSet, ...
            alabels(trainInd), ...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', cn);
            plabel3 = predict(c3,testSet); 
            cMat3{i} = confusionmat(alabels(testInd),plabel3,'Order',ord);            
        end
    
        c4 = fitcknn(...
        trainSet, ...
        alabels(trainInd), ...
        'Distance', 'Euclidean', ...
        'Exponent', [], ...
        'NumNeighbors', 1, ...
        'DistanceWeight', 'Equal', ...
        'Standardize', true, ...
        'ClassNames', cn);
        plabel4 = predict(c4,testSet);
        cMat4{i} = confusionmat(alabels(testInd),plabel4,'Order',ord); 

        if(totalSeq<=2000)
            subspaceDimension = max(1, min(74, length(trainSet) - 1));
            c5 = fitcensemble(...
            trainSet, ...
            alabels(trainInd), ...
            'Method', 'Subspace', ...
            'NumLearningCycles', 30, ...
            'Learners', 'discriminant', ...
            'NPredToSample', subspaceDimension, ...
            'ClassNames', cn);
            plabel5 = predict(c5,testSet);
            cMat5{i} = confusionmat(alabels(testInd),plabel5,'Order',ord); 

            subspaceDimension = max(1, min(74, length(trainSet) - 1));
            c6 = fitcensemble(...
            trainSet, ...
            alabels(trainInd), ...
            'Method', 'Subspace', ...
            'NumLearningCycles', 30, ...
            'Learners', 'knn', ...
            'NPredToSample', subspaceDimension, ...
            'ClassNames', cn);
            plabel6 = predict(c6,testSet);
            cMat6{i} = confusionmat(alabels(testInd),plabel6,'Order',ord); 
        end
    end

    cMat{1}=cMat1;
    cMat{2}=cMat2;
    cMat{3}=cMat3;
    cMat{4}=cMat4;
    if(totalSeq<=2000)
        cMat{5}=cMat5;
        cMat{6}=cMat6;
        clNames = {"LinearDiscriminant","LinearSVM","QuadraticSVM","FineKNN","SubspaceDiscriminant","SubspaceKNN","AverageAccuracy"};
    else
        clNames = {"LinearDiscriminant","LinearSVM","QuadraticSVM","FineKNN",'AverageAccuracy'};
    end
    accuracy=cell(1,length(cMat));
   
    for i=1:length(accuracy)
        cm = cMat{i};
        cMatrix = 0;
        for k = 1:folds
         cMatrix = cMatrix+cm{k};
        end
        accuracy{i} = round((trace(cMatrix)/totalSeq)*100,1);

    end    
    avg_acc = round(sum([accuracy{:}])/length(accuracy),1);
    
end

