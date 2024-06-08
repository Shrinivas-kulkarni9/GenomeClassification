close all;
clc ;
set(0, 'DefaultFigureRenderer', 'painters');
% can be uncommented for MATLAB2017 or earlier versions
%digits(4);


fprintf('Reading sequences .... \n');
[feature_mtx,labels,num_labels] = readData("DataBase/3classes/","Real");
fprintf('Generating numerical sequences, applying DFT, computing magnitude spectra .... \n');
siz = size(feature_mtx);
totalSeq = siz(2);
DFTmag = magSpec(feature_mtx);

%distance calculation by Pearson correlation coefficient
% change 'cor' to 'euc' for Euclidean
fprintf('Computing Distance matrix .... \n');
disMat = PCC(DFTmag);

%Multi-dimensional Scaling
fprintf('Performing Multi-dimensional scaling .... \n');
[Y,eigvals] = cmdscale(disMat);

%3D  plot
fprintf('Generating 3D plot .... \n');
index=1;
counter=1;
Cluster = zeros(1,totalSeq);
for i=1:totalSeq   
    Cluster(i)=index;
    if(counter==pointsPerCluster{index})
        index=index+1;
        counter=0;
    end
    counter= counter+1;
end
uniqueClusters  = unique(Cluster);
cmap = distinguishable_colors(numberOfClusters);
hf = figure;
hold on;
for h=1:numberOfClusters
    cIndex = Cluster == uniqueClusters(h);
    plot3(Y(cIndex,1),Y(cIndex,2),Y(cIndex,3),'.','markersize', 15, 'Color',cmap(h,:),'DisplayName',clusterNames{h});
end
view(3), axis vis3d, box on, datacursormode on
xlabel('x'), ylabel('y'), zlabel('z')
tname = strcat(selectedFolder,' (',int2str(totalSeq),' Sequences',')');
title(tname)
hdt = datacursormode(hf);
set(hdt,'UpdateFcn',{@myupdatefcn,Y,AcNmb})
legend('show');


%Phylogenetic Tree
% plotted for datasets with atmost 50 sequences
% parameter can be changed in following statement 
if(totalSeq<=50)
    fprintf('Creating Phylogenetic Tree .... \n');
    UPGMAtree = seqlinkage(disMat,'UPGMA',AcNmb);
    plot(UPGMAtree);
end

%Classification Code
clear a;
a=[];
for i=1:numberOfClusters
    for j=1:pointsPerCluster{i}
        a=[a; i];
    end
end
ATestlg = [disMat a];
rng(15,'twister');

alabels = a;
fprintf('Performing classification .... \n');
folds=10;
if (totalSeq<folds)
    folds = totalSeq;
end
[accuracy, avg_accuracy, clNames] = classificationCode(disMat,alabels, folds, totalSeq);
acc = [accuracy avg_accuracy];
s.ClassifierModel=cellstr(clNames.');
s.Accuracy=cell2mat(acc).';
ClassificationAccuracyScores = struct2table(s);

%new sequence classification example
%uncomment the following line and replace the accesion number 
%newSeqClassify('NC_028718.1', mLen, disMat, lg, alabels, clusterNames)

fprintf('**** Processing completed ****\n');