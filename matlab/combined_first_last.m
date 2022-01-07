%%
% Quick analysis of firing rate data of neurons from behaving animals
%

close all
clear all;

loadColors;

%% first days
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg_firstdays.mat');
% dca1=X_full_avg;
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg_firstdays.mat');
% vca1=X_full_avg;

%% last days
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg_lastdays.mat');
dca1=X_full_avg;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg_lastdays.mat');
vca1=X_full_avg;

% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg.mat');
% dca1=X_full_avg;
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg.mat');
% vca1=X_full_avg;

combinedCA1=[dca1;vca1];
times = -2:0.2:((31.*0.2)-2.2);
times(:,11) = [];

%% export data matrix for pca
% mstore=zscore(mres_all')';
% save('data.mat','mstore');
if(0)
%% PCA
[corr_comps_mall, corr_score, latent,tsquared,corr_explained_mall] = pca(zscore(mres_all')'); %Across trial average
disp(mean(mres_all'));
keyboard
pcomp = 3;
cs = copper(pcomp);
set(groot,'defaultAxesColorOrder',cs)

figure('Position', [1317        -331         953         901], 'Name', [cell2mat(id) ' rule: ' num2str(rule_switch)]);

% ht = [];
% subplot(2,2,2)
% for i=1:pcomp
%     h = plot(times, corr_comps_mall(:,i)./max(corr_comps_mall(:,i)), 'Linewidth', 2);
%     hold on;
% end
% title(['correct PCA across trial mean (' num2str(sum(corr_explained_mall(1:pcomp)),3) '% exp. var.)']);

% plot([0 0], [-1 1], ':k', 'Linewidth', 1);
% box off;
% ylim([-0.4 0.5]);
% ylabel('norm. activity [-]');

ht = [];
% subplot(2,2,1)
for i=1:pcomp
    h = plot(times, corr_comps_mall(:,i), 'Linewidth', 2);
    ht = [ht h];
    hold on;
end
title(['PCA across trial averaged data (' num2str(sum(corr_explained_mall(1:pcomp)),3) '% exp. var.)']);

legend(ht, {['1st comp. (' num2str(corr_explained_mall(1),2) '%)'], ['2nd comp. (' num2str(corr_explained_mall(2),2) '%)'], ['3rd comp. (' num2str(corr_explained_mall(3),2) '%)']})
plot([0 0], [-0.4 0.4], ':k', 'Linewidth', 1);
ylim([-0.4 0.5]);
box off;
ylabel('Norm. firing rate [-]');
xlabel('Time ( in seconds )');

end
% %% saving ventral data
%     save('vCA1_Xfull_Avg','X_full_avg');


if(1)
%% DPCA on decision, strategy, time and interaction of decision and strategy
Xfull=(combinedCA1);

% D = size(trialNum,1);
% minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
% meanFiringRate = mean(reshape(Xfull, D, []), 2);
% n = find(minN >= 2 & meanFiringRate < 50);
% t = 1:length(times);

% combinedParams = {{1, [1,2]},{2, [2,3]},{3, [3,4]}};
combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction D/S'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
times1 = -2:0.2:((31.*0.2)-2.2);
timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 12,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

explVar = dpca_explainedVariance(Xfull, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(Xfull, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', times,...
    'legendSubplot', 16,...
    'whichMarg', whichMarg, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours,...
    'timeEvents',timeEvents);
end
if(0)
%% DPCA on day, decision, strategy, time and interaction between decision and strategy
Xfull=(X_full_daywise);

%combinedParams = {{1, [1,2]},{2, [2,3]},{3, [3,4]}};
combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
%combinedParams ={{1}}; %gives blank why?
margNames = {'Day','Decision', 'Time', 'Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
decodingClasses = {[1 1; 2 2; 3 3; 4 4; 5 5; 6 6], [1 2; 1 2; 1 2; 1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6; 7 8; 9 10; 11 12]};

T=5;
time = (1:T) / 10;
timeEvents = time(round(length(time)/2));
keyboard
[W,V,whichMarg] = dpca(Xfull, 12,...
'combinedParams',combinedParams,'lambda', 0.0025); %2.919292602539062e-05

explVar = dpca_explainedVariance(Xfull, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(Xfull, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', times,...
    'legendSubplot', 16,...
    'whichMarg', whichMarg, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);
end
if(0)
%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(Xfull, @dpca_plot_default, ...
   'combinedParams', combinedParams);
end
if(0)
%% Step 1: PCA of the dataset
Xfull=(X_full_avg);

combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Decision', 'Strategy', 'Time', 'Interaction D/S'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
timeEvents = times1(11);

X = Xfull(:,:);
% X = bsxfun(@minus, X, mean(X,2));
keyboard
[W,~,~] = svd(X, 'econ');
W = W(:,1:12);

% minimal plotting
dpca_plot(Xfull, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(Xfull, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(Xfull, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', times,                        ...
    'marginalizationNames', margNames, ...
    'timeEvents', timeEvents,               ...
    'marginalizationColours', margColours);     

end
if(0)
%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).
Xfull=(X_full_avg);
combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Decision','Strategy', 'Time', 'Interaction D/S'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

ifSimultaneousRecording = false;
optimalLambda = dpca_optimizeLambda(Xfull, XFiringRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ... %     'numComps', 8,...
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'display', 'yes',...
    'filename', 'tmp_optimalLambdas.mat');
disp(optimalLambda);
Cnoise = dpca_getNoiseCovariance(Xfull, ...
    XFiringRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(Xfull, 12, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(Xfull, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(Xfull, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', times,                        ...
    'legendSubplot', 16);
end
if(0)
%% Bar graph of Optimal Lambda's

c = categorical({'All data','Sc08','Sc03','Sc11','JP28'});
optimalLamdas = [4.378938903808594e-05 1.477891880035400e-04 6.568408355712890e-05 3.325256730079651e-04 0.0025];
figure;
bar(c,optimalLamdas);
end

if(0)
%% Bar plot of firing rates of all first days for all animals
    figure;
    hist(corr_mres_all', 30);
end


if(0)
%% Optional: estimating "signal variance"
ifSimultaneousRecording = false;
Cnoise = dpca_getNoiseCovariance(Xfull, ...
    XFiringRates, trialNum, 'simultaneous', ifSimultaneousRecording);

explVar = dpca_explainedVariance(Xfull, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(Xfull, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', times,                        ...
    'legendSubplot', 16);
end
if(0)
    ifSimultaneousRecording=false;
    decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]};


    [accuracy,brier] = dpca_classificationAccuracy(Xfull, XFiringRates, trialNum, ...
    'lambda', 2.919292602539062e-05, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 100, ...        % increase to 100
    'verbose', 'yes',...
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], [], [], decodingClasses,'time', times,'timeEvents',timeEvents)

accuracyShuffle = dpca_classificationShuffled(XFiringRates, trialNum, ...
    'lambda', 2.919292602539062e-05, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 100, ...        % increase to 100
    'numShuffles', 100, ...  % increase to 100 (takes a lot of time)
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses,'time', times,'timeEvents',timeEvents);
end