clear all


load('/mnt/wwn-0x5000c500cc87eb78/dPCA/python/hidds_dni_all_epochs.npy_Xfull0.mat');

X = Xfull0(:, :);
X = bsxfun(@minus, X, mean(X,2));
% 

% 
% [W,~,~] = svd(X, 'econ');
% W = W(:,1:20);
% 
% % minimal plotting
% % dpca_plot(Xfull0, W, W, @dpca_plot_default);
% 
% 
% % computing explained variance
% explVar = dpca_explainedVariance(Xfull0, W, W, ...
%     'combinedParams', combinedParams);
% 
% % a bit more informative plotting
% dpca_plot(Xfull0, W, W, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'time', time,                        ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours);



% dpca
combinedParams = {{1, [1 2]}, {2}};
margNames = {'Cues', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
time=(1:10);

[W,V,whichMarg] = dpca(X, 5, ...
    'combinedParams', combinedParams);


explVar = dpca_explainedVariance(X, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(X, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);