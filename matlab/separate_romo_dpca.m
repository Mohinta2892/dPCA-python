%%
% Quick analysis of firing rate data of neurons from behaving animals
%

close all
clear all;

% loadColors;

tetrodes = [1:16]; %Tetrodes id used to record neural activity

%basepath = '/Users/ruicosta/MEGAsync/Bern_projects/anxiety_cells/code';
basepath = 'D:\Bristol\MSC dissertation\Stephanes data\';
%metadata = '/Users/ruicosta/MEGAsync/Bern_projects/anxiety_cells/code/fromCiocchi/data/database dCA1_detailed.csv';
metadata = 'D:\Bristol\MSC dissertation\Stephanes data\database dCA1_detailed.xls';
%[ndata, text, alldata] = xlsread(metadata);
alldata = readtable(metadata);
alldata = table2cell(alldata);

%Indexes (new) database
idx_dbn_rat = 1;
idx_dbn_day = 2;
idx_dbn_session = 3;
idx_dbn_srate = 4;
idx_dbn_tetrode = 5;
idx_dbn_neuron = 6;
idx_dbn_neurontype = 7;
idx_dbn_trial = 8;
idx_dbn_correct = 9;
idx_dbn_strat = 10;
idx_dbn_correctarm = 11;
idx_dbn_rewardzone = 12;
idx_dbn_startzone = 13;
idx_dbn_starttime = 14;
idx_dbn_endtime = 15;
idx_dbn_unk = 16;
idx_dbn_traj = 17;
idx_dbn_trialsomit = 18;
idx_dbn_hist_file = 19;
idx_dbn_beh_file = 20;
idx_dbn_track_file = 21;


plotOn = 0; %Turn on plotting of the firing rates
rule_switch = 0; %0: both 1: allocentric, 2: egocentric
inh_exc = 0; %0: both   1: exc   2: inh
ei_th = 7;

%% Select one of this IDs (each one is from a recording day)
% animals = {'SC03' 'SC08', 'SC11', 'JP28' }; %Animal id
% animals = {'JP28'}; %Animal id
% id = {[animals{1} '-0508-0117_resampled'],[animals{1} '-0511-0118_resampled'],[animals{1} '-0522-0119_resampled'],[animals{1} '-0524-0119_resampled']};
%  id = {[animals{1} '-0508-0117_resampled'], [animals{1} '-0511-0118_resampled'], [animals{1} '-0515-0116_resampled'], [animals{1} '-0517-0120_resampled'], [animals{1} '-0522-0119_resampled'], [animals{1} '-0524-0119_resampled']}; %SC08

%% Individual neurons
% id = {[animals{1} '-0524-0119_resampled']}; 

%% All other values
% id = {[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled']};  %SC03
%id = {[animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled']};  %SC03
% id = {[animals{1} '-0924-0119_resampled'], [animals{1} '-0926-0115_resampled'], [animals{1} '-0927-0115_resampled'], [animals{1} '-1002-0112_resampled']}; %SC11
%id = {[animals{1} '-0924-0119_resampled']}; %SC11
% id = {[animals{1} '-2011-0119_resampled'], [animals{1} '-2211-0115_resampled'], [animals{1} '-1012-0113_resampled']}; %JP28

%% alldata
% id={[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'],...
%     [animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled'], [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'], [animals{2} '-0522-0119_resampled'], [animals{2} '-0524-0119_resampled'],...
%     [animals{3} '-0924-0119_resampled'], [animals{3} '-0926-0115_resampled'], [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled'],...
%     [animals{4} '-2011-0119_resampled'], [animals{4} '-2211-0115_resampled'], [animals{4} '-1012-0113_resampled']};
% id = {[animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled'], [animals{2} '-0515-0116_resampled'], [animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'], [animals{3} '-0924-0119_resampled'], [animals{3} '-0926-0115_resampled'], [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled', [animals{4} '-2011-0119_resampled'], [animals{4} '-2211-0115_resampled'], [animals{4} '-1012-0113_resampled']]};  %SC08+SC03


%%First 3 days
% id={[animals{1} '-0508-0117_resampled'], [animals{1} '-0511-0118_resampled'], [animals{1} '-0515-0116_resampled'],...
%     [animals{2} '-0419-0120_resampled'], [animals{2} '-0422-0120_resampled'], [animals{2} '-0427-0117_resampled'],...
%     [animals{3} '-0924-0119_resampled'],[animals{3} '-0926-0115_resampled'],...
%     [animals{4} '-2011-0119_resampled']};

%% only first days
% id={[animals{1} '-0508-0117_resampled'], ...
%     [animals{2} '-0419-0120_resampled'], ...
%     [animals{3} '-0924-0119_resampled'],...
%     [animals{4} '-2011-0119_resampled']};
%% middle days
% id={[animals{1} '-0507-0117_resampled'], ...
%     [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'],...
%     [animals{3} '-0926-0115_resampled'],[animals{3} '-0927-0115_resampled'], ...
%     [animals{4} '-2211-0115_resampled']};
%% Only last day
% id={ [animals{1} '-0524-0119_resampled'],...
%      [animals{2} '-0531-0114_resampled'],...
%      [animals{3} '-1002-0112_resampled'],...
%      [animals{4} '-1012-0113_resampled']};
%%
% id = {[animals{1} '-0508-0117_resampled'], [animals{1} '-0511-0118_resampled'], [animals{1} '-0515-0116_resampled']}; %SC08 1st 3 days
% id = {[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled']}; %SC03 1st 3 days
% id = {[animals{1} '-0924-0119_resampled'],[animals{1} '-0926-0115_resampled']};
% id={[animals{1} '-2011-0119_resampled']}; %JP28
%%Last 3 days
% id={[animals{1} '-0517-0120_resampled'], [animals{1} '-0522-0119_resampled'], [animals{1} '-0524-0119_resampled'],...
%     [animals{2} '-0516-0112_resampled'], [animals{2} '-0523-0112_resampled'], [animals{2} '-0531-0114_resampled'],...
%     [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled'],...
%     [animals{4} '-1012-0113_resampled']};



% id = {[animals{1} '-0517-0120_resampled'], [animals{1} '-0522-0119_resampled'], [animals{1} '-0524-0119_resampled']}; %SC08 last 3 days
% id = {[animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled']}; % SC03
% id={[animals{1} '-0927-0115_resampled'], [animals{1} '-1002-0112_resampled']}; SC11
% id={[animals{1} '-1012-0113_resampled']}; %JP28


% id={[animals{1} '-2211-0115_resampled']};
%% 33% of data files from each animal for dca1
% animals = {'SC03' 'SC08', 'SC11', 'JP28' }; %Animal id
% id={[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'] ...
%     [animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled']...
%     [animals{3} '-0924-0119_resampled'],...
%     [animals{4} '-2011-0119_resampled']};
%% 33% of mid days for all animals dca1 data
% animals = {'SC03' 'SC08', 'SC11','JP28' }; %Animal id 
% id={[animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], ...
%     [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'], ...
%     [animals{3} '-0926-0115_resampled'], ...
%     [animals{4} '-2211-0115_resampled']}; 
%% 33% of last days for all animals dca1 data
animals = {'SC03' 'SC08', 'SC11', 'JP28' }; %Animal id 
id={[animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'], ...
    [animals{2} '-0522-0119_resampled'], [animals{2} '-0524-0119_resampled'], ...
    [animals{3} '-1002-0112_resampled'], ...
    [animals{4} '-1012-0113_resampled']}; 
%% Main code below
corr_res_all = [];
corr_mres_all = [];
corr_vres_all = [];

incorr_res_all = [];
incorr_mres_all = [];
incorr_vres_all = [];
mres_all=[];
across_trials_var=[];
across_trials=[];

unitcounter=0;
daycounter=0;
%basepath = '/Users/ruicosta/Documents/Bern/ciocchi_Data';
for i=1:size(id,2)
    animal = id{i};
    animal = animal(1:4);
    disp(id{i});
    %Filter by animal and day, and by correct vs incorrect trials
    filtered = find(strcmpi(alldata(:,idx_dbn_rat),animal) & strcmpi(alldata(:,idx_dbn_day), id{i}));
    daycounter=daycounter+1;
    fset = alldata(filtered, :);
%     fprintf("Count %d", i);
%     disp(fset);
    %Get tetrodes
    tet = unique(fset(:, idx_dbn_tetrode));

    %Build data from correct trials
    for j=1:size(tet,1) %Concatenate all the data into variables res_all, etc.
        t = tet{j};

%         if(plotOn)
%             figure('Name', ['tetrode ' num2str(t)])
%         end

        %Get neurons
        tidx = find(strcmpi(fset(:, idx_dbn_tetrode), t)); %gives the number of rows belonging to that tetrode value
        tset = fset(tidx,:); % gets all the rows and columns for only that tetrode

        z = 1;
        while z<=size(tset,1) %For this tetrode
            ni = tset{z,idx_dbn_neuron}; %get one neuron number at a time for that tetrode
            unitcounter=unitcounter+1;
            ni_filt = find(cell2mat(tset(:,idx_dbn_neuron))==ni); %get number of rows for that neuron number
            nset = tset(ni_filt,:); %get all rows and columns for that neuron number
%             disp(nset{1,idx_dbn_trialsomit})
            if(~strcmpi(nset{1,idx_dbn_trialsomit}, 'NaN'))
                nset(str2num(nset{1,idx_dbn_trialsomit}),:) = []; %remove trials to omit
            end
%             if(~isnan(nset{1,idx_dbn_trialsomit}))
%                 nset(str2num(nset{1,idx_dbn_trialsomit}),:) = []; %remove trials to omit
%             end
%             if(~isnan(nset{1,idx_dbn_trialsomit}))
%                 nset((nset{1,idx_dbn_trialsomit}),:) = []; %remove trials to omit
%             end


            % extract all trials for one neuron
            filtered_trials = find((cell2mat(nset(:, idx_dbn_correct)) == 1) | (cell2mat(nset(:, idx_dbn_correct)) == 0));
%             filtered_trials = find((cell2mat(nset(:, idx_dbn_correct)) == 0) );


            %separate the correct and incorrect outcomes from the rows
            %of the neuron subset
            if(rule_switch==0)
                 
                filtered_corr = find((cell2mat(nset(:, idx_dbn_correct)) == 1));
                filtered_incorr = find((cell2mat(nset(:, idx_dbn_correct)) == 0));
            elseif(rule_switch==1) %Allocentric
                filtered_corr = find((cell2mat(nset(:, idx_dbn_correct)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
                filtered_incorr = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
            else %Egocentric
                filtered_corr = find((cell2mat(nset(:, idx_dbn_correct)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) > 2));
                filtered_incorr = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_strat)) > 2));
            end 
%             keyboard
            %pull data for 2 decisions and 2 strategy
            %allocentric
                filt_corr_strat1 = find((cell2mat(nset(:, idx_dbn_correct)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
%                 filt_corr_strat2 = find((cell2mat(nset(:, idx_dbn_correct)) == (dec==1)) & (cell2mat(nset(:, idx_dbn_strat))==2));
            %egocentric
                filt_corr_strat3 = find((cell2mat(nset(:, idx_dbn_correct)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) > 2));
%                 filt_corr_strat4 = find((cell2mat(nset(:, idx_dbn_correct)) == (dec==1)) & (cell2mat(nset(:, idx_dbn_strat))==4));
            %allocentric
                filt_incorr_strat1 = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
               % filt_incorr_strat2 = find((cell2mat(nset(:, idx_dbn_correct)) == (dec==0)) & (cell2mat(nset(:, idx_dbn_strat))==2));
            %egocentric              
               filt_incorr_strat3 = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_strat)) > 2));
               % filt_incorr_strat4 = find((cell2mat(nset(:, idx_dbn_correct)) == (dec==0)) & (cell2mat(nset(:, idx_dbn_strat))==4));
         
            
    
            %try 
            %open the firing rate data file for that neuron, this file
            %contains data as trials vs timepoints, so row 1 is trial one
            %and there are 31 timepoints including the reward point at
            %timepoint 11
                file = [basepath '/data/' tset{1, idx_dbn_rat} '/' tset{1,idx_dbn_day} '/' nset{1, idx_dbn_hist_file}];
                load(file, '-mat');
            %pull out the corect trials from this firing rate data file
                corr_hist = histSave(filtered_corr,:);
            %pull out the incorrect trials from this firing rate data file
                incorr_hist = histSave(filtered_incorr,:);
                %generates 30 time points with a gap of 0.2s -2seconds
                %before and 4seconds after the reward
                times = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
                % adds values across columns and divides by the number of trials,
                %trial averages for each time point for each neuron -->
                %neuron-wise trial average
                corr_mres = mean(corr_hist,1); 
                incorr_mres = mean(incorr_hist,1);
                corr_vres = std(corr_hist,1)./sqrt(size(corr_hist,1));
                incorr_vres = std(incorr_hist,1)./sqrt(size(incorr_hist,1));
                %vres = std(Values_Reward_TIME{i,3},1);
                
                neuron_trials=histSave(filtered_trials,:);%filtered_trials
                neuron_trials(:,11)=[];
                
                across_trials=[across_trials; neuron_trials];
                
                acr_trial_var=std(neuron_trials,1)./sqrt(size(neuron_trials,1));
                across_trials_var=[across_trials_var;acr_trial_var];
                
                mres=mean(neuron_trials,1);
                mres_all=[mres_all;mres];
%                 mres_all( mres_all < 0 ) = 0;

%                 disp(zscore(corr_mres')');
%                 keyboard
                %delete the reward threshold column
                corr_mres(:,11)=[];
                incorr_mres(:,11)=[];
                %pull out the subset of trials with different decisions and
                %strategies
                corr_strat1=histSave(filt_corr_strat1,:); %dec 1
                %corr_strat2=histSave(filt_corr_strat2,:);
                corr_strat3=histSave(filt_corr_strat3,:);
%                 corr_strat4=histSave(filt_corr_strat4,:);              
                
                incorr_strat1=histSave(filt_incorr_strat1,:); %dec 0
%                 incorr_strat2=histSave(filt_corr_strat2,:);
                incorr_strat3=histSave(filt_incorr_strat3,:);
%                 incorr_strat4=histSave(filt_corr_strat4,:); 
                                
                %calculate trial averages for this neuron for each correct
                %and strategy combination
                corr_strat1(:,11)=[];
%                 corr_strat2(:,11)=[];
                corr_strat3(:,11)=[];
%                 corr_strat4(:,11)=[];
                
                corr_strat1_mres = mean(corr_strat1,1); 
%                 corr_strat2_mres = mean(corr_strat2,1); 
                corr_strat3_mres = mean(corr_strat3,1); 
%                 corr_strat4_mres = mean(corr_strat4,1); 
                
                %calculate trial averages for this neuron for each incorrect
                %and strategy combination
                incorr_strat1(:,11)=[];
%                 incorr_strat2(:,11)=[];
                incorr_strat3(:,11)=[];
%                 incorr_strat4(:,11)=[];
               
                incorr_strat1_mres = mean(incorr_strat1,1); 
%                 incorr_strat2_mres = mean(incorr_strat2,1); 
                incorr_strat3_mres = mean(incorr_strat3,1); 
%                 incorr_strat4_mres = mean(incorr_strat4,1); 

                % % Neurons*day*decision*strategy*time*trialaverages
                X_full_daywise(unitcounter,1,1,1:size(incorr_mres,2))=incorr_mres; %dec 0
%                 X_full_daywise(unitcounter,daycounter,1,2,1:size(incorr_strat3_mres,2))=incorr_strat3_mres; %dec 0
                X_full_daywise(unitcounter,1,2,1:size(corr_mres,2))=corr_mres; %dec 1
%                 X_full_daywise(unitcounter,daycounter,2,2,1:size(corr_strat3_mres,2))=corr_strat3_mres; %dec 1

               %%Neurons*dec*strategy*time*trialavg
                X_full_avg(unitcounter,1,1,1:size(incorr_strat1_mres,2))=incorr_strat1_mres; %dec 0
%                 X_full_avg(unitcounter,1,2,1:size(incorr_strat2_mres,2))=incorr_strat2_mres;
                trialNum(unitcounter, 1, 1) = size(incorr_strat1, 1);

                X_full_avg(unitcounter,1,2,1:size(incorr_strat3_mres,2))=incorr_strat3_mres;
%                 X_full_avg(unitcounter,1,4,1:size(incorr_strat4_mres,2))=incorr_strat4_mres;
                trialNum(unitcounter, 1, 2) = size(incorr_strat3, 1);

                X_full_avg(unitcounter,2,1,1:size(corr_strat1_mres,2))=corr_strat1_mres; %dec 1
%                 X_full_avg(unitcounter,2,2,1:size(corr_strat2_mres,2))=corr_strat2_mres;
                trialNum(unitcounter, 2, 1) = size(corr_strat1, 1);

                X_full_avg(unitcounter,2,2,1:size(corr_strat3_mres,2))=corr_strat3_mres;
%                 X_full_avg(unitcounter,2,4,1:size(corr_strat4_mres,2))=corr_strat4_mres;
                trialNum(unitcounter, 2, 2) = size(corr_strat3, 1);

                
                XFiringRates(unitcounter,1,1,1:size(incorr_strat1,2),1:size(incorr_strat1,1))=incorr_strat1';
                XFiringRates(unitcounter,1,2,1:size(incorr_strat3,2),1:size(incorr_strat3,1))=incorr_strat3';
                XFiringRates(unitcounter,2,1,1:size(corr_strat1,2),1:size(corr_strat1,1))=corr_strat1';
                XFiringRates(unitcounter,2,2,1:size(corr_strat3,2),1:size(corr_strat3,1))=corr_strat3';
                
                X_full_avg(isnan(X_full_avg)) = 0;
                X_full_daywise(isnan(X_full_daywise)) = 0;
                corr_res_all = [corr_res_all; corr_hist]; %simply adds one after another all the values for each correct trial for each neuron
                corr_mres_all = [corr_mres_all; corr_mres];
                corr_vres_all = [corr_vres_all; corr_vres];

                incorr_res_all = [incorr_res_all; incorr_hist];
                incorr_mres_all = [incorr_mres_all; incorr_mres];
                incorr_vres_all = [incorr_vres_all; incorr_vres];
                
            z=z+size(ni_filt, 1);
            %catch me
            %    continue;
            %end    
%             break
        end 
%         break
    end    
%     break
end
keyboard
corr_res_all(:,11) = []; %Remove NaN (reward)
incorr_res_all(:,11) = []; %Remove NaN (reward)
% corr_mres_all(:,11) = []; %Remove NaN (reward)
incorr_mres_all(:,11) = []; %Remove NaN (reward)
corr_vres_all(:,11) = []; %Remove NaN (reward)
incorr_vres_all(:,11) = []; %Remove NaN (reward)
times(:,11) = [];

%% export data matrix for pca
zmres_all=zscore(mres_all')';
save('dCA1_alldata_trialaveraged.mat','zmres_all');
%% non-negative

[Scores,PCS] = nnmf(mres_all,3);
figure;
htn=[];
for i =1:3
hn=plot(times,PCS(i,:),'Linewidth', 2);
htn=[htn hn];
hold on;
end
title(['NMF across trial averaged dCA1 data']);

plot([0 0], [-1 1], ':k', 'Linewidth', 1);
box off;
ylim([-0.4 0.5]);
ylabel('norm. activity [-]');
xlabel('Time(s)');
legend(htn, {['1st comp.'],['2nd comp.'],['3rd comp.']});

figure;
PC_T=PCS(1:2,:);
Scores_T=Scores(:,1:2);
biplot(PC_T','scores',Scores_T);
% axis([0 1.1 0 1.1])
xlabel('Column 1')
ylabel('Column 2')
if(1)
%% PCA
[corr_comps_mall, corr_score, latent,tsquared,corr_explained_mall] = pca(zscore(mres_all')'); %Across trial average
keyboard
[acr_trial_var_comps_all, acr_trial_var_score, acr_trial_var_latent,acr_trial_var_tsquared,acr_trial_var_explained_vall] = pca(zscore(across_trials_var')'); %Across trial variance
[acr_trial_comps_all, acr_trial__score, acr_trial_latent,acr_trial_tsquared,acr_trial_explained_all] = pca(zscore(across_trials')'); %Across trials
disp(mean(mres_all'));
keyboard
pcomp = 3;
cs = copper(pcomp);
set(groot,'defaultAxesColorOrder',cs)

figure('Position', [1317        -331         953         901], 'Name', [cell2mat(id) ' rule: ' num2str(rule_switch)]);

ht = [];
% subplot(3,3,1)
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

figure('Position', [1317        -331         953         901], 'Name', [cell2mat(id) ' rule: ' num2str(rule_switch)]);

ht1 = [];
% subplot(1,3,2)
for i=1:pcomp
    h = plot(times, acr_trial_var_comps_all(:,i), 'Linewidth', 2);
    ht1 = [ht1 h];
    hold on;
end
title(['PCA across trial-variability data (' num2str(sum(acr_trial_var_explained_vall(1:pcomp)),3) '% exp. var.)']);

legend(ht1, {['1st comp. (' num2str(acr_trial_var_explained_vall(1),2) '%)'], ['2nd comp. (' num2str(acr_trial_var_explained_vall(2),2) '%)'], ['3rd comp. (' num2str(acr_trial_var_explained_vall(3),2) '%)']})
plot([0 0], [-0.4 0.4], ':k', 'Linewidth', 1);
ylim([-0.4 0.5]);
box off;
ylabel('Norm. firing rate [-]');
xlabel('Time ( in seconds )');

figure('Position', [1317        -331         953         901], 'Name', [cell2mat(id) ' rule: ' num2str(rule_switch)]);
ht1 = [];
% subplot(1,3,3)
for i=1:pcomp
    h = plot(times, acr_trial_comps_all(:,i), 'Linewidth', 2);
    ht1 = [ht1 h];
    hold on;
end
title(['PCA across trials data (' num2str(sum(acr_trial_explained_all(1:pcomp)),3) '% exp. var.)']);

legend(ht1, {['1st comp. (' num2str(acr_trial_explained_all(1),2) '%)'], ['2nd comp. (' num2str(acr_trial_explained_all(2),2) '%)'], ['3rd comp. (' num2str(acr_trial_explained_all(3),2) '%)']})
plot([0 0], [-0.5 0.5], ':k', 'Linewidth', 1);
ylim([-0.5 0.5]);
box off;
ylabel('Norm. firing rate [-]');
xlabel('Time ( in seconds )');

end

%% saving ventral data
    save('dCA1_Xfull_Avg','X_full_avg');
% save('dCA1_Xfull_nostrat_FirstDay','X_full_daywise');
if(1)
%% DPCA on decision, strategy, time and interaction of decision and strategy
Xfull=(X_full_avg); %original matlab matrix
% load('D:\Bristol\MSC dissertation\TD learning\DRQN_Activations\NO_ROT_ACTIVATIONS\Xfull_drqn_2states_50neurons.mat');
% Xfull=Xfull; %python matrix

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
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
times1 = -2:0.2:((31*0.2)-2.2);
% times=-2:0.2:((30*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);
% timeEvents = times1(30);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda', 2.919292602539062e-05); %2.919292602539062e-05
% keyboard
explVar = dpca_explainedVariance(Xfull, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(Xfull, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', times1,...
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
    'numRep', 20, ...  % increase this number to ~10 for better accuracy
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
optimalLamdas = [0.000019462 0.000043789 0.000029193 1.477891880035400e-04 0.0025];
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