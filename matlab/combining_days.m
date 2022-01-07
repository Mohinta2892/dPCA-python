%%
% Quick analysis of firing rate data of neurons from behaving animals
%

close all
clear all;

% loadColors;

tetrodes = [1:16]; %Tetrodes id used to record neural activity

%basepath = '/Users/ruicosta/MEGAsync/Bern_projects/anxiety_cells/code';
basepath = 'D:\Bristol\MSC dissertation\Stephanes data\ventral data';
%metadata = '/Users/ruicosta/MEGAsync/Bern_projects/anxiety_cells/code/fromCiocchi/data/database dCA1_detailed.csv';
metadata = 'D:\Bristol\MSC dissertation\Stephanes data\database vCA1_detailed.xls';
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
% animals = {'SC03' 'SC08', 'SC11','SC12','JP28' }; %Animal id
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
% id={[animals{1} '-1105-0116_resampled'],[animals{1} '-1107-0115_resampled'],[animals{1} '-1113-0114_resampled']};
%% alldata
% id={[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'],...
%     [animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled'], [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'], [animals{2} '-0522-0119_resampled'], [animals{2} '-0524-0119_resampled'],...
%     [animals{3} '-0924-0119_resampled'], [animals{3} '-0926-0115_resampled'], [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled'],...
%     [animals{4} '-1105-0116_resampled'],[animals{4} '-1107-0115_resampled'],[animals{4} '-1113-0114_resampled'],...
%     [animals{5} '-2011-0119_resampled'], [animals{5} '-2211-0115_resampled'], [animals{5} '-1012-0113_resampled']};
% % id = {[animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled'], [animals{2} '-0515-0116_resampled'], [animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], [animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'], [animals{3} '-0924-0119_resampled'], [animals{3} '-0926-0115_resampled'], [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled', [animals{4} '-2011-0119_resampled'], [animals{4} '-2211-0115_resampled'], [animals{4} '-1012-0113_resampled']]};  %SC08+SC03
% id={[animals{1} '-1105-0116_resampled'],[animals{1} '-1107-0115_resampled'],[animals{1} '-1113-0114_resampled']}; %SC12

%%First 3 days
% id={[animals{1} '-0508-0117_resampled'], [animals{1} '-0511-0118_resampled'], [animals{1} '-0515-0116_resampled'],...
%     [animals{2} '-0419-0120_resampled'], [animals{2} '-0422-0120_resampled'], [animals{2} '-0427-0117_resampled'],...
%     [animals{3} '-0924-0119_resampled'],[animals{3} '-0926-0115_resampled'],...
%     [animals{4} '-1105-0116_resampled'], ...
%     [animals{5} '-2011-0119_resampled']};

%% only first days
% id={[animals{1} '-0508-0117_resampled'], ...
%     [animals{2} '-0419-0120_resampled'], ...
%     [animals{3} '-0924-0119_resampled'],...
%     [animals{4} '-1105-0116_resampled'],...
%     [animals{5} '-2011-0119_resampled']};
%% middle days
% id={[animals{1} '-0507-0117_resampled'], ...
%     [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'],...
%     [animals{3} '-0926-0115_resampled'],[animals{3} '-0927-0115_resampled'], ...
%       [animals{4} '-1107-0115_resampled'],...
%     [animals{5} '-2211-0115_resampled']};
%% Only last day
% id={ [animals{1} '-0524-0119_resampled'],...
%      [animals{2} '-0531-0114_resampled'],...
%      [animals{3} '-1002-0112_resampled'],...
%      [animals{4} '-1113-0114_resampled'],...
%      [animals{5} '-1012-0113_resampled']};
%%
% id = {[animals{1} '-0508-0117_resampled'], [animals{1} '-0511-0118_resampled'], [animals{1} '-0515-0116_resampled']}; %SC08 1st 3 days
% id = {[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], [animals{1} '-0427-0117_resampled']}; %SC03 1st 3 days
% id = {[animals{1} '-0924-0119_resampled'],[animals{1} '-0926-0115_resampled']};
% id={[animals{1} '-2011-0119_resampled']}; %JP28
%%Last 3 days
% id={[animals{1} '-0517-0120_resampled'], [animals{1} '-0522-0119_resampled'], [animals{1} '-0524-0119_resampled'],...
%     [animals{2} '-0516-0112_resampled'], [animals{2} '-0523-0112_resampled'], [animals{2} '-0531-0114_resampled'],...
%     [animals{3} '-0927-0115_resampled'], [animals{3} '-1002-0112_resampled'],...
%     [animals{4} '-1113-0114_resampled'],...
%     [animals{5} '-1012-0113_resampled']};



% id = {[animals{1} '-0517-0120_resampled'], [animals{1} '-0522-0119_resampled'], [animals{1} '-0524-0119_resampled']}; %SC08 last 3 days
% id = {[animals{1} '-0516-0112_resampled'], [animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled']}; % SC03
% id={[animals{1} '-0927-0115_resampled'], [animals{1} '-1002-0112_resampled']}; SC11
% id={[animals{1} '-1012-0113_resampled']}; %JP28


% id={[animals{1} '-2211-0115_resampled']};
%% 33% of first days for all animals vca1 data
% animals = {'SC03' 'SC08', 'SC11', 'SC12', 'JP28' }; %Animal id 
% id={[animals{1} '-0419-0120_resampled'], [animals{1} '-0422-0120_resampled'], ...
%     [animals{2} '-0508-0117_resampled'], [animals{2} '-0511-0118_resampled'], ...
%     [animals{3} '-0924-0119_resampled'], ...
%     [animals{4} '-1105-0116_resampled'], ...
%     [animals{5} '-2011-0119_resampled']}; 
% %% 33% of mid days for all animals vca1 data
% animals = {'SC03' 'SC08', 'SC11', 'SC12', 'JP28' }; %Animal id 
% id={[animals{1} '-0427-0117_resampled'], [animals{1} '-0507-0117_resampled'], ...
%     [animals{2} '-0515-0116_resampled'], [animals{2} '-0517-0120_resampled'], ...
%     [animals{3} '-0926-0115_resampled'], ...
%     [animals{4} '-1107-0115_resampled'], ...
%     [animals{5} '-2211-0115_resampled']}; 
%% 33% of mid days for all animals vca1 data
% animals = {'SC03' 'SC08', 'SC11', 'SC12', 'JP28' }; %Animal id 
% id={[animals{1} '-0523-0112_resampled'], [animals{1} '-0531-0114_resampled'], ...
%     [animals{2} '-0522-0119_resampled'], [animals{2} '-0524-0119_resampled'], ...
%     [animals{3} '-1002-0112_resampled'], ...
%     [animals{4} '-1113-0114_resampled'], ...
%     [animals{5} '-1012-0113_resampled']}; 
%% Main code starts here
corr_res_all = [];
corr_mres_all = [];
corr_vres_all = [];

incorr_res_all = [];
incorr_mres_all = [];
incorr_vres_all = [];
mres_all=[];
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
%         keyboard;
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
            
            %start zone startegy  based data layout
            filt_strat12_start1 = find((cell2mat(nset(:, idx_dbn_startzone)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
            filt_strat12_start2 = find((cell2mat(nset(:, idx_dbn_startzone)) == 2) & (cell2mat(nset(:, idx_dbn_strat)) <= 2));
            filt_strat34_start1 = find((cell2mat(nset(:, idx_dbn_startzone)) == 1) & (cell2mat(nset(:, idx_dbn_strat)) >  2));
            filt_strat34_start2 = find((cell2mat(nset(:, idx_dbn_startzone)) == 2) & (cell2mat(nset(:, idx_dbn_strat)) > 2));

            %decision start zone based data layout
            filt_corr_start1 = find((cell2mat(nset(:, idx_dbn_correct)) == 1) & (cell2mat(nset(:, idx_dbn_startzone)) == 1));
            filt_incorr_start1 = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_startzone)) == 1));
            filt_corr_start2 = find((cell2mat(nset(:, idx_dbn_correct)) == 1)  & (cell2mat(nset(:, idx_dbn_startzone)) == 2));
            filt_incorr_start2 = find((cell2mat(nset(:, idx_dbn_correct)) == 0) & (cell2mat(nset(:, idx_dbn_startzone)) == 2));

            
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

                %STRATEGY START ZONE 
                strat12_start1=histSave(filt_strat12_start1,:);
                strat12_start2=histSave(filt_strat12_start2,:);
                strat34_start1=histSave(filt_strat34_start1,:);
                strat34_start2=histSave(filt_strat34_start2,:);
                                
                
                %decision start zone
                corr_start1=histSave(filt_corr_start1,:);
                incorr_start1=histSave(filt_incorr_start1,:);
                corr_start2=histSave(filt_corr_start2,:);
                incorr_start2=histSave(filt_incorr_start2,:);
                
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

                %STRATEGY START ZONE
                strat12_start1(:,11)=[];
                strat12_start2(:,11)=[];
                strat34_start1(:,11)=[];
                strat34_start2(:,11)=[];
                
                                %decision start zone
                corr_start1(:,11)=[];
                incorr_start1(:,11)=[];
                corr_start2(:,11)=[];
                incorr_start2(:,11)=[];
               
                incorr_strat1_mres = mean(incorr_strat1,1); 
%                 incorr_strat2_mres = mean(incorr_strat2,1); 
                incorr_strat3_mres = mean(incorr_strat3,1); 
%                 incorr_strat4_mres = mean(incorr_strat4,1); 

                %strategy start zone
                strat12_start1_mres = mean(strat12_start1,1);
                strat12_start2_mres = mean(strat12_start2,1); 
                strat34_start1_mres = mean(strat34_start1,1); 
                strat34_start2_mres = mean(strat34_start2,1); 
                
                
                %decision start zone
                corr_start1_mres = mean(corr_start1,1);
                corr_start2_mres = mean(corr_start2,1);
                incorr_start1_mres = mean(incorr_start1,1);
                incorr_start2_mres = mean(incorr_start2,1);

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

                %Neurons x strategy x start zone x time x trial average
                X_full_avg_strategy_startzone(unitcounter,1,1,1:size(strat12_start1_mres,2))=strat12_start1_mres; %allo start 1
                X_full_avg_strategy_startzone(unitcounter,1,2,1:size(strat12_start2_mres,2))=strat12_start2_mres; %allo start 2
                X_full_avg_strategy_startzone(unitcounter,2,1,1:size(strat34_start1_mres,2))=strat34_start1_mres; %ego start 2
                X_full_avg_strategy_startzone(unitcounter,2,2,1:size(strat34_start2_mres,2))=strat34_start2_mres; %ego start 2

               %Neurons x dec x start zone x time x trial average
                X_full_avg_dec_startzone(unitcounter,1,1,1:size(corr_start1_mres,2))=corr_start1_mres; %allo start 1
                X_full_avg_dec_startzone(unitcounter,1,2,1:size(corr_start2_mres,2))=corr_start2_mres; %allo start 2
                X_full_avg_dec_startzone(unitcounter,2,1,1:size(incorr_start1_mres,2))=incorr_start1_mres; %ego start 2
                X_full_avg_dec_startzone(unitcounter,2,2,1:size(incorr_start2_mres,2))=incorr_start2_mres; %ego start 2

                
                XFiringRates(unitcounter,1,1,1:size(incorr_strat1,2),1:size(incorr_strat1,1))=incorr_strat1';
                XFiringRates(unitcounter,1,2,1:size(incorr_strat3,2),1:size(incorr_strat3,1))=incorr_strat3';
                XFiringRates(unitcounter,2,1,1:size(corr_strat1,2),1:size(corr_strat1,1))=corr_strat1';
                XFiringRates(unitcounter,2,2,1:size(corr_strat3,2),1:size(corr_strat3,1))=corr_strat3';
                
                X_full_avg(isnan(X_full_avg)) = 0;
                X_full_daywise(isnan(X_full_daywise)) = 0;
                X_full_avg_strategy_startzone(isnan(X_full_avg_strategy_startzone)) = 0;
                X_full_avg_dec_startzone(isnan(X_full_avg_dec_startzone))=0;

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
% % corr_mres_all(:,11) = []; %Remove NaN (reward)
incorr_mres_all(:,11) = []; %Remove NaN (reward)
corr_vres_all(:,11) = []; %Remove NaN (reward)
incorr_vres_all(:,11) = []; %Remove NaN (reward)
times(:,11) = [];

%% export data matrix for pca
mstore=zscore(mres_all')';
save('data.mat','mstore');

%% export data vca1 strategy start zone
save('vCA1_dec_startzone_trialavg.mat','X_full_avg_dec_startzone');

if(1)
%% PCA
[corr_comps_mall, corr_score, latent,tsquared,corr_explained_mall] = pca(zscore(mres_all')'); %Across trial average
disp(mean(mres_all'));
% keyboard
pcomp = 4;
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

legend(ht, {['1st comp. (' num2str(corr_explained_mall(1),2) '%)'], ['2nd comp. (' num2str(corr_explained_mall(2),2) '%)'], ['3rd comp. (' num2str(corr_explained_mall(3),2) '%)'], ['4th comp. (' num2str(corr_explained_mall(4),2) '%)']})
plot([0 0], [-0.4 0.4], ':k', 'Linewidth', 1);
ylim([-0.4 0.5]);
box off;
ylabel('Norm. firing rate [-]');
xlabel('Time ( in seconds )');

end
%% saving ventral data
    save('vCA1_Xfull_nostrat_FirstDay','X_full_daywise');

%% loading matrix for combining
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg_firstdays');
dca1_F=X_full_avg; %%X_full_daywise
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg_middays');
dca1_M=X_full_avg;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg_lastdays');
dca1_L=X_full_avg;

load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg_firstdays');
vca1_F=X_full_avg;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg_middays');
vca1_M=X_full_avg;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg_lastdays');
vca1_L=X_full_avg;

load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg');
dca1=X_full_avg;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg');
vca1=X_full_avg;

% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\JP28_dXfull_avg_lastdays');
% dca1=X_full_avg;
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\SC12_vXfull_avg_lastdays');
% vca1=X_full_avg;


combo_F=[dca1_F;vca1_F];
combo_M=[dca1_M;vca1_M];
combo_L=[dca1_L;vca1_L];
CA1=zeros(278,3,2,30);

CA1(1:39,1,1:size(combo_F,3),1:size(combo_F,4))=combo_F;
CA1(40:232,2,1:size(combo_M,3),1:size(combo_M,4))=combo_M;
CA1(233:278,3,1:size(combo_L,3),1:size(combo_L,4))=combo_L;

All_data=[dca1;vca1];

% load('D:\Bristol\MSC dissertation\TD learning\DRQN_Activations\NO_ROT_ACTIVATIONS\Xfull_drqn_2states');
Xfull=Xfull;


if(0)
%% DPCA on decision, strategy, time and interaction of decision and strategy
% Xfull=CA1;
% Xfull=All_data;
% Xfull=vca1;

% Xfull=combo_F;
Xfull=X_full_avg;

% D = size(trialNum,1);
% minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
% meanFiringRate = mean(reshape(Xfull, D, []), 2);
% n = find(minN >= 2 & meanFiringRate < 50);
% t = 1:length(times);

% combinedParams = {{1, [1,2]},{2, [2,3]},{3, [3,4]}};
combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
times1 = -2:0.2:((31*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

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

%% PCA on full data

load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_alldata_trialaverage');
dca1=zmres_all;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_alldata_trialaverage');
vca1=zmres_all;

all_data=[dca1;vca1];

[corr_comps_mall, corr_score, latent,tsquared,corr_explained_mall] = pca(all_data); %Across trial average

pcomp = 4;
cs = copper(pcomp);
set(groot,'defaultAxesColorOrder',cs)

% figure('Position', [1317        -331         953         901], 'Name', [cell2mat(id) ' rule: ' num2str(rule_switch)]);
figure;
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
times = -2:0.2:((31*0.2)-2.2);
times(:,11)=[];
ht = [];
% subplot(2,2,1)
for i=1:pcomp
    h = plot(times, corr_comps_mall(:,i), 'Linewidth', 2);
    ht = [ht h];
    hold on;
end
title(['PCA on trial-averaged CA1 data (' num2str(sum(corr_explained_mall(1:pcomp)),3) '% exp. var.)']);

h  = plot([0 0], [-0.4 0.4], ':k', 'Linewidth', 1);
ht = [ht h];
ylim([-0.4 0.5]);
box off;
ylabel('morm. activity');
xlabel('time (s)');
legend(ht, {['1st comp. (' num2str(corr_explained_mall(1),2) '%)'], ['2nd comp. (' num2str(corr_explained_mall(2),2) '%)'], ['3rd comp. (' num2str(corr_explained_mall(3),2) '%)'], ['4th comp. (' num2str(corr_explained_mall(4),2) '%)'], ['reward zone']})



%% plot 30 neurons

load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_alldata_trialaverage');
dca1=zmres_all;

figure;

times = -2:0.2:((31*0.2)-2.2);
times(:,11)=[];
ht = [];

% subplot(2,2,1)
for i=30:60
    h = plot(times, dca1(i,:), 'Linewidth', 2);
    ht = [ht h];
    hold on;
end
title(['Trial-averaged firing rates of 30 neurons']);

% legend(ht, {['1st comp. (' num2str(corr_explained_mall(1),2) '%)'], ['2nd comp. (' num2str(corr_explained_mall(2),2) '%)'], ['3rd comp. (' num2str(corr_explained_mall(3),2) '%)'], ['4th comp. (' num2str(corr_explained_mall(4),2) '%)']})
plot([0 0], [-3,5],':k', 'Linewidth', 1);
% ylim([-0.4 0.5]);
box off;
ylabel('Norm. firing rate [-]');
xlabel('Time ( in seconds )');


%% DPCA on full data with decision strategy and start zone

load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dca1_dec_strat_start');
d_X_full_all=X_full_all;
load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vca1_dec_strat_start');
v_X_full_all=X_full_all;

all_data=[d_X_full_all;v_X_full_all];

Xfull=(all_data);

% D = size(trialNum,1);
% minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
% meanFiringRate = mean(reshape(Xfull, D, []), 2);
% n = find(minN >= 2 & meanFiringRate < 50);
% t = 1:length(times);

% combinedParams = {{1, [1,2]},{2, [2,3]},{3, [3,4]}};
% combinedParams ={{1, [1 3], [1 4], [1 5]}, {2, [2 3], [2 4], [2 5]}, {3, [3 4], [3 5]}, {4, [4 5]}, {5}};
combinedParams ={{1, [1 3], [1 4]}, {2, [2 3], [2 4]}, {3}, {4}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Start Zone', 'Time'};
% margNames = {'Decision','Strategy', 'Start zone', 'Reward Zone', 'Time'};
% margNames = {'Decision', 'Time'};
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171, ; 110 60 171 ; 50 50 50]/256;
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

load('tmp_optimalLambdas.mat', 'optimalLambdas');
times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 8,...
    'combinedParams',combinedParams,'lambda', 2.919292602539062e-05); %2.919292602539062e-05

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


%% combining first 33% of days for all animals in vca1 and dca1

load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_firstday_dca1\Xfull_dca1_firstdays');
dca1_Xfull=Xfull;
load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_firstday_vca1\Xfull_vca1_firstdays');
vca1_Xfull=Xfull;

Xfull=[dca1_Xfull;vca1_Xfull];

combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
times1 = -2:0.2:((31*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

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

%% combining mid 33% of days for all animals in vca1 and dca1

load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_midday_dca1\Xfull_dca1_middays');
dca1_Xfull=Xfull;
load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_midday_vca1\Xfull_vca1_middays');
vca1_Xfull=X_full_avg;

Xfull=[dca1_Xfull;vca1_Xfull];

combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
times1 = -2:0.2:((31*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

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

%% combining mid 33% of days for all animals in vca1 and dca1

load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_lastday_dca1\Xfull_dca1_lastdays');
dca1_Xfull=Xfull;
load('D:\Bristol\Code_NatureComm\Mat\All_data_animals\33_percent_lastday_vca1\Xfull_vca1_lastdays');
vca1_Xfull=Xfull;

Xfull=[dca1_Xfull;vca1_Xfull];

combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
times1 = -2:0.2:((31*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

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

%% dpca on model activations

load('/home/samia/Documents/dgx2/bristol_continual_learning/experiment_data/Partial_ENV/partial_env3_dabal/new_5x5_with_act/activations_hcdqn_neuron_data/dpca_mat_files_aroundreward_smooth0/Xfull_CA1_allseeds.mat');
dca1_Xfull=Xfull;
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\dCA1_Xfull_Avg');
% dca1_Xfull=X_full_avg;
% load('D:\Bristol\MSC dissertation\elife2016dpca-master\StaticComponents\vCA1_Xfull_Avg');
% vca1_Xfull=X_full_avg;

% Xfull=[dca1_Xfull;vca1_Xfull];


combinedParams ={{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% combinedParams ={{1, [1 3]}}; %gives blank why?
margNames = {'Decision','Strategy', 'Time', 'Interaction'};
% margNames = {'Decision', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
load('tmp_optimalLambdas.mat', 'optimalLambdas');
% times1 = -2:0.2:((size(corr_hist,2).*0.2)-2.2);
% times1 = -2:0.2:((31*0.2)-2.2);
times1 = -2:0.2:((201*0.2)-2.2);
times1(:,11) = [];

timeEvents = times1(11);

[W,V,whichMarg] = dpca(Xfull, 16,...
    'combinedParams',combinedParams,'lambda',2.919292602539062e-05); %2.919292602539062e-05 6.568408355712890e-05

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
