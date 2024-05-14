% This script performs the plotting and statistical tests from the
% manuscript "A paradigm for characterizing motion misperception in people
% with typical vision and low vision"

%% set up
close all; clear all;

% which group to analyze
grp = 'nv'; % typical vision
%grp = 'lv'; % low vision

% plot each individual participant's data?
plotIndividuals = 1;

% add path to required toolboxes
addpath(genpath('helper_functions'))

% Note: to run this code, you will need to download the CircStats toolbox
% and add it to your path:
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

% add path to data files
addpath(genpath('data'));

% data structure for each subject is as follows:
%
% dat.stim columns = contrast, speed, orientation, true direction
% dat.resp = response: direction axis 1-12 or 0 (no motion)
%
% for orientation,
% +45 = tilted top right (clockwise), -45 = tilted top left (counter
% clockwise)
%
% for true direction,
% 90 = up/down, 0 = right/left, 45 = down left to up right, -45 = down
% right to up left

% angular value associated with each possible response of 1-12, used to
% convert responses into angular motion directions
angles = [90 75 60 45 30 15 0 -15 -30 -45 -60 -75];

% mirrored angles for visualization of back-and-forth motion
angles_mirror = mod(angles+180,360);

% edges of angular histogram bins
angles_edges = -7:15:360-7;

% look up tables for converting responses when combining across different conditions
rspConverterMirrorOrientation = [1 12 11 10 9 8 7 6 5 4 3 2]; % mirror responses over the horizontal axis

% colors for low contrast (red) and high contrast (cyan) in polar histograms
ccolor = [1 0 0; 0 1 1];

% get file names for the group we're going to analyze
switch grp
    case 'nv'
        % load subjects with normal vision,
        filenames = {'LV_NBLV1.mat' 'LV_NBLV2.mat' 'LV_NBLV3.mat' 'LV_NBLV4.mat' 'LV_NBLV5.mat'};
        subjTags = {'NV1' 'NV2' 'NV3' 'NV4' 'NV5'};
    case 'lv'
        % load subjects with low vision,
        filenames = { 'LV_BLV3.mat' 'LV_BLV1.mat' 'LV_BLV2.mat' 'LV_BLV5.mat' 'LV_BLV4.mat'};
        subjTags = {'LV1' 'LV2' 'LV3' 'LV4' 'LV5'};
        subjTags = {'OA' 'OAG' 'OC1' 'OC2' 'SD'};
    otherwise
        error('no group')
end

% initiate structure where we'll store all data for statistics
testConditionTrials = [];

%% calculate angular errors (and make individual polar plots if requested)

% figure for polar plots
fig_polars = figure;

for s = 1:length(filenames)

    % load in the data
    disp(['subj ' num2str(s)]);
    filename = filenames{s};
    load(filename);
    data = [dat.stim, dat.resp];

    if plotIndividuals
        % open a figure for this subject
        figure; hold on;
        sgtitle(['Subj ' subjTags{s} '(blue = high contrast, red = low contrast)'])
        set(gcf,'Position',[1 1 1400 1000]);
    end

    % fix error in NV5 in which 110 was pressed instead of 11
    if strcmp(filename,'LV_NBLV5.mat')
        data(data(:,5) == 110,5) = 11;
    end

    %each subject has their own contrast threshold, calculate and report it
    lc_moving = round(mean(dat.dynamicThres),3);
    disp(['contrast threshold: ' num2str(lc_moving)])

    % report # of no motion trials
    disp(['no motion responses: ' num2str(sum(dat.resp==0)) '/192 -- ' num2str(100*sum(dat.resp==0)/192,3) '%' ])
    fprintf('\n')

    % store the two contrast levels used for this subject
    contrast = [lc_moving, 1];

    % mirror reflect the responses for the -45 orientation

    % get all negative rhombus orientations and code them as positive
    indNegOrientation = data(:,3) < 0;
    data(indNegOrientation,3) =  -data(indNegOrientation,3);

    % recode the motion direction so it is relative to the shape orientation
    data(indNegOrientation,4) =  -data(indNegOrientation,4);

    % -90 and +90 are the same, so convert any -90's back to +90
    data(abs(data(:,4)+90)<0.001,4) = 90;

    % for trials when the rhombus had a positive rotation, mirror over the
    % horizontal axis to align them with trials with a negative orientation
    data(indNegOrientation & data(:,5)>0,5) = rspConverterMirrorOrientation(data(indNegOrientation & data(:,5)>0,5));

    % collapse across slow and fast trials
    data(:,2) = 1;

    % counters
    indCounter = [];
    subdataAll = [];
    splotCntr  = 1;

    % calculate errors
    for speed = unique(data(:,2))' % for each speed
        for orientation = unique(data(:,3))' % for each rhombus orientation
            for motion_dir = [0 90 45 -45] % for each motion direction

                % set radius maximum for plotting
                max_rad = 0.5;

                % bayesian prediction is orthogonal to rhombus orientation
                bayesianPrediction = -orientation;

                % reference lines and titles for individual plots
                if plotIndividuals

                    % open plot
                    subplot(4,4,splotCntr);

                    %plot the rhombus orientation
                    polarplot(deg2rad([orientation orientation+90 orientation+180 orientation+270 orientation]), [0.4 0.07 0.4 0.07 0.4], 'color', [0.5 0.5 0.5], 'LineWidth', 1);  % Adjust the color and line properties as needed
                    hold on;

                    %plot the ground truth motion
                    polarplot(deg2rad([motion_dir,0,mod(motion_dir+180,360)]),[max_rad 0 max_rad],'k-','Linewidth',1.5);

                    %plot the Bayesian prediction
                    polarplot(deg2rad([bayesianPrediction,0,mod(bayesianPrediction+180,360)]),[max_rad 0 max_rad],'b:','Linewidth',1.5);

                    % plot title
                    title(['Motion Dir ' num2str(motion_dir)]);

                end

                % for each contrast level
                for c = 1:2

                    % get trial indices
                    ind = find(abs(data(:,1)-contrast(c))<0.001 & abs(data(:,2)-speed)<0.001 & ...
                        abs(data(:,3)-orientation)<0.001 & abs(data(:,4)-motion_dir)<0.001);
                    indCounter = [indCounter; ind];

                    % grab data just for these indicates
                    subdata = data(ind,:);

                    %remove no motion response trials
                    subdata(subdata(:,5)==0,:) = [];

                    % concatenate with other conditions for statistical analysis
                    subdataAll = [subdataAll; subdata];

                    % convert responses to angles in radians
                    thetaDeg = angles(subdata(:,5));
                    theta = deg2rad(thetaDeg);

                    % mirror angles for polar plot
                    theta2 = deg2rad(angles_mirror(subdata(:,5)));

                    % calculate error magnitude in radians
                    error = [];
                    for r = 1:length(theta)

                        % calculate circular distance from closest axis of ground truth motion

                        % positive = clockwise
                        % negative = counter clockwise

                        error_1 = circ_dist(deg2rad(motion_dir),theta(r));
                        error_2 = circ_dist(deg2rad(mod(motion_dir+180,360)),theta(r));

                        if abs(error_1) < abs(error_2)
                            error(r) = error_1;
                        else
                            error(r) = error_2;
                        end

                    end

                    %plot the data
                    if plotIndividuals
                        polarhistogram([theta theta2],'BinEdges',round(deg2rad(angles_edges),3),'FaceColor',ccolor(c,:),'FaceAlpha',0.3,'Normalization','probability');
                    end

                    % store these trials in a matrix with all subjects together
                    % column: subject, speed, motion dir, orientation, contrast, response, response mirrored, error, repeat number
                    repeats = 1:length(theta);
                    testConditionTrials = [testConditionTrials; [s.*ones([length(theta) 1]) speed.*ones([length(theta) 1]) motion_dir.*ones([length(theta) 1]) orientation.*ones([length(theta) 1]) c.*ones([length(theta) 1]) theta' theta2' error' repeats']];

                end

                % clean up plot
                if plotIndividuals
                    ax = gca; ax.FontSize = 8; rticklabels(ax, []);
                    thetaticklabels(ax,[]);
                    splotCntr  = splotCntr + 1;
                end

            end
        end

    end

    % save the plot
    if plotIndividuals
        print(gcf,['./plots/' grp '_' subjTags{s} '_polar_plots.pdf'],'-dpdf','-fillpage')
    end
end

%% plot average polar histogram across all subjects in this group

figure; hold on;
sgtitle(['Group ' grp '(blue = high contrast, red = low contrast)'])
set(gcf,'Position',[1 1 1400 1000]);

splotCntr = 1;

for speed = unique(testConditionTrials(:,2))' % for each speed
    for orientation = unique(testConditionTrials(:,4))' % for each orientation
        for motion_dir = [0 90 45 -45] % for each motion direction

            % get subplot
            subplot(4,4,splotCntr)

            % set radius maximum for plotting
            max_rad = 0.5;

            %plot the rhombus orientation
            polarplot(deg2rad([orientation orientation+90 orientation+180 orientation+270 orientation]), [0.4 0.07 0.4 0.07 0.4], 'color', [0.5 0.5 0.5], 'LineWidth', 1);  % Adjust the color and line properties as needed
            hold on;

            %plot the ground truth motion
            polarplot(deg2rad([motion_dir,0,mod(motion_dir+180,360)]),[max_rad 0 max_rad],'k-','Linewidth',1.5);

            %plot the Bayesian prediction
            if motion_dir==0 || motion_dir==90
                bayesianPrediction = -orientation;
                polarplot(deg2rad([bayesianPrediction,0,mod(bayesianPrediction+180,360)]),[max_rad 0 max_rad],'b:','Linewidth',1.5);
            end

            % set title
            title(['Motion Dir ' num2str(motion_dir)]);

            % for each contrast level
            for c = 1:2

                % get trial indices
                ind = find(abs(testConditionTrials(:,5)-c)<0.001 & abs(testConditionTrials(:,2)-speed)<0.001 & ...
                    abs(testConditionTrials(:,4)-orientation)<0.001 & abs(testConditionTrials(:,3)-motion_dir)<0.001);

                % grab data just for these indicates
                subdata = testConditionTrials(ind,:);

                %plot the data
                polarhistogram([subdata(:,6) subdata(:,7)],'BinEdges',round(deg2rad(angles_edges),3),'FaceColor',ccolor(c,:),'FaceAlpha',0.3,'Normalization','probability');

            end

            % plot clean up
            ax = gca; ax.FontSize = 8; rticklabels(ax, []);
            splotCntr  = splotCntr + 1;

        end
    end
end

print(gcf,['./plots/' grp '_ALL_polar_plots.pdf'],'-dpdf','-fillpage')

%% plot mean error to examine whether there is a systematic bias

% spacing between subjects for plotting
spacing = linspace(-0.25,0.25,5);

for motion_dir = [0 90 -45 45] % for each motion direction

    if motion_dir == 0
        figCard = figure; hold on;
        set(gcf,'Position',[1 1 650 200]);
        subplot(1,3,1); hold on;
        title('Horizontal motion')

        % holder for stats
        statsDataCardinal_lowC{s} = [];
        statsDataCardinal_highC{s} = [];

    elseif motion_dir == -45
        figDiag = figure; hold on;
        set(gcf,'Position',[1 1 650 200]);
        subplot(1,3,1); hold on;
        title('Short axis motion')

        % holder for stats
        statsDataDiagonal_lowC{s} = [];
        statsDataDiagonal_highC{s} = [];

    elseif motion_dir == 90
        figure(figCard)
        subplot(1,3,2); hold on;
        title('Vertical motion')
    elseif motion_dir == 45
        figure(figDiag)
        subplot(1,3,2); hold on;
        title('Long axis motion')
    end


    for s = 1:length(filenames) % for each subject

        for speed = unique(testConditionTrials(:,2))' % for each speed

            for c = 1:2 % for each contrast

                % column: subject, speed, motion dir, orientation, contrast, response, response mirrored, error

                % get relevant trials
                ind =   testConditionTrials(:,3) == motion_dir & ...
                    testConditionTrials(:,1) == s & ...
                    testConditionTrials(:,2)==speed & ...
                    testConditionTrials(:,5)==c;

                these_data          = testConditionTrials(ind,8);

                % calculate mean and stdev
                meanEst             = rad2deg(mean(these_data));
                stdEst              = rad2deg(std(these_data));

                % store data for statistical analysis
                % flip sign of two motions to align
                if c == 1
                    if motion_dir == 0
                        statsDataCardinal_lowC{s} = [statsDataCardinal_lowC{s}; these_data];
                    elseif motion_dir == 90
                        statsDataCardinal_lowC{s} = [statsDataCardinal_lowC{s}; -these_data];
                    elseif motion_dir == -45
                        statsDataDiagonal_lowC{s} = [statsDataDiagonal_lowC{s}; these_data];
                    elseif motion_dir == 45
                        statsDataDiagonal_lowC{s} = [statsDataDiagonal_lowC{s}; -these_data];
                    end
                elseif c == 2
                    if motion_dir == 0
                        statsDataCardinal_highC{s} = [statsDataCardinal_highC{s}; these_data];
                    elseif motion_dir == 90
                        statsDataCardinal_highC{s} = [statsDataCardinal_highC{s}; -these_data];
                    elseif motion_dir == -45
                        statsDataDiagonal_highC{s} = [statsDataDiagonal_highC{s}; these_data];
                    elseif motion_dir == 45
                        statsDataDiagonal_highC{s} = [statsDataDiagonal_highC{s}; -these_data];
                    end
                end

                speed_x = 2;

                if c == 1 % low contrast unfilled
                    h(1) = plot(speed_x+spacing(s),meanEst,'o','Color',ColorIt(s),'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',7);
                    errorbar(speed_x+spacing(s),meanEst,1.96*stdEst/sqrt(sum(ind)),'k','linewidth',0.5,'CapSize',0);
                else % high contrast filled
                    h(2) = plot(speed_x+spacing(s)+0.025,meanEst,'o','Color',ColorIt(s),'MarkerFaceColor',ColorIt(s),'LineWidth',1.5,'MarkerSize',7);
                    errorbar(speed_x+spacing(s)+0.025,meanEst,1.96*stdEst/sqrt(sum(ind)),'k','linewidth',0.5,'CapSize',0);
                end

            end
        end
    end
    set(gca,'FontSize',10);
    xlim([1.5 2.5]);
    set(gca,'xtick',[]);
    ylim([-60 60]);

    ylabel('Mean angular estimate (deg)');
    plot(xlim,[0 0],'k-');
    box on;

    if motion_dir == 0 || motion_dir == -45
        legend(h,'low contrast','high contrast','location','southeast')
    end

end

%% add combined absolute error plot and run stats

for s = 1:length(filenames) % for each subject

    % cardinals
    figure(figCard)
    subplot(1,3,3); hold on;
    title('Combined');

    meanEst_lowC = rad2deg(mean(statsDataCardinal_lowC{s}));
    stdEst_lowC  = rad2deg(std(statsDataCardinal_lowC{s}));

    meanEst_highC = rad2deg(mean(statsDataCardinal_highC{s}));
    stdEst_highC  = rad2deg(std(statsDataCardinal_highC{s}));

    plot(speed_x+spacing(s),meanEst_lowC,'o','Color',ColorIt(s),'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',7);
    errorbar(speed_x+spacing(s),meanEst_lowC,1.96*stdEst_lowC/sqrt(numel(statsDataCardinal_lowC{s})),'k','linewidth',0.5,'CapSize',0);

    plot(speed_x+spacing(s)+0.025,meanEst_highC,'o','Color',ColorIt(s),'MarkerFaceColor',ColorIt(s),'LineWidth',1.5,'MarkerSize',7);
    errorbar(speed_x+spacing(s)+0.025,meanEst_highC,1.96*stdEst_highC/sqrt(numel(statsDataCardinal_highC{s})),'k','linewidth',0.5,'CapSize',0);

    xlim([1.5 2.5]);
    set(gca,'xtick',[]);
    ylim([-60 60]);
    ylabel('Mean angular estimate (deg)');
    plot(xlim,[0 0],'k-');
    box on;
    set(gca,'FontSize',10);

    % stats
    [~,p_lowC,~,stats_lowC] = ttest(statsDataCardinal_lowC{s});
    [~,p_highC,~,stats_highC] = ttest(statsDataCardinal_highC{s});
    [~,p_highVlow,~,stats_highVlow] = ttest2(statsDataCardinal_highC{s},statsDataCardinal_lowC{s});

    % effect sizes
    d_lowC = meanEffectSize(statsDataCardinal_lowC{s},'effect','cohen');
    d_highC = meanEffectSize(statsDataCardinal_highC{s},'effect','cohen');
    d_highVlow = meanEffectSize(statsDataCardinal_highC{s},statsDataCardinal_lowC{s},'effect','cohen','paired','off','variancetype','unequal');

    fprintf('\n')
    fprintf('\n')
    disp(['subj ' num2str(s) ' STATS: Cardinals']);
    disp(['lowC vs 0: t(' num2str(stats_lowC.df) ') = ' num2str(stats_lowC.tstat,3) '; p = ' num2str(p_lowC,4) '; D = ' num2str(abs(d_lowC.Effect),2)])
    disp(['highC vs 0: t(' num2str(stats_highC.df) ') = ' num2str(stats_highC.tstat,3) '; p = ' num2str(p_highC,4) '; D = ' num2str(abs(d_highC.Effect),2)])
    disp(['highC vs lowC: t(' num2str(stats_highVlow.df) ') = ' num2str(stats_highVlow.tstat,3) '; p = ' num2str(p_highVlow,4) '; D = ' num2str(abs(d_highVlow.Effect),2)])

    % diagonals
    figure(figDiag)
    subplot(1,3,3); hold on;
    title('Combined');

    meanEst_lowC = rad2deg(circ_mean(statsDataDiagonal_lowC{s}));
    stdEst_lowC  = rad2deg(std(statsDataDiagonal_lowC{s}));

    meanEst_highC = rad2deg(circ_mean(statsDataDiagonal_highC{s}));
    stdEst_highC  = rad2deg(std(statsDataDiagonal_highC{s}));

    plot(speed_x+spacing(s),meanEst_lowC,'o','Color',ColorIt(s),'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',7);
    errorbar(speed_x+spacing(s),meanEst_lowC,1.96*stdEst_lowC/sqrt(numel(statsDataDiagonal_lowC{s})),'k','linewidth',0.5,'CapSize',0);

    plot(speed_x+spacing(s)+0.025,meanEst_highC,'o','Color',ColorIt(s),'MarkerFaceColor',ColorIt(s),'LineWidth',1.5,'MarkerSize',7);
    errorbar(speed_x+spacing(s)+0.025,meanEst_highC,1.96*stdEst_highC/sqrt(numel(statsDataDiagonal_highC{s})),'k','linewidth',0.5,'CapSize',0);

    xlim([1.5 2.5]);
    set(gca,'xtick',[]);
    ylim([-60 60]);
    ylabel('Mean angular estimate (deg)');
    plot(xlim,[0 0],'k-');
    box on;
    set(gca,'FontSize',10);

    % stats
    [~,p_lowC,~,stats_lowC] = ttest(statsDataDiagonal_lowC{s});
    [~,p_highC,~,stats_highC] = ttest(statsDataDiagonal_highC{s});
    [~,p_highVlow,~,stats_highVlow] = ttest2(statsDataDiagonal_highC{s},statsDataDiagonal_lowC{s});

    % effect sizes
    d_lowC = meanEffectSize(statsDataDiagonal_lowC{s},'effect','cohen');
    d_highC = meanEffectSize(statsDataDiagonal_highC{s},'effect','cohen');
    d_highVlow = meanEffectSize(statsDataDiagonal_highC{s},statsDataCardinal_lowC{s},'effect','cohen','paired','off','variancetype','unequal');

    fprintf('\n')
    disp(['subj ' num2str(s) ' STATS: Diagonals']);
    disp(['lowC vs 0: t(' num2str(stats_lowC.df) ') = ' num2str(stats_lowC.tstat,3) '; p = ' num2str(p_lowC,4) '; D = ' num2str(abs(d_lowC.Effect),2)])
    disp(['highC vs 0: t(' num2str(stats_highC.df) ') = ' num2str(stats_highC.tstat,3) '; p = ' num2str(p_highC,4) '; D = ' num2str(abs(d_highC.Effect),2)])
    disp(['highC vs lowC: t(' num2str(stats_highVlow.df) ') = ' num2str(stats_highVlow.tstat,3) '; p = ' num2str(p_highVlow,4) '; D = ' num2str(abs(d_highVlow.Effect),2)])

end

% save plots
saveas(figCard,['./plots/' grp '_mean_error_cardinal_orientations.png'])
saveas(figDiag,['./plots/' grp '_mean_error_diagonal_orientations.png'])
