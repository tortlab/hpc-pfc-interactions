<<<<<<< HEAD
%%% Scritps for plotting the Figure 1 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% B

load('B_data.mat')

% plotAuxEx - Mean binned (time and position) spectrogram data, 2 is PFC
% FA - Vector of frequencies for the spectrograms
% msXCombined - Mean speeds for time and position bins

% HPC spectrogram
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
imagesc(1:30,FA, plotAuxEx)
set(gca,'fontsize',14)
xticks(1:1:30)
xticklabels(string(1:0.1:2))
xticklabels([-5:-1 string(1:0.05:1.95) 1:5])
set(gca,'Ydir','normal')
title('CA1','fontsize',13)
ylim([0 20])
xlabel('Linearized Position')
ylabel('Frequency (Hz)')
colormap parula
colorbar
line([5.5 5.5],[0 20],'Color', 'k', 'LineWidth', 3);
line([25.5 25.5],[0 20],'Color', 'k', 'LineWidth', 3);
yyaxis right
set(gca,'YColor', 'r')
ylim([1 20])
hold on
plot(smooth(mean(msXCombined)), 'r', 'LineWidth', 4)
yticklabels(string(yticks*10))
ylabel('Speed cm/s', 'Color', 'r')
hold off

% PFC spectrogram
subplot(2,1,2)
imagesc(1:30,FA, plotAuxEx2)
set(gca,'fontsize',14)
title('PFC','fontsize',13)
xticks(1:1:30)
xticklabels(string(1:0.1:2))
xticklabels([-5:-1 string(1:0.05:1.95) 1:5])
ylim([0 20])
set(gca,'Ydir','normal')
xlabel('Linearized Position')
ylabel('Frequency (Hz)')
colormap parula
colorbar
line([5.5 5.5],[0 20],'Color', 'k', 'LineWidth', 3);
line([25.5 25.5],[0 20],'Color', 'k', 'LineWidth', 3);
yyaxis right
set(gca,'YColor', 'r')
ylim([1 20])
hold on
plot(smooth(mean(msXCombined)), 'r', 'LineWidth', 4)
yticklabels(string(yticks*10))
ylabel('Speed cm/s', 'Color', 'r')
hold off

%% C

load('C_data.mat')

% plotAux_l - Mean binned spectrogram data for left trials, 2 is PFC
% plotAux_r - Mean binned spectrogram data for right trials, 2 is PFC


colorL = mean(plotAux_l(:,F>6 & F<10),2); % This line for HPC
colorR = mean(plotAux_r(:,F>6 & F<10),2); % This line for HPC

% colorL = mean(plotAux_l2(:,F>0.5 & F<4),2); % This line for PFC
% colorR = mean(plotAux_r2(:,F>0.5 & F<4),2); % This line for PFC

% Making custom colormaps from mean power (theta or delta) values
colL = colorL/max(colorL);
[~, idxL] = sort(colL);
[~, idxL] = sort(idxL);
cL = parula(length(colL));

colR = colorR/max(colorR);
[~, idxR] = sort(colR);
[~, idxR] = sort(idxR);
cR = parula(length(colR));


for trial = [1 2]
    if(trial == 2)
        
        whlrl_trial(:,2) = whlrl_trial(:,2)+40;
        bin=1;
        
        for i=1:0.05:1.95
            cond = find(whlrl_trial(:,7)>i & whlrl_trial(:,7)<i+0.05,1,...
                'first');
            plot(whlrl_trial(cond,1),whlrl_trial(cond,2),'o',...
                'MarkerFaceColor',[cR(idxR(bin),1), cR(idxR(bin),2), ...
                cR(idxR(bin),3)],'MarkerEdgeColor','k','MarkerSize',50)
            hold on
            bin = bin+1;
        end
        
    else
        
        whlrl_trial = removerows(whlrl_trial,whlrl_trial(:,7)<1.0382);
        bin=1;
        
        for i=1:0.05:1.95
            cond = find(whlrl_trial(:,7)>i & whlrl_trial(:,7)<i+0.05,1,...
                'first');
            plot(whlrl_trial(cond,1),whlrl_trial(cond,2),'o',...
                'MarkerFaceColor',[cL(idxL(bin),1), cL(idxL(bin),2), ...
                cL(idxL(bin),3)],'MarkerEdgeColor','k','MarkerSize',50)
            hold on
            bin = bin+1;
        end
        
    end
end

colormap parula
h = colorbar
ylabel(h,'Power')
caxis([min(color) max(color)])
set(gca, 'Ydir', 'normal')
set(gca,'xtick',[],'ytick',[])
ylim([0 500])
xlim([150 650])
axis square
hold off
camroll(-90)

%% D

load('D_data.mat')

figure()
set(gca,'fontsize',20)
ylabel('Power')
xticklabels([{'Choice Point'}, {'Turn Start'}, {'Turn End'}])
boxplot([nanmean(cptrials,2), nanmean(tstrials,2), nanmean(tetrials,2)])







=======
%%% Scritps for plotting the Figure 1 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% B

load('B_data.mat')

% plotAuxEx - Mean binned (time and position) spectrogram data, 2 is PFC
% FA - Vector of frequencies for the spectrograms
% msXCombined - Mean speeds for time and position bins

% HPC spectrogram
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
imagesc(1:30,FA, plotAuxEx)
set(gca,'fontsize',14)
xticks(1:1:30)
xticklabels(string(1:0.1:2))
xticklabels([-5:-1 string(1:0.05:1.95) 1:5])
set(gca,'Ydir','normal')
title('CA1','fontsize',13)
ylim([0 20])
xlabel('Linearized Position')
ylabel('Frequency (Hz)')
colormap parula
colorbar
line([5.5 5.5],[0 20],'Color', 'k', 'LineWidth', 3);
line([25.5 25.5],[0 20],'Color', 'k', 'LineWidth', 3);
yyaxis right
set(gca,'YColor', 'r')
ylim([1 20])
hold on
plot(smooth(mean(msXCombined)), 'r', 'LineWidth', 4)
yticklabels(string(yticks*10))
ylabel('Speed cm/s', 'Color', 'r')
hold off

% PFC spectrogram
subplot(2,1,2)
imagesc(1:30,FA, plotAuxEx2)
set(gca,'fontsize',14)
title('PFC','fontsize',13)
xticks(1:1:30)
xticklabels(string(1:0.1:2))
xticklabels([-5:-1 string(1:0.05:1.95) 1:5])
ylim([0 20])
set(gca,'Ydir','normal')
xlabel('Linearized Position')
ylabel('Frequency (Hz)')
colormap parula
colorbar
line([5.5 5.5],[0 20],'Color', 'k', 'LineWidth', 3);
line([25.5 25.5],[0 20],'Color', 'k', 'LineWidth', 3);
yyaxis right
set(gca,'YColor', 'r')
ylim([1 20])
hold on
plot(smooth(mean(msXCombined)), 'r', 'LineWidth', 4)
yticklabels(string(yticks*10))
ylabel('Speed cm/s', 'Color', 'r')
hold off

%% C

load('C_data.mat')

% plotAux_l - Mean binned spectrogram data for left trials, 2 is PFC
% plotAux_r - Mean binned spectrogram data for right trials, 2 is PFC


colorL = mean(plotAux_l(:,F>6 & F<10),2); % This line for HPC
colorR = mean(plotAux_r(:,F>6 & F<10),2); % This line for HPC

% colorL = mean(plotAux_l2(:,F>0.5 & F<4),2); % This line for PFC
% colorR = mean(plotAux_r2(:,F>0.5 & F<4),2); % This line for PFC

% Making custom colormaps from mean power (theta or delta) values
colL = colorL/max(colorL);
[~, idxL] = sort(colL);
[~, idxL] = sort(idxL);
cL = parula(length(colL));

colR = colorR/max(colorR);
[~, idxR] = sort(colR);
[~, idxR] = sort(idxR);
cR = parula(length(colR));


for trial = [1 2]
    if(trial == 2)
        
        whlrl_trial(:,2) = whlrl_trial(:,2)+40;
        bin=1;
        
        for i=1:0.05:1.95
            cond = find(whlrl_trial(:,7)>i & whlrl_trial(:,7)<i+0.05,1,...
                'first');
            plot(whlrl_trial(cond,1),whlrl_trial(cond,2),'o',...
                'MarkerFaceColor',[cR(idxR(bin),1), cR(idxR(bin),2), ...
                cR(idxR(bin),3)],'MarkerEdgeColor','k','MarkerSize',50)
            hold on
            bin = bin+1;
        end
        
    else
        
        whlrl_trial = removerows(whlrl_trial,whlrl_trial(:,7)<1.0382);
        bin=1;
        
        for i=1:0.05:1.95
            cond = find(whlrl_trial(:,7)>i & whlrl_trial(:,7)<i+0.05,1,...
                'first');
            plot(whlrl_trial(cond,1),whlrl_trial(cond,2),'o',...
                'MarkerFaceColor',[cL(idxL(bin),1), cL(idxL(bin),2), ...
                cL(idxL(bin),3)],'MarkerEdgeColor','k','MarkerSize',50)
            hold on
            bin = bin+1;
        end
        
    end
end

colormap parula
h = colorbar
ylabel(h,'Power')
caxis([min(color) max(color)])
set(gca, 'Ydir', 'normal')
set(gca,'xtick',[],'ytick',[])
ylim([0 500])
xlim([150 650])
axis square
hold off
camroll(-90)

%% D

load('D_data.mat')

figure()
set(gca,'fontsize',20)
ylabel('Power')
xticklabels([{'Choice Point'}, {'Turn Start'}, {'Turn End'}])
boxplot([nanmean(cptrials,2), nanmean(tstrials,2), nanmean(tetrials,2)])







>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
