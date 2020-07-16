%%% Scritps for plotting the Figure 2 graphics of Hippocampal-Prefrontal
%%% Interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

load('A_data.mat')

% plotAuxEx - Mean binned (time and position) coherogram data
% F - Vector of frequencies

figure()
imagesc(1:30,F,plotAux')
axis xy
ylim([5 20])
set(gca,'fontsize',14)
xticks(1:1:30)
xticklabels(string(1:0.1:2))
xticklabels([-5:-1 string(1:0.05:1.95) 1:5])
xlabel('Linearized Position')
title('CA1 - PFC coherence','fontsize',15)
ylabel('Frequency (Hz)')
colormap parula
colorbar
line([5.5 5.5],[0 20],'Color', 'k', 'LineWidth', 3); % Time bins division
line([25.5 25.5],[0 20],'Color', 'k', 'LineWidth', 3); % Time bins division

%% B

load('B_data.mat')

% plotAux_l - Mean binned spectrogram data for left trials, 2 is PFC
% plotAux_r - Mean binned spectrogram data for right trials, 2 is PFC

% Making custom colormaps from mean power (theta or delta) values
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

for trial = [12 18]

    if(trial ==18)
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


%% C

load('C_data.mat')

% avg_theta_pow - Average theta power in referred linearized positions
% avg_theta_ppow - Avg peak theta power in referred linPos
% avg_theta_coh - Avg theta coherence in referred linPos
% avg_theta_pcoh - Avg peak theta coherence in referred linPos

figure
subplot(221)
boxplot([avg_theta_pow(:,1),avg_theta_pow(:,2)])

subplot(222)
boxplot([avg_theta_coh(:,1),avg_theta_coh(:,2)])

subplot(223)
boxplot([avg_theta_ppow(:,1),avg_theta_ppow(:,2)])

subplot(224)
boxplot([avg_theta_pcoh(:,1),avg_theta_pcoh(:,2)])


%% D

load('D_data.mat')

% mTaskEv - Mean task events coherence
% F - Vector of frequencies
% Tcoh - Time vector
% p - Time before and after event points (250ms)

eventlabel{1} = 'Trial start';
eventlabel{2} = 'Choice Point';
eventlabel{3} = 'Turn';
eventlabel{4} = 'Reward';

figure()
for event = 1:4
        subplot(1,4,event)
        imagesc(Tcoh-p,F,squeeze(mTaskEv(event,:,:))')
        
        title(eventlabel{event}, 'fontsize',13)
        set(gca,'fontsize',12)
        xlabel('Time (s)')
        if event == 1
            ylabel('Frequency (Hz)')
        else
            set(gca,'ytick',[])
        end
        axis xy
        caxis([min(min(min(mTaskEv(:,:,F>0 & F<20))))...
            max(max(max(mTaskEv(:,:,F>0 & F<20))))])
        ylim([5 20])
end

colorbar

%% E

load('E_data.mat')

% ThetaCoh - Mean theta coherence in each task event

eventlabel{1} = 'Trial start';
eventlabel{2} = 'Choice Point';
eventlabel{3} = 'Turn';
eventlabel{4} = 'Reward';

figure()
plot(mean(ThetaCoh))
xticklabels(eventlabel)
hold on
errorbar(mean(ThetaCoh),std(ThetaCoh)/sqrt(count))
hold off
xlim([0 5])


