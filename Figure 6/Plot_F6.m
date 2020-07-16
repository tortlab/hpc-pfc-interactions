%%% Scritps for plotting the Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort
%% A

load('A_data.mat')

% xHH - Avg MVL values for HPC LFP - HPC spikes
% xPP - Avg MVL values for PFC LFP - PFC spikes
% xPH - Avg MVL values for PFC LFP - HPC spikes
% xHP - Avg MVL values for HPC LFP - PFC spikes
% eHH, ePP, ePH, eHP - Errorbars for the above MVL combinations

figure
subplot(211)
boundedline(1:50,xHH(1:50),eHH(1:50), 'b', ...
    1:50,xPH(1:50),ePH(1:50),'r', 'alpha')
set(gca,'fontsize',20)
xlabel('Frequency (Hz)')
ylabel('MVL')
title('HPC Cells')
legend('HPC->HPC', 'PFC->HPC')
yl = ylim();
xlim([0 30])

subplot(212)
boundedline(1:50,xHP(1:50),eHP(1:50), 'b', ...
    1:50,xPP(1:50),ePP(1:50),'r', 'alpha')
set(gca,'fontsize',20)
xlabel('Frequency (Hz)')
ylabel('MVL')
title('PFC Cells')
legend('HPC->PFC', 'PFC->PFC')
xlim([0 30])

% A2
figure
bar(1:4,pkTheta, 'edgecolor', [0, 0, 0.2], 'linew', 1.3)
hold on
errorbar(1:4,pkTheta,SEM,'.','color','k', 'linew', 1.5)
hold off
ylabel('Peak theta MVL')
xticklabels({'HPC-HPC', 'PFC-HPC', 'HPC-PFC', 'PFC-PFC'})
set(gca, 'fontsize', 15)


%% B

load('B_data.mat')

% poolHH - Mean MVL on spatial bins for HPC LFP - HPC spikes
% poolPP - Mean MVL on spatial bins for PFC LFP - PFC spikes
% poolPH - Mean MVL on spatial bins for PFC LFP - HPC spikes
% poolHP - Mean MVL on spatial bins for HPC LFP - PFC spikes

comb_labels = {'PFC->PFC', 'PFC->HPC', 'HPC->PFC', 'HPC->HPC'};

figure
for comb=1:4
    
    subplot(1,4,comb)
    imagesc(1:5,1:50,squeeze(nanmean(pool)))
    xlabel('Bin')
    ylabel('Frequency (Hz)')
    axis xy
    xticks(1:5)
    set(gca, 'fontsize', 15)
    cb = colorbar
    set(cb, 'YTick',caxis())
    axis square
    title(comb_labels{comb})
    
end


%% C

load('C_data.mat')

% poolHH - Mean theta MVL on spatial bins for HPC LFP - HPC spikes
% poolPP - Mean delta MVL on spatial bins for PFC LFP - PFC spikes
% poolPH - Mean theta MVL on spatial bins for PFC LFP - HPC spikes
% poolHP - Mean theta MVL on spatial bins for HPC LFP - PFC spikes

plotter{1} = poolHH;
plotter{2} = poolPP;
plotter{3} = poolPH;
plotter{4} = poolHP;

figure
for j = 1:4
    subplot(1,4,j)
    errorbar(mean(plotter{j}),...
        std(plotter{j})/sqrt(length(plotter{j})),'bo-','markerf','w')
    set(gca,'xtick',1:5)
    xlim([0.5 5.5])
    axis square
    box off
    yl = ylim();
    ym = (yl(1)+yl(2))/2;
    yticks([yl(1), ym, yl(2)])
end
axis square