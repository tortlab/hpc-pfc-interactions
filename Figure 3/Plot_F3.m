%%% Scritps for plotting the Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

load('A_data.mat')

bound = squeeze(SEM_hpc_ts)';
bound2 = squeeze(SEM_pfc_ts)';

figure()
subplot(1,4,1)
boundedline(freqs,squeeze(ts_GC(1,:)),bound,'-b')
hold on
boundedline(freqs,squeeze(ts_GC(2,:)),bound2,'-r')
ylim([0 maxi])
box off
title('Start')
ylabel('Granger Causality')
xlabel('Frequency (Hz)')
hold off
xlim([0 40])

bound = squeeze(SEM_hpc_cp)';
bound2 = squeeze(SEM_pfc_cp)';

subplot(1,4,2)
boundedline(freqs,squeeze(cp_GC(1,:)),bound,'-b')
hold on
boundedline(freqs,squeeze(cp_GC(2,:)),bound2,'-r')
ylim([0 maxi])
box off
title('Choice Point')
ylabel('Granger Causality')
xlabel('Frequency (Hz)')
hold off
xlim([0 40])

bound = squeeze(SEM_hpc_tu)';
bound2 = squeeze(SEM_pfc_tu)';

subplot(1,4,3)
boundedline(freqs,squeeze(tu_GC(1,:)),bound,'-b')
hold on
boundedline(freqs,squeeze(tu_GC(2,:)),bound2,'-r')
ylim([0 maxi])
box off
title('Turn')
ylabel('Granger Causality')
xlabel('Frequency (Hz)')
hold off
xlim([0 40])

bound = squeeze(SEM_hpc_re)';
bound2 = squeeze(SEM_pfc_re)';

subplot(1,4,4)
boundedline(freqs,squeeze(re_GC(1,:)),bound,'-b')
hold on
boundedline(freqs,squeeze(re_GC(2,:)),bound2,'-r')
ylim([0 maxi])
box off
title('Reward')
ylabel('Granger Causality')
xlabel('Frequency (Hz)')
hold off
xlim([0 40])
%% B

load('B_data.mat')

figure
subplot(2,1,1)
imagesc(1:10, freqs, squeeze(F_mean(:,1,:))')
axis xy
xlabel("Bin")
ylabel("Frequency (Hz)")
title("HPC - PFC")
h = colorbar;
title(h,'Granger Causality')
ylim([0 20])

subplot(2,1,2)
imagesc(1:10, freqs, squeeze(F_mean(:,2,:))')
axis xy
xlabel("Bin")
ylabel("Frequency (Hz)")
title("PFC - HPC")
h = colorbar;
title(h,'Granger Causality')
ylim([0 20])

%% C

load('C_data.mat')

colorL = nanmean(hpc_gcL);
colorR = nanmean(hpc_gcR);

% Making custom colormaps from mean GC (theta or delta) values
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
        
        for i=1:0.1:1.95
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
        
        for i=1:0.1:1.95
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
set(h, 'YTick',caxis())
set(gca, 'Ydir', 'normal')
set(gca,'xtick',[],'ytick',[])
ylim([0 500])
xlim([150 650])
axis square
hold off
camroll(-90)