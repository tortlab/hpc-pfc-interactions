%%% Scritps for plotting the Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort
%% A

load('A_data.mat')

% PhaseFreqVector - Frequency vector for phase-extracted oscillations
% PhaseFreq_BandWidth - Size of bandwith for phase vector above
% AmpFreqVector - Frequency vector for amplitude-extracted oscillations
% AmpFreq_BandWidth - Size of bandwith for amplitude vector above
% Comodulogram - Matrix of comodulation values

figure
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,...
    AmpFreqVector+AmpFreq_BandWidth/2,...
    squeeze(Comodulogram(:,:))',30,'lines','none')
set(gca,'fontsize',12)
xlabel('Phase Freq (Hz)')
a = xlim;
b = ylim;
xlim([a(1) 20])
ylim([b(1) 100])
ylabel('Amplitude Frequency (Hz)')
h = colorbar
clims = caxis()
caxis([0.0001 clims(2)])
set(h,'Ticks',[0.0001 clims(2)])
xticks([2,5,10, 15, 20])
yticks([20, 40, 60, 80, 100])
set(gca,'fontsize',20)

%% B

load('B_data.mat')

% nbins - Number of spatial bins
% AmpFreqVector - Frequency vector for amplitude-extracted oscillations
% AmpFreq_BandWidth - Size of bandwith for amplitude vector above
% FreqBinComod - Matrix of comodulation values for spatial bins

figure
contourf(1:nbins,...
    AmpFreqVector+AmpFreq_BandWidth/2,FreqBinComod',30,'lines','none')
axis xy
set(gca,'fontsize',12)
xlabel('Bin')
xticks([1 2 3 4 5])
ylabel('Amplitude Frequency (Hz)')
h = colorbar
clims = caxis()
caxis([0.0005 clims(2)])
set(h,'Ticks',[0.0005 clims(2)])



%% C

load('C_data.mat')

% ca1P_ca1A - Matrix of comodulation values from ca1 phases and ca1 amps
% pfcP_pfcA - Matrix of comodulation values from pfc phases and pfc amps
% ca1P_pfcA - Matrix of comodulation values from ca1 phases and pfc amps
% pfcP_ca1A - Matrix of comodulation values from pfc phases and ca1 amps

comods(1,:,:) = ca1P_ca1A;
comods(2,:,:) = pfcP_pfcA;
comods(3,:,:) = ca1P_pfcA;
comods(4,:,:) = pfcP_ca1A;

figure
for i = 1:4
    clear comod
    comod = squeeze(comods(i,:,:));
    comod = zscore(comod,0,2);
    
    subplot(1,4,i)
    errorbar(mean(comod),...
        std(comod)/sqrt(13),'bo-','markerf','w')
    set(gca,'xtick',1:5)
    xlim([0.5 5.5])
    ylim([-1.15 1.3])
    axis square
    box off
end