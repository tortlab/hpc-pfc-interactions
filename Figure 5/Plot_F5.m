<<<<<<< HEAD
%%% Scritps for plotting the Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort
%%

load('F5_data.mat')

% slow_vector - Frequency vector of slow oscillations
% slow_BandWidth - Size of bandwidth for the slow vector
% fast_vector - Frequency vector of fast oscillations
% fast_BandWidth - Size of bandwidth for the fast vector
% pfcSW_m, pfcSW_n - Actual and surrogate distributions for PFC slow waves
% ca1SW_m, ca1SW_n - Actual and surrogate distributions for HPC slow waves

figure
subplot(2,1,1)
contourf(slow_vector+slow_BandWidth/2,fast_vector+fast_BandWidth/2,...
    dist,50,'edgecolor','none')
set(gca,'Fontsize',12)
ylabel('PLV Frequency (Hz)','Fontsize',20)
xlabel('Phase Frequency (Hz)','Fontsize',20)
h = colorbar;
ylabel(h, 'Mod Index','Fontsize',15)

subplot(2,1,2)
contourf(slow_vector+slow_BandWidth/2,fast_vector+fast_BandWidth/2,...
    pfcSW,50,'edgecolor','none')
set(gca,'Fontsize',12)
ylabel('PLV Frequency (Hz)','Fontsize',20)
xlabel('Phase Frequency (Hz)','Fontsize',20)
h = colorbar;
ylabel(h, 'Mod Index','Fontsize',15)


figure
subplot(2,1,1)
boxplot([pfcSW_m' pfcSW_n'])
subplot(2,1,2)
boxplot([ca1SW_m' ca1SW_n'])





=======
%%% Scritps for plotting the Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort
%%

load('F5_data.mat')

% slow_vector - Frequency vector of slow oscillations
% slow_BandWidth - Size of bandwidth for the slow vector
% fast_vector - Frequency vector of fast oscillations
% fast_BandWidth - Size of bandwidth for the fast vector
% pfcSW_m, pfcSW_n - Actual and surrogate distributions for PFC slow waves
% ca1SW_m, ca1SW_n - Actual and surrogate distributions for HPC slow waves

figure
subplot(2,1,1)
contourf(slow_vector+slow_BandWidth/2,fast_vector+fast_BandWidth/2,...
    dist,50,'edgecolor','none')
set(gca,'Fontsize',12)
ylabel('PLV Frequency (Hz)','Fontsize',20)
xlabel('Phase Frequency (Hz)','Fontsize',20)
h = colorbar;
ylabel(h, 'Mod Index','Fontsize',15)

subplot(2,1,2)
contourf(slow_vector+slow_BandWidth/2,fast_vector+fast_BandWidth/2,...
    pfcSW,50,'edgecolor','none')
set(gca,'Fontsize',12)
ylabel('PLV Frequency (Hz)','Fontsize',20)
xlabel('Phase Frequency (Hz)','Fontsize',20)
h = colorbar;
ylabel(h, 'Mod Index','Fontsize',15)


figure
subplot(2,1,1)
boxplot([pfcSW_m' pfcSW_n'])
subplot(2,1,2)
boxplot([ca1SW_m' ca1SW_n'])





>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
