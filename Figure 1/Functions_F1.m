%%% Scritps for generating the data for Figure 1 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% B

plotAux = [];
plotAux2 = [];
linPos_Power(30).pow =[];
linPos_Power2(30).pow =[];

for sess=1:13
    [srate_lfp, dt, lfpPFC, lfpCA1, srate_bhv, whlrld, whlrl_speed,...
        SessionNP] = LoadSessionData(paths{sess});
    
    n_trials = length(unique(whlrld(:,6)))-1;
    
    for trial=1:n_trials
        
        ts = SessionNP(trial,2); % last nose-poking
        p = (find(whlrld(:,6)==trial,1,'first')/srate_bhv);
        p2 = (find(whlrld(:,6)==trial,1,'last')/srate_bhv);
        k = 5*srate_lfp;
        period = round(p*srate_lfp-k:p2*srate_lfp+k);
        
        ts = SessionNP(trial,2) - SessionNP(trial,1)-k/srate_lfp;
        te = SessionNP(trial,3) - SessionNP(trial,1)+k/srate_lfp;
        
        lfpCA1t=lfpCA1(:,period);
        lfpPFCt=lfpPFC(:,period);
        
        % Spectrogram parameters
        WindowLength = 0.5; % Seconds
        WindowLength = WindowLength*srate_lfp;
        Overlap = round(0.45*srate_lfp);
        NFFT = 2^13;
        
        [SA, FA, TA, PA] = spectrogram(lfpCA1t(chCA1(sess),:),...
            WindowLength,Overlap,NFFT, srate_lfp);
        
        [SA2, FA2, TA2, PA2] = spectrogram(lfpPFCt(chPFC,:),...
            WindowLength,Overlap,NFFT, srate_lfp);
  
        k = (k/srate_lfp)*srate_bhv;
        aux = whlrld(round(p*srate_bhv+TA*srate_bhv-k),7);
        tmp = 1:0.05:2;
        
        for posi = 1:5
            ind = [ts-6+posi ts-5+posi];
            linPos_Power(posi).pow = horzcat(linPos_Power(posi).pow, ...
                PA(:,TA>ind(1)&TA<ind(2)));
            linPos_Power2(posi).pow = horzcat(linPos_Power2(posi).pow, ...
                PA2(:,TA>ind(1)&TA<ind(2)));
        end
        
        for posi=6:25
            linPos_Power(posi).pow = horzcat(linPos_Power(posi).pow, ...
                PA(:,aux > tmp(posi-5) & aux < tmp(posi-4)));
            linPos_Power2(posi).pow = horzcat(linPos_Power2(posi).pow, ...
                PA2(:,aux > tmp(posi-5) & aux < tmp(posi-4)));
        end
        
       
        for posi = 26:30
            ind = [te+posi-26 te+posi-25];
            linPos_Power(posi).pow = horzcat(linPos_Power(posi).pow, ...
                PA(:,TA>ind(1)&TA<ind(2)));
            linPos_Power2(posi).pow = horzcat(linPos_Power2(posi).pow, ...
                PA2(:,TA>ind(1)&TA<ind(2)));
        end
        
    end
 
end


for posi=1:30
        linPos_Power(posi).pow = mean(linPos_Power(posi).pow,2);
        linPos_Power2(posi).pow = mean(linPos_Power2(posi).pow,2);
        plotAux = horzcat(plotAux, linPos_Power(posi).pow);
        plotAux2 = horzcat(plotAux2, linPos_Power2(posi).pow);
end

%% D


% controlData4 =
% 1 = theta (6-10)
% 2 = theta/delta
% 3 = speed
% 4 = delta (0.5 - 4)

cp = 3:5;
ts = 9:11;
te = 12:14;


rat = 1:508;

while(p<0.8)
clear cptrials tstrials tetrials

allowSpeed = [4 5];

cptrials(:,1) = nanmean(controlData3(nanmean(controlData3(rat,cp,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,cp,3),2)<=allowSpeed(2),cp,3),2);
cptrials(:,2) = nanmean(controlData3(nanmean(controlData3(rat,cp,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,cp,3),2)<=allowSpeed(2),cp,1),2);
cptrials(:,3) = nanmean(controlData3(nanmean(controlData3(rat,cp,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,cp,3),2)<=allowSpeed(2),cp,4),2);


tstrials(:,1) = nanmean(controlData3(nanmean(controlData3(rat,ts,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,ts,3),2)<=allowSpeed(2),ts,3),2);
tstrials(:,2) = nanmean(controlData3(nanmean(controlData3(rat,ts,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,ts,3),2)<=allowSpeed(2),ts,1),2);
tstrials(:,3) = nanmean(controlData3(nanmean(controlData3(rat,ts,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,ts,3),2)<=allowSpeed(2),ts,4),2);

tetrials(:,1) = nanmean(controlData3(nanmean(controlData3(rat,te,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,te,3),2)<=allowSpeed(2),te,3),2);
tetrials(:,2) = nanmean(controlData3(nanmean(controlData3(rat,te,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,te,3),2)<=allowSpeed(2),te,1),2);
tetrials(:,3) = nanmean(controlData3(nanmean(controlData3(rat,te,3),2)>=allowSpeed(1)...
    & nanmean(controlData3(rat,te,3),2)<=allowSpeed(2),te,4),2);


size = min([length(cptrials), length(tstrials), length(tetrials)]);
cptrials = datasample(cptrials, size, 'Replace', false);
tstrials = datasample(tstrials, size, 'Replace', false);
tetrials = datasample(tetrials, size, 'Replace', false);

clear s_comp
s_comp = [nanmean(cptrials(:,1),2), nanmean(tstrials(:,1),2),...
    nanmean(tetrials(:,1),2)];

clear t_comp
t_comp = [nanmean(cptrials(:,2),2), nanmean(tstrials(:,2),2),...
    nanmean(tetrials(:,2),2)];

clear d_comp
d_comp = [nanmean(cptrials(:,3),2), nanmean(tstrials(:,3),2),...
    nanmean(tetrials(:,3),2)];

p = kruskalwallis(s_comp);

end
