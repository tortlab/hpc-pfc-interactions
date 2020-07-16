<<<<<<< HEAD
%%% Scritps for generating the data for Figure 2 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

k = 7;
plotAux = [];
clear linPos_Coh
linPos_Coh(30).coh =[];
linPosCoh = zeros(13,30,5001);

for sess=1:13
    
    [srate_lfp, dt, lfpPFC, lfpCA1, srate_bhv, whlrld, whlrl_speed,...
        SessionNP] = LoadSessionData(paths{sess});
    
    plotAux = [];
    clear linPos_Coh
    linPos_Coh(30).coh =[];
    
    if(sess>10)
        chPFC =11;
    end
    
    for trial=1:trials(sess)
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
%         if(whlrld(ts_idx,5)~=2) % 1 = R, 2 = L
%             continue
%         end

        totTrial = sum(trials(1:sess-1)) + trial
        [ids{sess} ' ' num2str(trial)]
        start_time = (SessionNP(trial,1)-k)*srate_lfp;
        end_time = (SessionNP(trial,3)+k)*srate_lfp;
        ts = SessionNP(trial,2) - SessionNP(trial,1)+k;
        te = SessionNP(trial,3) - SessionNP(trial,1)+k;
      
        
        lfpCA1t=lfpCA1(:,round(start_time):round(end_time));
        lfpPFCt=lfpPFC(:,round(start_time):round(end_time));
        
       
        windowlength = 4*srate_lfp;
        stepsize = round(srate_lfp/8);
        
        Nwindow = (length(lfpCA1t)-windowlength)/stepsize+1;
        
        clear Cxy_T T
        for nwin = 1:Nwindow
            win = (1:windowlength) + (nwin-1)* stepsize;
            [Cxy, F]=mscohere(lfpCA1t(chCA1(sess),win),lfpPFCt(chPFC,win),1*srate_lfp,...
                [],10000,srate_lfp);
            Cxy_T(nwin,:)=Cxy;
            T(nwin) = median(win)*dt;
        end
       
        p  = round(SessionNP(trial,1)*srate_bhv);
        p2  = round(SessionNP(trial,3)*srate_bhv);
        aux = whlrld(round(p:4.875:p2),7);
        tmp = 1:0.05:2;
        
        Cxy_trial = Cxy_T(find(T<k,1,'last')+1:find(T>te,1,'first'),:);
        if size(Cxy_trial,1)<size(aux,1)
           aux = removerows(aux,1:size(aux,1)-size(Cxy_trial,1));
        end
        
        for posi = 1:5
            ind = [ts-6+posi ts-5+posi];
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_T(T>ind(1)&T<ind(2),:));
        end
        
        for posi = 6:25
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_trial(aux > tmp(posi-5) & aux < tmp(posi-4),:));              
        end
        
        for posi = 26:30
            ind = [te+posi-26 te+posi-25];
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_T(T>ind(1)&T<ind(2),:));
        end
        
    end
    for j=1:30
        if(size(linPos_Coh(j).coh,1)>1)
            linPos_Coh(j).coh = mean(linPos_Coh(j).coh);
        end
        plotAux = vertcat(plotAux, linPos_Coh(j).coh);
    end
    linPosCoh(sess,:,:) = plotAux;
end


plotAux = [];

for j=1:30
    if(size(linPos_Coh(j).coh,1)>1)
        linPos_Coh(j).coh = mean(linPos_Coh(j).coh);
    end
    plotAux = vertcat(plotAux, linPos_Coh(j).coh);
end

%% C

before_bif = 4:6;
after_bif = 12:14;

clear meanthetapow
meanthetapow = zeros(13,3,2);
for i = 1:13
    meanthetapow(i,:,1) = nanmean(controlData3(st(i,1):st(i,2),before_bif,1));
    meanthetapow(i,:,2) = nanmean(controlData3(st(i,1):st(i,2),after_bif,1));
end

avg_theta_pow = squeeze(mean(meanthetapow,2));

clear meanthetacoh
meanthetacoh = zeros(13,3,2);
for i = 1:13
    meanthetacoh(i,:,1) = nanmean(controldataCoh(st(i,1):st(i,2),before_bif,1));
    meanthetacoh(i,:,2) = nanmean(controldataCoh(st(i,1):st(i,2),after_bif,1));
end

avg_theta_coh = squeeze(mean(meanthetacoh,2));

clear meanthetapcoh
meanthetapcoh = zeros(13,3,2);
for i = 1:13
    meanthetapcoh(i,:,1) = nanmean(peakCohData(st(i,1):st(i,2),before_bif));
    meanthetapcoh(i,:,2) = nanmean(peakCohData(st(i,1):st(i,2),after_bif));
end

avg_theta_pcoh = squeeze(mean(meanthetapcoh,2));

clear meanthetappow
meanthetappow = zeros(13,3,2);
for i = 1:13
    meanthetappow(i,:,1) = nanmean(peakPowData(st(i,1):st(i,2),before_bif));
    meanthetappow(i,:,2) = nanmean(peakPowData(st(i,1):st(i,2),after_bif));
end

avg_theta_ppow = squeeze(mean(meanthetappow,2));
%% D 

clear all
load('D:\Github\Coherogram_PFC_CA1\session_datapaths.mat');
paths = values(session_datapaths);
ids = keys(session_datapaths);
load('D:\Dataset_PFC2\Project data\TDRatio.mat');

chPFC = 1;
p=0.5;

for event = 1:4
    lfpAllCA1=[];
    lfpAllPFC=[];
    for sess=1:13

        
        [srate_lfp, dt, lfpPFC, lfpCA1, srate_bhv, whlrld, whlrl_speed,...
            SessionNP] = LoadSessionData(paths{sess});
        
        ids{sess}
        
            if(sess>10)
                chPFC=11;
            end
        
        n_trials = length(unique(whlrld(:,6)))-1;
        
        for trial=1:n_trials
            start_time = SessionNP(trial,1)*1250;
            end_time = SessionNP(trial,3)*1250;
            trial_duration = (end_time - start_time)/1250;
            
                      
            ts = SessionNP(trial,2); % last nose-poking
            ts_idx = find(whlrld(:,6)==trial,1,'first');
            if(whlrld(ts_idx,5)==2)
                cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.30,1,'first')...
                    /srate_bhv);
            elseif (whlrld(ts_idx,5)==1)
                cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.24,1,'first')...
                    /srate_bhv);
            end
%             cp2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.38,1,'first')...
%                 /srate_bhv);
            turn = (find(whlrld(:,6)==trial & whlrld(:,7)>1.46,1,'first')...
                /srate_bhv);
            turn2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.58,1,'first')...
                /srate_bhv);
            re = SessionNP(trial,3);
            
%             cpM = (cp+cp2)/2;
            tuM = (turn+turn2)/2;
            
            if(isempty(re))
                continue
            end
            
            if event == 1
                epoch = round((ts-p)*srate_lfp):round((ts+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 2
                epoch = round((cpM-p)*srate_lfp):round((cpM+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 3
                epoch = round((tuM-p)*srate_lfp):round((tuM+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 4
                epoch = round((re-p)*srate_lfp):round((re+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            
            
        end
    end
    
    step = 125;
    count=0;
    subepoch =[1:(srate_lfp/2)]+count*step;
    
    clear CxyAll
    
    while subepoch(end)<=length(lfpAllCA1)
        temp1 = reshape(lfpAllCA1(:,subepoch)',[],1);
        temp2 = reshape(lfpAllPFC(:,subepoch)',[],1);
        
        [Cxy F] = mscohere(temp1,temp2,0.5*srate_lfp,0,2^13,srate_lfp);
        CxyAll(1+count,:) = Cxy;
        
        Tcoh(1+count) = mean(subepoch)*dt;
        count = count+1;
        subepoch = [1:(srate_lfp/2)]+count*step;
        
    end
    CohEvent(event,:,:) = CxyAll;
end

%% E

figure()
plot(mean(ThetaCoh))
hold on
errorbar(mean(ThetaCoh),std(ThetaCoh)/sqrt(count))
hold off
xlim([0 5])

=======
%%% Scritps for generating the data for Figure 2 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

k = 7;
plotAux = [];
clear linPos_Coh
linPos_Coh(30).coh =[];
linPosCoh = zeros(13,30,5001);

for sess=1:13
    
    [srate_lfp, dt, lfpPFC, lfpCA1, srate_bhv, whlrld, whlrl_speed,...
        SessionNP] = LoadSessionData(paths{sess});
    
    plotAux = [];
    clear linPos_Coh
    linPos_Coh(30).coh =[];
    
    if(sess>10)
        chPFC =11;
    end
    
    for trial=1:trials(sess)
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
%         if(whlrld(ts_idx,5)~=2) % 1 = R, 2 = L
%             continue
%         end

        totTrial = sum(trials(1:sess-1)) + trial
        [ids{sess} ' ' num2str(trial)]
        start_time = (SessionNP(trial,1)-k)*srate_lfp;
        end_time = (SessionNP(trial,3)+k)*srate_lfp;
        ts = SessionNP(trial,2) - SessionNP(trial,1)+k;
        te = SessionNP(trial,3) - SessionNP(trial,1)+k;
      
        
        lfpCA1t=lfpCA1(:,round(start_time):round(end_time));
        lfpPFCt=lfpPFC(:,round(start_time):round(end_time));
        
       
        windowlength = 4*srate_lfp;
        stepsize = round(srate_lfp/8);
        
        Nwindow = (length(lfpCA1t)-windowlength)/stepsize+1;
        
        clear Cxy_T T
        for nwin = 1:Nwindow
            win = (1:windowlength) + (nwin-1)* stepsize;
            [Cxy, F]=mscohere(lfpCA1t(chCA1(sess),win),lfpPFCt(chPFC,win),1*srate_lfp,...
                [],10000,srate_lfp);
            Cxy_T(nwin,:)=Cxy;
            T(nwin) = median(win)*dt;
        end
       
        p  = round(SessionNP(trial,1)*srate_bhv);
        p2  = round(SessionNP(trial,3)*srate_bhv);
        aux = whlrld(round(p:4.875:p2),7);
        tmp = 1:0.05:2;
        
        Cxy_trial = Cxy_T(find(T<k,1,'last')+1:find(T>te,1,'first'),:);
        if size(Cxy_trial,1)<size(aux,1)
           aux = removerows(aux,1:size(aux,1)-size(Cxy_trial,1));
        end
        
        for posi = 1:5
            ind = [ts-6+posi ts-5+posi];
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_T(T>ind(1)&T<ind(2),:));
        end
        
        for posi = 6:25
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_trial(aux > tmp(posi-5) & aux < tmp(posi-4),:));              
        end
        
        for posi = 26:30
            ind = [te+posi-26 te+posi-25];
            linPos_Coh(posi).coh = vertcat(linPos_Coh(posi).coh,...
                Cxy_T(T>ind(1)&T<ind(2),:));
        end
        
    end
    for j=1:30
        if(size(linPos_Coh(j).coh,1)>1)
            linPos_Coh(j).coh = mean(linPos_Coh(j).coh);
        end
        plotAux = vertcat(plotAux, linPos_Coh(j).coh);
    end
    linPosCoh(sess,:,:) = plotAux;
end


plotAux = [];

for j=1:30
    if(size(linPos_Coh(j).coh,1)>1)
        linPos_Coh(j).coh = mean(linPos_Coh(j).coh);
    end
    plotAux = vertcat(plotAux, linPos_Coh(j).coh);
end

%% C

before_bif = 4:6;
after_bif = 12:14;

clear meanthetapow
meanthetapow = zeros(13,3,2);
for i = 1:13
    meanthetapow(i,:,1) = nanmean(controlData3(st(i,1):st(i,2),before_bif,1));
    meanthetapow(i,:,2) = nanmean(controlData3(st(i,1):st(i,2),after_bif,1));
end

avg_theta_pow = squeeze(mean(meanthetapow,2));

clear meanthetacoh
meanthetacoh = zeros(13,3,2);
for i = 1:13
    meanthetacoh(i,:,1) = nanmean(controldataCoh(st(i,1):st(i,2),before_bif,1));
    meanthetacoh(i,:,2) = nanmean(controldataCoh(st(i,1):st(i,2),after_bif,1));
end

avg_theta_coh = squeeze(mean(meanthetacoh,2));

clear meanthetapcoh
meanthetapcoh = zeros(13,3,2);
for i = 1:13
    meanthetapcoh(i,:,1) = nanmean(peakCohData(st(i,1):st(i,2),before_bif));
    meanthetapcoh(i,:,2) = nanmean(peakCohData(st(i,1):st(i,2),after_bif));
end

avg_theta_pcoh = squeeze(mean(meanthetapcoh,2));

clear meanthetappow
meanthetappow = zeros(13,3,2);
for i = 1:13
    meanthetappow(i,:,1) = nanmean(peakPowData(st(i,1):st(i,2),before_bif));
    meanthetappow(i,:,2) = nanmean(peakPowData(st(i,1):st(i,2),after_bif));
end

avg_theta_ppow = squeeze(mean(meanthetappow,2));
%% D 

clear all
load('D:\Github\Coherogram_PFC_CA1\session_datapaths.mat');
paths = values(session_datapaths);
ids = keys(session_datapaths);
load('D:\Dataset_PFC2\Project data\TDRatio.mat');

chPFC = 1;
p=0.5;

for event = 1:4
    lfpAllCA1=[];
    lfpAllPFC=[];
    for sess=1:13

        
        [srate_lfp, dt, lfpPFC, lfpCA1, srate_bhv, whlrld, whlrl_speed,...
            SessionNP] = LoadSessionData(paths{sess});
        
        ids{sess}
        
            if(sess>10)
                chPFC=11;
            end
        
        n_trials = length(unique(whlrld(:,6)))-1;
        
        for trial=1:n_trials
            start_time = SessionNP(trial,1)*1250;
            end_time = SessionNP(trial,3)*1250;
            trial_duration = (end_time - start_time)/1250;
            
                      
            ts = SessionNP(trial,2); % last nose-poking
            ts_idx = find(whlrld(:,6)==trial,1,'first');
            if(whlrld(ts_idx,5)==2)
                cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.30,1,'first')...
                    /srate_bhv);
            elseif (whlrld(ts_idx,5)==1)
                cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.24,1,'first')...
                    /srate_bhv);
            end
%             cp2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.38,1,'first')...
%                 /srate_bhv);
            turn = (find(whlrld(:,6)==trial & whlrld(:,7)>1.46,1,'first')...
                /srate_bhv);
            turn2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.58,1,'first')...
                /srate_bhv);
            re = SessionNP(trial,3);
            
%             cpM = (cp+cp2)/2;
            tuM = (turn+turn2)/2;
            
            if(isempty(re))
                continue
            end
            
            if event == 1
                epoch = round((ts-p)*srate_lfp):round((ts+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 2
                epoch = round((cpM-p)*srate_lfp):round((cpM+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 3
                epoch = round((tuM-p)*srate_lfp):round((tuM+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            if event == 4
                epoch = round((re-p)*srate_lfp):round((re+p)*srate_lfp);
                lfpAllCA1 = [lfpAllCA1;lfpCA1(chCA1(sess),epoch)];
                lfpAllPFC = [lfpAllPFC;lfpPFC(chPFC,epoch)];
            end
            
            
        end
    end
    
    step = 125;
    count=0;
    subepoch =[1:(srate_lfp/2)]+count*step;
    
    clear CxyAll
    
    while subepoch(end)<=length(lfpAllCA1)
        temp1 = reshape(lfpAllCA1(:,subepoch)',[],1);
        temp2 = reshape(lfpAllPFC(:,subepoch)',[],1);
        
        [Cxy F] = mscohere(temp1,temp2,0.5*srate_lfp,0,2^13,srate_lfp);
        CxyAll(1+count,:) = Cxy;
        
        Tcoh(1+count) = mean(subepoch)*dt;
        count = count+1;
        subepoch = [1:(srate_lfp/2)]+count*step;
        
    end
    CohEvent(event,:,:) = CxyAll;
end

%% E

figure()
plot(mean(ThetaCoh))
hold on
errorbar(mean(ThetaCoh),std(ThetaCoh)/sqrt(count))
hold off
xlim([0 5])

>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
