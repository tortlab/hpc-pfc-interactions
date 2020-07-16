<<<<<<< HEAD
%%% Scritps for generating the data for Figure 3 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

sessSpec = [];
sess_ts_GC = [];
sess_cp_GC = [];
sess_tu_GC = [];
sess_re_GC = [];

for sess=1:13
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    p = round(0.25*srate_lfp);
    ts_lfp = [];
    cp_lfp = [];
    tu_lfp = [];
    re_lfp = [];
    for trial = 1:n_trials
        disp([sess trial])
        ts = SessionNP(trial,2); % last nose-poking
        re = SessionNP(trial,3);
        
        if(isempty(re))
            continue
        end
        
        % Getting the epoch of the event inside the trial
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
        if(whlrld(ts_idx,5)==2)
            cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.30 ...
                & whlrld(:,7)<1.50,1,'first')/srate_bhv);
        elseif (whlrld(ts_idx,5)==1)
            cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.24 ...
                & whlrld(:,7)<1.50,1,'first')/srate_bhv);
        end
        
        turn = (find(whlrld(:,6)==trial & whlrld(:,7)>1.46...
            & whlrld(:,7)<1.7,1,'first')/srate_bhv);
        turn2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.58...
            & whlrld(:,7)<1.7,1,'first')/srate_bhv);
        
        tuM = (turn+turn2)/2;
        for event = 1:4
            if event == 1
                trial_epoch = round(ts*srate_lfp-p):round(ts*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                ts_lfp = [ts_lfp, ds_lfp];
      
            elseif event == 2
                trial_epoch = round(cpM*srate_lfp-p):round(cpM*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                cp_lfp = [cp_lfp, ds_lfp];

            elseif event == 3
                trial_epoch = round(tuM*srate_lfp-p):round(tuM*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                tu_lfp = [tu_lfp, ds_lfp];
            else
                
                trial_epoch = round(re*srate_lfp-p):round(re*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                re_lfp = [re_lfp, ds_lfp];
            end
        end
    end
    
    X1 = zscore(ts_lfp')';
    X2 = zscore(cp_lfp')';
    X3 = zscore(tu_lfp')';
    X4 = zscore(re_lfp')';
    morder = 15;
    
    [A2,Sig,E2]= tsdata_to_var([X1],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    
    ts_GC(1,:) = squeeze(F_spect(1,2,:));
    ts_GC(2,:) = squeeze(F_spect(2,1,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [A2,Sig,E2]= tsdata_to_var([X2],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    cp_GC(1,:) = squeeze(F_spect(1,2,:));
    cp_GC(2,:) = squeeze(F_spect(2,1,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [A2,Sig,E2]= tsdata_to_var([X3],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    
    tu_GC(1,:) = squeeze(F_spect(1,2,:));
    tu_GC(2,:) = squeeze(F_spect(2,1,:));
    
    [A2,Sig,E2]= tsdata_to_var([X4],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    re_GC(1,:) = squeeze(F_spect(1,2,:));
    re_GC(2,:) = squeeze(F_spect(2,1,:));
    
    sess_ts_GC(sess,:,:) = ts_GC;
    sess_cp_GC(sess,:,:) = cp_GC;
    sess_tu_GC(sess,:,:) = tu_GC;
    sess_re_GC(sess,:,:) = re_GC;
    
end

%% B

tmp = 1:0.1:2; % The spatial bin boundaries
nbins = length(tmp)-1;
bin_lfp = cell(13,10); % Initializing the cell which will contain all segments

for sess=1:13
    
    % Loading the data for this session
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    for trial = 1:n_trials
        
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
        if(whlrld(ts_idx,5)~=1) % 1 = Right, 2 = Left
            continue
        end
        
        disp([sess trial])
        ts = SessionNP(trial,2); % Last nose-poking
        re = SessionNP(trial,3); % End of trial
        % Getting the trial LFP segment
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        trial_lfpCA1 = lfpCA1(1,trial_epoch);
        trial_lfpPFC = lfpPFC(1,trial_epoch);
        trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
        % Trial timing with the behavior sampling rate
        trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
        
        for bin = 1:nbins
            % Defining the segment of the position data
            bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                whlrld(trial_bhv,7)<tmp(bin+1)));
            % Converting the segment time from behavior to LFP sampling rate
            bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
            bin_epoch = bin_epoch(bin_epoch>0); % Negatives can happen in the first value
            % Sometimes the position value jumps a bin (fast movements)
            if(isempty(bin_epoch))
                continue
            end
            % Extracting the bin-associated segment in the signal and concatenating it
            temp_lfp = trial_lfp(:,bin_epoch(1):bin_epoch(end));
            bin_lfp{sess,bin} = [bin_lfp{sess,bin}, temp_lfp];
        end
        
    end
end

ds_bin_lfp = cell(13,10);
for sess = 1:13
    for bin = 1:10
        ds_bin_lfp{sess,bin}(1,:) = downsample(bin_lfp{sess,bin}(1,:),10);
        ds_bin_lfp{sess,bin}(2,:) = downsample(bin_lfp{sess,bin}(2,:),10);
    end
end

morder = 15;
srate=125;
freqs = sfreqs(1000,srate);
F_spect = nan(13,10,2,2,length(freqs));

for sess = 1:13
    for bin = 1:10
        disp([sess bin])
        clear X
        X = ds_bin_lfp{sess,bin};
        [A2,Sig,E2]= tsdata_to_var(X,morder);
        [G,info] = var_to_autocov(A2,Sig);
        [F_spect(sess,bin,:,:,:),fres] = autocov_to_spwcgc(G,1000,[]);
    end
end

F_spectR = F_spect;
=======
%%% Scritps for generating the data for Figure 3 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

sessSpec = [];
sess_ts_GC = [];
sess_cp_GC = [];
sess_tu_GC = [];
sess_re_GC = [];

for sess=1:13
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    p = round(0.25*srate_lfp);
    ts_lfp = [];
    cp_lfp = [];
    tu_lfp = [];
    re_lfp = [];
    for trial = 1:n_trials
        disp([sess trial])
        ts = SessionNP(trial,2); % last nose-poking
        re = SessionNP(trial,3);
        
        if(isempty(re))
            continue
        end
        
        % Getting the epoch of the event inside the trial
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
        if(whlrld(ts_idx,5)==2)
            cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.30 ...
                & whlrld(:,7)<1.50,1,'first')/srate_bhv);
        elseif (whlrld(ts_idx,5)==1)
            cpM = (find(whlrld(:,6)==trial & whlrld(:,7)>1.24 ...
                & whlrld(:,7)<1.50,1,'first')/srate_bhv);
        end
        
        turn = (find(whlrld(:,6)==trial & whlrld(:,7)>1.46...
            & whlrld(:,7)<1.7,1,'first')/srate_bhv);
        turn2 = (find(whlrld(:,6)==trial & whlrld(:,7)>1.58...
            & whlrld(:,7)<1.7,1,'first')/srate_bhv);
        
        tuM = (turn+turn2)/2;
        for event = 1:4
            if event == 1
                trial_epoch = round(ts*srate_lfp-p):round(ts*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                ts_lfp = [ts_lfp, ds_lfp];
      
            elseif event == 2
                trial_epoch = round(cpM*srate_lfp-p):round(cpM*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                cp_lfp = [cp_lfp, ds_lfp];

            elseif event == 3
                trial_epoch = round(tuM*srate_lfp-p):round(tuM*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                tu_lfp = [tu_lfp, ds_lfp];
            else
                
                trial_epoch = round(re*srate_lfp-p):round(re*srate_lfp+p);
                
                trial_lfpCA1 = lfpCA1(1,trial_epoch);
                trial_lfpPFC = lfpPFC(1,trial_epoch);
                trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
                
                clear ds_lfp
                ds_lfp(1,:) = downsample(trial_lfp(1,:),10);
                ds_lfp(2,:) = downsample(trial_lfp(2,:),10);
                
                re_lfp = [re_lfp, ds_lfp];
            end
        end
    end
    
    X1 = zscore(ts_lfp')';
    X2 = zscore(cp_lfp')';
    X3 = zscore(tu_lfp')';
    X4 = zscore(re_lfp')';
    morder = 15;
    
    [A2,Sig,E2]= tsdata_to_var([X1],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    
    ts_GC(1,:) = squeeze(F_spect(1,2,:));
    ts_GC(2,:) = squeeze(F_spect(2,1,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [A2,Sig,E2]= tsdata_to_var([X2],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    cp_GC(1,:) = squeeze(F_spect(1,2,:));
    cp_GC(2,:) = squeeze(F_spect(2,1,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [A2,Sig,E2]= tsdata_to_var([X3],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    
    tu_GC(1,:) = squeeze(F_spect(1,2,:));
    tu_GC(2,:) = squeeze(F_spect(2,1,:));
    
    [A2,Sig,E2]= tsdata_to_var([X4],morder);
    
    [G,info] = var_to_autocov(A2,Sig);
    
    srate=125;
    freqs = sfreqs(1000,srate);
    
    [F_spect,fres] = autocov_to_spwcgc(G,1000,[]);
    re_GC(1,:) = squeeze(F_spect(1,2,:));
    re_GC(2,:) = squeeze(F_spect(2,1,:));
    
    sess_ts_GC(sess,:,:) = ts_GC;
    sess_cp_GC(sess,:,:) = cp_GC;
    sess_tu_GC(sess,:,:) = tu_GC;
    sess_re_GC(sess,:,:) = re_GC;
    
end

%% B

tmp = 1:0.1:2; % The spatial bin boundaries
nbins = length(tmp)-1;
bin_lfp = cell(13,10); % Initializing the cell which will contain all segments

for sess=1:13
    
    % Loading the data for this session
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    for trial = 1:n_trials
        
        ts_idx = find(whlrld(:,6)==trial,1,'first');
        
        if(whlrld(ts_idx,5)~=1) % 1 = Right, 2 = Left
            continue
        end
        
        disp([sess trial])
        ts = SessionNP(trial,2); % Last nose-poking
        re = SessionNP(trial,3); % End of trial
        % Getting the trial LFP segment
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        trial_lfpCA1 = lfpCA1(1,trial_epoch);
        trial_lfpPFC = lfpPFC(1,trial_epoch);
        trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
        % Trial timing with the behavior sampling rate
        trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
        
        for bin = 1:nbins
            % Defining the segment of the position data
            bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                whlrld(trial_bhv,7)<tmp(bin+1)));
            % Converting the segment time from behavior to LFP sampling rate
            bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
            bin_epoch = bin_epoch(bin_epoch>0); % Negatives can happen in the first value
            % Sometimes the position value jumps a bin (fast movements)
            if(isempty(bin_epoch))
                continue
            end
            % Extracting the bin-associated segment in the signal and concatenating it
            temp_lfp = trial_lfp(:,bin_epoch(1):bin_epoch(end));
            bin_lfp{sess,bin} = [bin_lfp{sess,bin}, temp_lfp];
        end
        
    end
end

ds_bin_lfp = cell(13,10);
for sess = 1:13
    for bin = 1:10
        ds_bin_lfp{sess,bin}(1,:) = downsample(bin_lfp{sess,bin}(1,:),10);
        ds_bin_lfp{sess,bin}(2,:) = downsample(bin_lfp{sess,bin}(2,:),10);
    end
end

morder = 15;
srate=125;
freqs = sfreqs(1000,srate);
F_spect = nan(13,10,2,2,length(freqs));

for sess = 1:13
    for bin = 1:10
        disp([sess bin])
        clear X
        X = ds_bin_lfp{sess,bin};
        [A2,Sig,E2]= tsdata_to_var(X,morder);
        [G,info] = var_to_autocov(A2,Sig);
        [F_spect(sess,bin,:,:,:),fres] = autocov_to_spwcgc(G,1000,[]);
    end
end

F_spectR = F_spect;
>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
disp('end')