<<<<<<< HEAD
%%% Scritps for generating the data for Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%%

srate = 1250; % sampling rate
slow_vector = 0:1:12; % vector of frequencies to filter the slow waves
fast_vector = 20:5:100; % vector of frequencies to filter the fast waves
slow_BandWidth = 2; % filter bandwidth (slow waves)
fast_BandWidth = 20; % filter bandwidth (fast waves)
numbin = 18; % number of phase bins (slow waves)

for sess = 1:13
    
    disp(sess)
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    SlowPhase = [];
    SlowPhase2 = [];
    FastPhase1 = [];
    FastPhase2 = [];
    clear lfp_Fast1 lfp_Fast2 lfp_Slow lfp_Slow2
    

    for ii=1:length(fast_vector)
        Ff1 = fast_vector(ii); % selecting frequency (low cut)
        Ff2=Ff1+fast_BandWidth; % selecting frequency (high cut)
        
        lfp_Fast1(:,ii) = eegfilt(lfpCA1,srate,Ff1,Ff2);
        lfp_Fast2(:,ii) = eegfilt(lfpPFC,srate,Ff1,Ff2);
    end
    
    for jj=1:length(slow_vector)
        Sf1 = slow_vector(jj); % selecting frequency (low cut)
        Sf2 = Sf1 + slow_BandWidth; % selecting frequency (high cut)
        
        lfp_Slow(:,jj) = eegfilt(lfpCA1,srate,Sf1,Sf2);
    end
    
    for jj=1:length(slow_vector)
        Sf1 = slow_vector(jj); % selecting frequency (low cut)
        Sf2 = Sf1 + slow_BandWidth; % selecting frequency (high cut)
        
        lfp_Slow2(:,jj) = eegfilt(lfpPFC,srate,Sf1,Sf2);
    end

    for trial = 1:n_trials
        disp([num2str(sess) ' ' num2str(trial)])
        ts = SessionNP(trial,1); % last nose-poking
        re = SessionNP(trial,3);
        
        % Getting the epoch of the event inside the trial
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        
        clear FastF1 FastF2 FP1 FP2 SP SP2 SlowF
        % Obtaining the fast frequency phase time-series
        for ii=1:length(fast_vector)
            FP1(ii, :) = angle(hilbert(lfp_Fast1(trial_epoch,ii))); % getting the intantaneous phase of lfp1
            FP2(ii, :) = angle(hilbert(lfp_Fast2(trial_epoch,ii))); % getting the intantaneous phase of lfp2
        end
        
        FastPhase1 = [FastPhase1,FP1];
        FastPhase2 = [FastPhase2,FP2];
        
        % Obtaining the slow frequency phase time-series
        for jj=1:length(slow_vector)
            SP(jj, :) = angle(hilbert(lfp_Slow(trial_epoch,jj)));% filtering lfp with the slow wave reference
        end
        
        for jj=1:length(slow_vector)
            SP2(jj, :) = angle(hilbert(lfp_Slow2(trial_epoch,jj)));% filtering lfp with the slow wave reference
        end
        
        SlowPhase = [SlowPhase, SP];
        SlowPhase2 = [SlowPhase2, SP2];
    end
    
    % Loop through the frequencies and compute the plv_modindex comodulogram
    plv_modindex_comodulogram = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase,1))); % pre-allocating
    for i = 1:size(SlowPhase,1) % loop through slow frequencies
        for j = 1:size(FastPhase1,1) % loop through fast frequencies
            plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
                FastPhase2(j,:)',SlowPhase(i,:)',numbin); % plv_modindex calculation
            plv_modindex_comodulogram(j,i) = plv_phase_modindex; % storing the results in the variable plv_modindex_comodulogram
        end
    end
    
    plv_modindex_comodulogram2 = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase2,1))); % pre-allocating
    for i = 1:size(SlowPhase2,1) % loop through slow frequencies
        for j = 1:size(FastPhase1,1) % loop through fast frequencies
            plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
                FastPhase2(j,:)',SlowPhase2(i,:)',numbin); % plv_modindex calculation
            plv_modindex_comodulogram2(j,i) = plv_phase_modindex; % storing the results in the variable plv_modindex_comodulogram
        end
    end
    
    plv_MI_ca1(sess,:,:) = plv_modindex_comodulogram;
    plv_MI_pfc(sess,:,:) = plv_modindex_comodulogram2;
    
    %%%%%%% Surrogate analysis %%%%%%%%
%     dist = squeeze(mean(gather(plv_MI_ca1)));
%     dist2 =squeeze(mean(gather(plv_MI_pfc)));

%     for shuffle = 1:100
%         disp(shuffle)
%         shift = randi(10,1);
%         NullPhase = circshift(SlowPhase, shift*srate_lfp, 2);
%         shift2 = randi(10,1);
%         NullPhase2 = circshift(SlowPhase2, shift2*srate_lfp, 2);
%         
%         plv_modindex_comodulogram = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase,1)));
%         for i = 1:size(SlowPhase,1)
%             for j = 1:size(FastPhase1,1)
%                 if(isnan(mat(j,i)))
%                     continue
%                 end
%                 plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
%                     FastPhase2(j,:)',NullPhase(i,:)',numbin);
%                 plv_modindex_comodulogram(j,i) = plv_phase_modindex;
%             end
%         end
%         
%         plv_modindex_comodulogram2 = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase2,1)));
%         for i = 1:size(SlowPhase2,1)
%             for j = 1:size(FastPhase1,1)
%                 if(isnan(mat2(j,i)))
%                     continue
%                 end
%                 plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
%                     FastPhase2(j,:)',NullPhase2(i,:)',numbin);
%                 plv_modindex_comodulogram2(j,i) = plv_phase_modindex;
%             end
%         end
%         
%         nulldist(lst+shuffle,:,:) = plv_modindex_comodulogram;
%         nulldist2(lst+shuffle,:,:) = plv_modindex_comodulogram2;
%     end
%     lst = lst+shuffle;
%     
%     ndist = gather(nulldist);
%     ndist2 = gather(nulldist2);
%     
%     ndist(ndist==0) = nan;
%     ndist2(ndist2==0) = nan;
%     
%     clear test test2
%     for i = 1:17
%         for j = 1:13
%             test(i,j) = gather(length(find(dist(i,j)>ndist(:,i,j))));
%         end
%     end
%     
%     for i = 1:17
%         for j = 1:13
%             test2(i,j) = gather(length(find(dist2(i,j)>ndist2(:,i,j))));
%         end
%     end
%     
%     threshold = lst-65;
%     mat = test;
%     mat(mat<threshold)= nan;
%     
%     
%     mat2 = test2;
%     mat2(mat2<threshold)= nan;
end

=======
%%% Scritps for generating the data for Figure 5 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%%

srate = 1250; % sampling rate
slow_vector = 0:1:12; % vector of frequencies to filter the slow waves
fast_vector = 20:5:100; % vector of frequencies to filter the fast waves
slow_BandWidth = 2; % filter bandwidth (slow waves)
fast_BandWidth = 20; % filter bandwidth (fast waves)
numbin = 18; % number of phase bins (slow waves)

for sess = 1:13
    
    disp(sess)
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    n_trials = length(unique(whlrld(:,6)))-1;
    SlowPhase = [];
    SlowPhase2 = [];
    FastPhase1 = [];
    FastPhase2 = [];
    clear lfp_Fast1 lfp_Fast2 lfp_Slow lfp_Slow2
    

    for ii=1:length(fast_vector)
        Ff1 = fast_vector(ii); % selecting frequency (low cut)
        Ff2=Ff1+fast_BandWidth; % selecting frequency (high cut)
        
        lfp_Fast1(:,ii) = eegfilt(lfpCA1,srate,Ff1,Ff2);
        lfp_Fast2(:,ii) = eegfilt(lfpPFC,srate,Ff1,Ff2);
    end
    
    for jj=1:length(slow_vector)
        Sf1 = slow_vector(jj); % selecting frequency (low cut)
        Sf2 = Sf1 + slow_BandWidth; % selecting frequency (high cut)
        
        lfp_Slow(:,jj) = eegfilt(lfpCA1,srate,Sf1,Sf2);
    end
    
    for jj=1:length(slow_vector)
        Sf1 = slow_vector(jj); % selecting frequency (low cut)
        Sf2 = Sf1 + slow_BandWidth; % selecting frequency (high cut)
        
        lfp_Slow2(:,jj) = eegfilt(lfpPFC,srate,Sf1,Sf2);
    end

    for trial = 1:n_trials
        disp([num2str(sess) ' ' num2str(trial)])
        ts = SessionNP(trial,1); % last nose-poking
        re = SessionNP(trial,3);
        
        % Getting the epoch of the event inside the trial
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        
        clear FastF1 FastF2 FP1 FP2 SP SP2 SlowF
        % Obtaining the fast frequency phase time-series
        for ii=1:length(fast_vector)
            FP1(ii, :) = angle(hilbert(lfp_Fast1(trial_epoch,ii))); % getting the intantaneous phase of lfp1
            FP2(ii, :) = angle(hilbert(lfp_Fast2(trial_epoch,ii))); % getting the intantaneous phase of lfp2
        end
        
        FastPhase1 = [FastPhase1,FP1];
        FastPhase2 = [FastPhase2,FP2];
        
        % Obtaining the slow frequency phase time-series
        for jj=1:length(slow_vector)
            SP(jj, :) = angle(hilbert(lfp_Slow(trial_epoch,jj)));% filtering lfp with the slow wave reference
        end
        
        for jj=1:length(slow_vector)
            SP2(jj, :) = angle(hilbert(lfp_Slow2(trial_epoch,jj)));% filtering lfp with the slow wave reference
        end
        
        SlowPhase = [SlowPhase, SP];
        SlowPhase2 = [SlowPhase2, SP2];
    end
    
    % Loop through the frequencies and compute the plv_modindex comodulogram
    plv_modindex_comodulogram = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase,1))); % pre-allocating
    for i = 1:size(SlowPhase,1) % loop through slow frequencies
        for j = 1:size(FastPhase1,1) % loop through fast frequencies
            plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
                FastPhase2(j,:)',SlowPhase(i,:)',numbin); % plv_modindex calculation
            plv_modindex_comodulogram(j,i) = plv_phase_modindex; % storing the results in the variable plv_modindex_comodulogram
        end
    end
    
    plv_modindex_comodulogram2 = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase2,1))); % pre-allocating
    for i = 1:size(SlowPhase2,1) % loop through slow frequencies
        for j = 1:size(FastPhase1,1) % loop through fast frequencies
            plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
                FastPhase2(j,:)',SlowPhase2(i,:)',numbin); % plv_modindex calculation
            plv_modindex_comodulogram2(j,i) = plv_phase_modindex; % storing the results in the variable plv_modindex_comodulogram
        end
    end
    
    plv_MI_ca1(sess,:,:) = plv_modindex_comodulogram;
    plv_MI_pfc(sess,:,:) = plv_modindex_comodulogram2;
    
    %%%%%%% Surrogate analysis %%%%%%%%
%     dist = squeeze(mean(gather(plv_MI_ca1)));
%     dist2 =squeeze(mean(gather(plv_MI_pfc)));

%     for shuffle = 1:100
%         disp(shuffle)
%         shift = randi(10,1);
%         NullPhase = circshift(SlowPhase, shift*srate_lfp, 2);
%         shift2 = randi(10,1);
%         NullPhase2 = circshift(SlowPhase2, shift2*srate_lfp, 2);
%         
%         plv_modindex_comodulogram = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase,1)));
%         for i = 1:size(SlowPhase,1)
%             for j = 1:size(FastPhase1,1)
%                 if(isnan(mat(j,i)))
%                     continue
%                 end
%                 plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
%                     FastPhase2(j,:)',NullPhase(i,:)',numbin);
%                 plv_modindex_comodulogram(j,i) = plv_phase_modindex;
%             end
%         end
%         
%         plv_modindex_comodulogram2 = gpuArray(zeros(size(FastPhase1,1),size(SlowPhase2,1)));
%         for i = 1:size(SlowPhase2,1)
%             for j = 1:size(FastPhase1,1)
%                 if(isnan(mat2(j,i)))
%                     continue
%                 end
%                 plv_phase_modindex = plv_modindex(FastPhase1(j,:)',...
%                     FastPhase2(j,:)',NullPhase2(i,:)',numbin);
%                 plv_modindex_comodulogram2(j,i) = plv_phase_modindex;
%             end
%         end
%         
%         nulldist(lst+shuffle,:,:) = plv_modindex_comodulogram;
%         nulldist2(lst+shuffle,:,:) = plv_modindex_comodulogram2;
%     end
%     lst = lst+shuffle;
%     
%     ndist = gather(nulldist);
%     ndist2 = gather(nulldist2);
%     
%     ndist(ndist==0) = nan;
%     ndist2(ndist2==0) = nan;
%     
%     clear test test2
%     for i = 1:17
%         for j = 1:13
%             test(i,j) = gather(length(find(dist(i,j)>ndist(:,i,j))));
%         end
%     end
%     
%     for i = 1:17
%         for j = 1:13
%             test2(i,j) = gather(length(find(dist2(i,j)>ndist2(:,i,j))));
%         end
%     end
%     
%     threshold = lst-65;
%     mat = test;
%     mat(mat<threshold)= nan;
%     
%     
%     mat2 = test2;
%     mat2(mat2<threshold)= nan;
end

>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
