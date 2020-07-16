<<<<<<< HEAD
%%% Scritps for generating the data for Figure 4 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

% 1 = mPFC
% 2 = CA1
%
% amplitude, phase

comb = [2,2; 1,1; 1,2; 2,1];

for order = 1:1
    PhaseFreqVector = 0:0.5:20;
    AmpFreqVector   = 10:5:100;
    PhaseFreq_BandWidth = 4;
    AmpFreq_BandWidth   = 10;
    AllComod = [];
    
    for sess=1:13
        AmpFreqCat = [];
        PhaseFreqCat = [];
        lfpPFC = lfpPFC_all{sess};
        lfpCA1 = lfpCA1_all{sess};
        whlrl_speed = whlrl_speed_all{sess};
        whlrld = whlrld_all{sess};
        SessionNP = SessionNP_all{sess};
        lfp_amp = zeros(length(lfpCA1),length(AmpFreqVector));
        lfp_phase = zeros(length(lfpPFC),length(PhaseFreqVector));
        
        parfor ii=1:length(AmpFreqVector)
            Af1 = AmpFreqVector(ii);
            Af2 = Af1+AmpFreq_BandWidth;
            lfp_amp(:,ii) = eegfilt(lfpCA1,srate_lfp,Af1,Af2);
        end
        
        parfor jj=1:length(PhaseFreqVector)
            Pf1 = PhaseFreqVector(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidth;
            lfp_phase(:,jj) = eegfilt(lfpPFC,srate_lfp,Pf1,Pf2);
        end
        
        n_trials = length(unique(whlrld(:,6)))-1;
        
        for trial = 1:n_trials
            disp([num2str(sess) ' ' num2str(trial)])
            ts = SessionNP(trial,2); % last nose-poking
            re = SessionNP(trial,3);
            
            if(isempty(re))
                continue
            end
            
            % Getting the epoch of the event inside the trial
            trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);

            clear AmpFreqTransformed PhaseFreqTransformed
            
            parfor ii=1:length(AmpFreqVector)
                AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
            end
            
            AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
            
            parfor jj=1:length(PhaseFreqVector)
                PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
            end
            
            PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
        end
        
        nbin     = 18;
        position = zeros(1,nbin);
        winsize  = 2*pi/nbin;
        
        for j=1:nbin
            position(j) = -pi+(j-1)*winsize;
        end
        
        Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        
        counter1 = 0;
        for ii = 1:length(PhaseFreqVector)
            counter1 = counter1+1;
            Pf1      = PhaseFreqVector(ii);
            Pf2      = Pf1+PhaseFreq_BandWidth;
            
            counter2 = 0;
            for jj = 1:length(AmpFreqVector)
                counter2 = counter2+1;
                
                Af1 = AmpFreqVector(jj);
                Af2 = Af1+AmpFreq_BandWidth;
                [MI,MeanAmp] = ModIndex_v2(PhaseFreqCat(ii, :), AmpFreqCat(jj, :), position);
                Comodulogram(counter1,counter2) = MI;
            end
        end
        
        AllComod(sess,:,:) = Comodulogram;
        
    end
end

%% B

% 1 = mPFC
% 2 = CA1
PhaseFreqVector = 0:0.5:20;
AmpFreqVector   = 10:5:100;
PhaseFreq_BandWidth = 4;
AmpFreq_BandWidth   = 10;

% order = 1-CA1p-CA1a, 2:PFCp-PFCa, 3:CA1p-PFCa, 4:PFCp-CA1a
order = 4;

for sess=1:13
    AmpFreqCat = [];
    PhaseFreqCat = [];
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    tmp = [1,1.2,1.4,1.6,1.8,2];
    nbins = length(tmp)-1;
    AmpFreqBin = cell(nbins,length(AmpFreqVector));
    PhaseFreqBin = cell(nbins,length(PhaseFreqVector));
    n_trials = length(unique(whlrld(:,6)))-1;
    
    lfp_amp = zeros(length(lfpCA1),length(AmpFreqVector));
    lfp_phase = zeros(length(lfpPFC),length(PhaseFreqVector));
    
    parfor ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2 = Af1+AmpFreq_BandWidth;
        lfp_amp(:,ii) = eegfilt(lfpCA1,srate_lfp,Af1,Af2);
    end
    
    parfor jj=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(jj);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        lfp_phase(:,jj) = eegfilt(lfpPFC,srate_lfp,Pf1,Pf2);
    end
    
    
    for trial = 1:n_trials
        disp([sess trial])
        ts = SessionNP(trial,2); % last nose-poking
        re = SessionNP(trial,3);
        
        if(isempty(re))
            continue
        end
        
        % Getting the epoch of the event inside the trial
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        
        clear AmpFreqTransformed PhaseFreqTransformed
        
        parfor ii=1:length(AmpFreqVector)
            AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
        end
        
        AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
        
        parfor jj=1:length(PhaseFreqVector)
            PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
        end
        
        PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
        
        trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
        ovlp = 0;
        for bin = 1:5
            bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                whlrld(trial_bhv,7)<tmp(bin+1)));
            bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
            bin_epoch = bin_epoch(bin_epoch>0);
            binSpeed(sess,bin) = mean(whlrl_speed(bin_epoch_bhv,7));
            
            if(isempty(bin_epoch))
                continue
            end
            
            for ii = 1:length(AmpFreqVector)
                AmpFreqBin{bin,ii} = horzcat(AmpFreqBin{bin,ii}, ...
                    AmpFreqTransformed(ii,bin_epoch(1):bin_epoch(end)));
            end
            
            for jj = 1:length(PhaseFreqVector)
                PhaseFreqBin{bin,jj} = horzcat(PhaseFreqBin{bin,jj}, ...
                    PhaseFreqTransformed(jj,bin_epoch(1):bin_epoch(end)));
            end
        end
    end
    
    
    nbin     = 18;
    position = zeros(1,nbin);
    winsize  = 2*pi/nbin;
    
    for j=1:nbin
        position(j) = -pi+(j-1)*winsize;
    end
    
    for bin = 1:nbins
        AmpBinMat = reshape(cell2mat(AmpFreqBin(bin,:)),...
            length(AmpFreqBin{bin,1}), length(AmpFreqVector))';
        
        PhaseBinMat = reshape(cell2mat(PhaseFreqBin(bin,:)),...
            length(PhaseFreqBin{bin,1}),length(PhaseFreqVector))';
        
        Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        
        counter1 = 0;
        for ii = 1:length(PhaseFreqVector)
            counter1 = counter1+1;
            Pf1      = PhaseFreqVector(ii);
            Pf2      = Pf1+PhaseFreq_BandWidth;
            
            counter2 = 0;
            for jj = 1:length(AmpFreqVector)
                counter2 = counter2+1;
                
                Af1 = AmpFreqVector(jj);
                Af2 = Af1+AmpFreq_BandWidth;
                [MI,MeanAmp] = ModIndex_v2(PhaseBinMat(ii, :), ...
                    AmpBinMat(jj, :), position);
                Comodulogram(counter1,counter2) = MI;
            end
        end
        BinComodulogram(bin,:,:) = Comodulogram;
    end
    
    LinPosComod(sess,:,:,:) = BinComodulogram;
    
    theta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=7 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
    
    delta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=4 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=6);
    
    hTheta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=8.5 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
    
    if order == 1
        bandfreq = theta;
    elseif order == 2
        bandfreq = delta;
    elseif order == 3
        bandfreq = delta;
    elseif order == 4
        bandfreq = hTheta;
    end
    
    FreqBinComodulogram = squeeze(mean(BinComodulogram(:,bandfreq,:),2));
end

%% C

% 1 = mPFC
% 2 = CA1
%
% amplitude, phase

% phasebands = [6; 2; 0; 8.5];
% ampbands = [60,40,50; 80,80,50; 0,0,0; 65,50,60];
% phasebandw = [4,4,0,1.5];

lfps = [{'lfpPFC'}; {'lfpCA1'}];

phasebands = [6; 2; 2; 6];

ampbands = [60,45,50; 80,80,80; 80,80,80; 65,45,65];
phasebandw = [4,4,4,4];


comb = [2,2; 1,1; 1,2; 2,1];
chPFC = [33,33,17,49,1,9,57,31,41,25,49,25,49];

for order = 1:4
    
    PhaseFreqVector = phasebands(order);
    AmpFreqVector   = ampbands(order,:);
    
    PhaseFreq_BandWidth = phasebandw(order);
    AmpFreq_BandWidth   = 10;
    
    for sess=1:13
        AmpFreqCat = [];
        PhaseFreqCat = [];
        lfpPFC = lfpPFC_all{sess};
        lfpCA1 = lfpCA1_all{sess};
        whlrl_speed = whlrl_speed_all{sess};
        whlrld = whlrld_all{sess};
        SessionNP = SessionNP_all{sess};
        tmp = [1,1.2,1.4,1.6,1.8,2];
        nbins = length(tmp)-1;
        AmpFreqBin = cell(nbins,length(AmpFreqVector));
        PhaseFreqBin = cell(nbins,length(PhaseFreqVector));
        n_trials = length(unique(whlrld(:,6)))-1;
        
        lfp_amp = zeros(length(eval(lfps{comb(order,1)})),length(AmpFreqVector));
        lfp_phase = zeros(length(eval(lfps{comb(order,2)})),length(PhaseFreqVector));

        for ii=1:length(AmpFreqVector)
            Af1 = AmpFreqVector(ii);
            Af2 = Af1+AmpFreq_BandWidth;
            lfp_amp(:,ii) = eegfilt(eval(lfps{comb(order,1)}),srate_lfp,Af1,Af2);
        end
        
        for jj=1:length(PhaseFreqVector)
            Pf1 = PhaseFreqVector(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidth;
            lfp_phase(:,jj) = eegfilt(eval(lfps{comb(order,2)}),srate_lfp,Pf1,Pf2);
        end
        
        for trial = 1:n_trials
            disp([sess trial])
            ts = SessionNP(trial,2); % last nose-poking
            re = SessionNP(trial,3);
            
            if(isempty(re))
                continue
            end
            
            % Getting the epoch of the event inside the trial
            trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
            
            trial_lfpCA1 = lfpCA1(1,trial_epoch);
            trial_lfpPFC = lfpPFC(1,trial_epoch);
            trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
            
            clear AmpFreqTransformed PhaseFreqTransformed
            
            parfor ii=1:length(AmpFreqVector)
                AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
            end
            
            AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
            
            parfor jj=1:length(PhaseFreqVector)
                PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
            end
            
            PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
            
            
            trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
            ovlp = 0;
            for bin = 1:5
                bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                    whlrld(trial_bhv,7)<tmp(bin+1)));
                bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
                bin_epoch = bin_epoch(bin_epoch>0);
                binSpeed(sess,bin) = mean(whlrl_speed(bin_epoch_bhv,7));
                
                if(isempty(bin_epoch))
                    continue
                end
                
                for ii = 1:length(AmpFreqVector)
                    AmpFreqBin{bin,ii} = horzcat(AmpFreqBin{bin,ii}, ...
                        AmpFreqTransformed(ii,bin_epoch(1):bin_epoch(end)));
                end
                
                for jj = 1:length(PhaseFreqVector)
                    PhaseFreqBin{bin,jj} = horzcat(PhaseFreqBin{bin,jj}, ...
                        PhaseFreqTransformed(jj,bin_epoch(1):bin_epoch(end)));
                end
            end
        end
        
        
        nbin     = 18;
        position = zeros(1,nbin);
        winsize  = 2*pi/nbin;
        
        for j=1:nbin
            position(j) = -pi+(j-1)*winsize;
        end
        
        for bin = 1:nbins
            AmpBinMat = reshape(cell2mat(AmpFreqBin(bin,:)),...
                length(AmpFreqBin{bin,1}), length(AmpFreqVector))';
            
            PhaseBinMat = reshape(cell2mat(PhaseFreqBin(bin,:)),...
                length(PhaseFreqBin{bin,1}),length(PhaseFreqVector))';
            
            Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
            
            counter1 = 0;
            for ii = 1:length(PhaseFreqVector)
                counter1 = counter1+1;
                Pf1      = PhaseFreqVector(ii);
                Pf2      = Pf1+PhaseFreq_BandWidth;
                
                counter2 = 0;
                for jj = 1:length(AmpFreqVector)
                    counter2 = counter2+1;
                    
                    Af1 = AmpFreqVector(jj);
                    Af2 = Af1+AmpFreq_BandWidth;
                    [MI,MeanAmp] = ModIndex_v2(PhaseBinMat(ii, :), ...
                        AmpBinMat(jj, :), position);
                    Comodulogram(counter1,counter2) = MI;
                end
            end
            BinComodulogram(bin,:,:) = Comodulogram;
        end
        
        LinPosComod(sess,:,:,:) = BinComodulogram;
        
        theta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=7 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
        
        delta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=4 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=6);
        
        hTheta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=8.5 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
        
        if order == 1
            bandfreq = theta;
        elseif order == 2
            bandfreq = delta;
        elseif order == 3
            bandfreq = delta;
        elseif order == 4
            bandfreq = hTheta;
        end
        
        FreqBinComodulogram = squeeze(mean(BinComodulogram(:,bandfreq,:),2)); 
    end
=======
%%% Scritps for generating the data for Figure 4 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort

%% A

% 1 = mPFC
% 2 = CA1
%
% amplitude, phase

comb = [2,2; 1,1; 1,2; 2,1];

for order = 1:1
    PhaseFreqVector = 0:0.5:20;
    AmpFreqVector   = 10:5:100;
    PhaseFreq_BandWidth = 4;
    AmpFreq_BandWidth   = 10;
    AllComod = [];
    
    for sess=1:13
        AmpFreqCat = [];
        PhaseFreqCat = [];
        lfpPFC = lfpPFC_all{sess};
        lfpCA1 = lfpCA1_all{sess};
        whlrl_speed = whlrl_speed_all{sess};
        whlrld = whlrld_all{sess};
        SessionNP = SessionNP_all{sess};
        lfp_amp = zeros(length(lfpCA1),length(AmpFreqVector));
        lfp_phase = zeros(length(lfpPFC),length(PhaseFreqVector));
        
        parfor ii=1:length(AmpFreqVector)
            Af1 = AmpFreqVector(ii);
            Af2 = Af1+AmpFreq_BandWidth;
            lfp_amp(:,ii) = eegfilt(lfpCA1,srate_lfp,Af1,Af2);
        end
        
        parfor jj=1:length(PhaseFreqVector)
            Pf1 = PhaseFreqVector(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidth;
            lfp_phase(:,jj) = eegfilt(lfpPFC,srate_lfp,Pf1,Pf2);
        end
        
        n_trials = length(unique(whlrld(:,6)))-1;
        
        for trial = 1:n_trials
            disp([num2str(sess) ' ' num2str(trial)])
            ts = SessionNP(trial,2); % last nose-poking
            re = SessionNP(trial,3);
            
            if(isempty(re))
                continue
            end
            
            % Getting the epoch of the event inside the trial
            trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);

            clear AmpFreqTransformed PhaseFreqTransformed
            
            parfor ii=1:length(AmpFreqVector)
                AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
            end
            
            AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
            
            parfor jj=1:length(PhaseFreqVector)
                PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
            end
            
            PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
        end
        
        nbin     = 18;
        position = zeros(1,nbin);
        winsize  = 2*pi/nbin;
        
        for j=1:nbin
            position(j) = -pi+(j-1)*winsize;
        end
        
        Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        
        counter1 = 0;
        for ii = 1:length(PhaseFreqVector)
            counter1 = counter1+1;
            Pf1      = PhaseFreqVector(ii);
            Pf2      = Pf1+PhaseFreq_BandWidth;
            
            counter2 = 0;
            for jj = 1:length(AmpFreqVector)
                counter2 = counter2+1;
                
                Af1 = AmpFreqVector(jj);
                Af2 = Af1+AmpFreq_BandWidth;
                [MI,MeanAmp] = ModIndex_v2(PhaseFreqCat(ii, :), AmpFreqCat(jj, :), position);
                Comodulogram(counter1,counter2) = MI;
            end
        end
        
        AllComod(sess,:,:) = Comodulogram;
        
    end
end

%% B

% 1 = mPFC
% 2 = CA1
PhaseFreqVector = 0:0.5:20;
AmpFreqVector   = 10:5:100;
PhaseFreq_BandWidth = 4;
AmpFreq_BandWidth   = 10;

% order = 1-CA1p-CA1a, 2:PFCp-PFCa, 3:CA1p-PFCa, 4:PFCp-CA1a
order = 4;

for sess=1:13
    AmpFreqCat = [];
    PhaseFreqCat = [];
    lfpPFC = lfpPFC_all{sess};
    lfpCA1 = lfpCA1_all{sess};
    whlrl_speed = whlrl_speed_all{sess};
    whlrld = whlrld_all{sess};
    SessionNP = SessionNP_all{sess};
    tmp = [1,1.2,1.4,1.6,1.8,2];
    nbins = length(tmp)-1;
    AmpFreqBin = cell(nbins,length(AmpFreqVector));
    PhaseFreqBin = cell(nbins,length(PhaseFreqVector));
    n_trials = length(unique(whlrld(:,6)))-1;
    
    lfp_amp = zeros(length(lfpCA1),length(AmpFreqVector));
    lfp_phase = zeros(length(lfpPFC),length(PhaseFreqVector));
    
    parfor ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2 = Af1+AmpFreq_BandWidth;
        lfp_amp(:,ii) = eegfilt(lfpCA1,srate_lfp,Af1,Af2);
    end
    
    parfor jj=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(jj);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        lfp_phase(:,jj) = eegfilt(lfpPFC,srate_lfp,Pf1,Pf2);
    end
    
    
    for trial = 1:n_trials
        disp([sess trial])
        ts = SessionNP(trial,2); % last nose-poking
        re = SessionNP(trial,3);
        
        if(isempty(re))
            continue
        end
        
        % Getting the epoch of the event inside the trial
        trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
        
        clear AmpFreqTransformed PhaseFreqTransformed
        
        parfor ii=1:length(AmpFreqVector)
            AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
        end
        
        AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
        
        parfor jj=1:length(PhaseFreqVector)
            PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
        end
        
        PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
        
        trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
        ovlp = 0;
        for bin = 1:5
            bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                whlrld(trial_bhv,7)<tmp(bin+1)));
            bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
            bin_epoch = bin_epoch(bin_epoch>0);
            binSpeed(sess,bin) = mean(whlrl_speed(bin_epoch_bhv,7));
            
            if(isempty(bin_epoch))
                continue
            end
            
            for ii = 1:length(AmpFreqVector)
                AmpFreqBin{bin,ii} = horzcat(AmpFreqBin{bin,ii}, ...
                    AmpFreqTransformed(ii,bin_epoch(1):bin_epoch(end)));
            end
            
            for jj = 1:length(PhaseFreqVector)
                PhaseFreqBin{bin,jj} = horzcat(PhaseFreqBin{bin,jj}, ...
                    PhaseFreqTransformed(jj,bin_epoch(1):bin_epoch(end)));
            end
        end
    end
    
    
    nbin     = 18;
    position = zeros(1,nbin);
    winsize  = 2*pi/nbin;
    
    for j=1:nbin
        position(j) = -pi+(j-1)*winsize;
    end
    
    for bin = 1:nbins
        AmpBinMat = reshape(cell2mat(AmpFreqBin(bin,:)),...
            length(AmpFreqBin{bin,1}), length(AmpFreqVector))';
        
        PhaseBinMat = reshape(cell2mat(PhaseFreqBin(bin,:)),...
            length(PhaseFreqBin{bin,1}),length(PhaseFreqVector))';
        
        Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        
        counter1 = 0;
        for ii = 1:length(PhaseFreqVector)
            counter1 = counter1+1;
            Pf1      = PhaseFreqVector(ii);
            Pf2      = Pf1+PhaseFreq_BandWidth;
            
            counter2 = 0;
            for jj = 1:length(AmpFreqVector)
                counter2 = counter2+1;
                
                Af1 = AmpFreqVector(jj);
                Af2 = Af1+AmpFreq_BandWidth;
                [MI,MeanAmp] = ModIndex_v2(PhaseBinMat(ii, :), ...
                    AmpBinMat(jj, :), position);
                Comodulogram(counter1,counter2) = MI;
            end
        end
        BinComodulogram(bin,:,:) = Comodulogram;
    end
    
    LinPosComod(sess,:,:,:) = BinComodulogram;
    
    theta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=7 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
    
    delta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=4 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=6);
    
    hTheta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=8.5 ...
        & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
    
    if order == 1
        bandfreq = theta;
    elseif order == 2
        bandfreq = delta;
    elseif order == 3
        bandfreq = delta;
    elseif order == 4
        bandfreq = hTheta;
    end
    
    FreqBinComodulogram = squeeze(mean(BinComodulogram(:,bandfreq,:),2));
end

%% C

% 1 = mPFC
% 2 = CA1
%
% amplitude, phase

% phasebands = [6; 2; 0; 8.5];
% ampbands = [60,40,50; 80,80,50; 0,0,0; 65,50,60];
% phasebandw = [4,4,0,1.5];

lfps = [{'lfpPFC'}; {'lfpCA1'}];

phasebands = [6; 2; 2; 6];

ampbands = [60,45,50; 80,80,80; 80,80,80; 65,45,65];
phasebandw = [4,4,4,4];


comb = [2,2; 1,1; 1,2; 2,1];
chPFC = [33,33,17,49,1,9,57,31,41,25,49,25,49];

for order = 1:4
    
    PhaseFreqVector = phasebands(order);
    AmpFreqVector   = ampbands(order,:);
    
    PhaseFreq_BandWidth = phasebandw(order);
    AmpFreq_BandWidth   = 10;
    
    for sess=1:13
        AmpFreqCat = [];
        PhaseFreqCat = [];
        lfpPFC = lfpPFC_all{sess};
        lfpCA1 = lfpCA1_all{sess};
        whlrl_speed = whlrl_speed_all{sess};
        whlrld = whlrld_all{sess};
        SessionNP = SessionNP_all{sess};
        tmp = [1,1.2,1.4,1.6,1.8,2];
        nbins = length(tmp)-1;
        AmpFreqBin = cell(nbins,length(AmpFreqVector));
        PhaseFreqBin = cell(nbins,length(PhaseFreqVector));
        n_trials = length(unique(whlrld(:,6)))-1;
        
        lfp_amp = zeros(length(eval(lfps{comb(order,1)})),length(AmpFreqVector));
        lfp_phase = zeros(length(eval(lfps{comb(order,2)})),length(PhaseFreqVector));

        for ii=1:length(AmpFreqVector)
            Af1 = AmpFreqVector(ii);
            Af2 = Af1+AmpFreq_BandWidth;
            lfp_amp(:,ii) = eegfilt(eval(lfps{comb(order,1)}),srate_lfp,Af1,Af2);
        end
        
        for jj=1:length(PhaseFreqVector)
            Pf1 = PhaseFreqVector(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidth;
            lfp_phase(:,jj) = eegfilt(eval(lfps{comb(order,2)}),srate_lfp,Pf1,Pf2);
        end
        
        for trial = 1:n_trials
            disp([sess trial])
            ts = SessionNP(trial,2); % last nose-poking
            re = SessionNP(trial,3);
            
            if(isempty(re))
                continue
            end
            
            % Getting the epoch of the event inside the trial
            trial_epoch = round(ts*srate_lfp):round(re*srate_lfp);
            
            trial_lfpCA1 = lfpCA1(1,trial_epoch);
            trial_lfpPFC = lfpPFC(1,trial_epoch);
            trial_lfp = [lfpPFC(1,trial_epoch); lfpCA1(1,trial_epoch)];
            
            clear AmpFreqTransformed PhaseFreqTransformed
            
            parfor ii=1:length(AmpFreqVector)
                AmpFreqTransformed(ii,:) = abs(hilbert(lfp_amp(trial_epoch,ii)));
            end
            
            AmpFreqCat = [AmpFreqCat, AmpFreqTransformed];
            
            parfor jj=1:length(PhaseFreqVector)
                PhaseFreqTransformed(jj,:) = angle(hilbert(lfp_phase(trial_epoch,jj)));
            end
            
            PhaseFreqCat = [PhaseFreqCat, PhaseFreqTransformed];
            
            
            trial_bhv = round(ts*srate_bhv):round(re*srate_bhv);
            ovlp = 0;
            for bin = 1:5
                bin_epoch_bhv = trial_bhv(find(whlrld(trial_bhv,7)>tmp(bin) & ...
                    whlrld(trial_bhv,7)<tmp(bin+1)));
                bin_epoch = round((bin_epoch_bhv/srate_bhv)*srate_lfp - ts*srate_lfp);
                bin_epoch = bin_epoch(bin_epoch>0);
                binSpeed(sess,bin) = mean(whlrl_speed(bin_epoch_bhv,7));
                
                if(isempty(bin_epoch))
                    continue
                end
                
                for ii = 1:length(AmpFreqVector)
                    AmpFreqBin{bin,ii} = horzcat(AmpFreqBin{bin,ii}, ...
                        AmpFreqTransformed(ii,bin_epoch(1):bin_epoch(end)));
                end
                
                for jj = 1:length(PhaseFreqVector)
                    PhaseFreqBin{bin,jj} = horzcat(PhaseFreqBin{bin,jj}, ...
                        PhaseFreqTransformed(jj,bin_epoch(1):bin_epoch(end)));
                end
            end
        end
        
        
        nbin     = 18;
        position = zeros(1,nbin);
        winsize  = 2*pi/nbin;
        
        for j=1:nbin
            position(j) = -pi+(j-1)*winsize;
        end
        
        for bin = 1:nbins
            AmpBinMat = reshape(cell2mat(AmpFreqBin(bin,:)),...
                length(AmpFreqBin{bin,1}), length(AmpFreqVector))';
            
            PhaseBinMat = reshape(cell2mat(PhaseFreqBin(bin,:)),...
                length(PhaseFreqBin{bin,1}),length(PhaseFreqVector))';
            
            Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
            
            counter1 = 0;
            for ii = 1:length(PhaseFreqVector)
                counter1 = counter1+1;
                Pf1      = PhaseFreqVector(ii);
                Pf2      = Pf1+PhaseFreq_BandWidth;
                
                counter2 = 0;
                for jj = 1:length(AmpFreqVector)
                    counter2 = counter2+1;
                    
                    Af1 = AmpFreqVector(jj);
                    Af2 = Af1+AmpFreq_BandWidth;
                    [MI,MeanAmp] = ModIndex_v2(PhaseBinMat(ii, :), ...
                        AmpBinMat(jj, :), position);
                    Comodulogram(counter1,counter2) = MI;
                end
            end
            BinComodulogram(bin,:,:) = Comodulogram;
        end
        
        LinPosComod(sess,:,:,:) = BinComodulogram;
        
        theta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=7 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
        
        delta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=4 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=6);
        
        hTheta = find(PhaseFreqVector+PhaseFreq_BandWidth/2>=8.5 ...
            & PhaseFreqVector+PhaseFreq_BandWidth/2<=10);
        
        if order == 1
            bandfreq = theta;
        elseif order == 2
            bandfreq = delta;
        elseif order == 3
            bandfreq = delta;
        elseif order == 4
            bandfreq = hTheta;
        end
        
        FreqBinComodulogram = squeeze(mean(BinComodulogram(:,bandfreq,:),2)); 
    end
>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
end