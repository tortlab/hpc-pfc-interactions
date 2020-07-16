%%% Scritps for generating the data for Figure 6 of Hippocampal-Prefrontal
%%% interactions during Decision-Making, https://doi.org/10.1101/2020.06.24.168732
%%% @author Lucas CS Tavares under the supervision of Adriano BL Tort
%% A

for sess = 1:13
    disp(sess)
    SessionNP = SessionNP_all{sess};
    ntrials = size(SessionNP,1);
    clear spikeind_pfc spikeind_ca1
    spikeind_pfc = spikeind_PFC{sess};
    spikeind_ca1 = spikeind_CA1{sess};
    
    lfpCA1 = lfpCA1_all{sess};
    lfpPFC = lfpPFC_all{sess};
    nCells_pfc = length(unique(spikeind_pfc));
    nCells_ca1 = length(unique(spikeind_ca1));
    phase_ca1 = zeros(50,length(lfpCA1));
    phase_pfc = zeros(50,length(lfpPFC));
    whlrld = whlrld_all{sess};
    epoch = [];
    
    % Start the phase vectors
    parfor ii = 1:50
        phase_ca1(ii,:) = angle(hilbert(eegfilt(lfpCA1,srate_lfp,ii,ii+1)));
        phase_pfc(ii,:) = angle(hilbert(eegfilt(lfpPFC,srate_lfp,ii,ii+1)));
    end
    
    % Get spike timings in seconds for PFC pyr cells
    for ncell = 1:nCells_pfc
        cell_spks_pfc{ncell} = vertcat(spkts_pfc_t{sess}{:,ncell});
    end
    
    for ncell = 1:nCells_ca1
        cell_spks_ca1{ncell} = vertcat(spkts_ca1_t{sess}{:,ncell});
    end
    
    clear pall kappall RallHH RallHP RallPH RallPP
    
    for trial = 1:ntrials
        epoch = [epoch,round(SessionNP(trial,2)*srate_lfp):...
            round(SessionNP(trial,3)*srate_lfp)];
    end
    
    % Compute PLV for each bin
    for j = 1:nCells_ca1
        spks_indices=round(cell_spks_ca1{1,j}*srate_lfp);
        
        if(isempty(spks_indices))
            continue
        end
        
        if(length(intersect(spks_indices,epoch))<10)
            continue
        end
        
        for ii = 1:50
            [meanangle,meanvectorlength,angleSD,CI,kappa]= ...
                anglemean(phase_ca1...
                (ii,intersect(spks_indices,epoch)));
            
            [meanangle2,meanvectorlength2,angleSD2,CI2,kappa2]= ...
                anglemean(phase_pfc...
                (ii,intersect(spks_indices,epoch)));
            
            RallHH(j,ii)=meanvectorlength;
            RallPH(j,ii)=meanvectorlength2;
        end
    end
    
    Rsess_all_HH{sess} = RallHH;
    Rsess_all_PH{sess} = RallPH;
    
    for j = 1:nCells_pfc
        spks_indices=round(cell_spks_pfc{1,j}*srate_lfp);
        
        if(isempty(spks_indices))
            continue
        end
        
        if(length(intersect(spks_indices,epoch))<10)
            continue
        end
        
        for ii = 1:50
            [meanangle,meanvectorlength,angleSD,CI,kappa]= ...
                anglemean(phase_ca1...
                (ii,intersect(spks_indices,epoch)));
            
            [meanangle2,meanvectorlength2,angleSD2,CI2,kappa2]= ...
                anglemean(phase_pfc...
                (ii,intersect(spks_indices,epoch)));
            
            RallHP(j,ii)=meanvectorlength;
            RallPP(j,ii)=meanvectorlength2;
        end
    end
    Rsess_all_HP{sess} = RallHP;
    Rsess_all_PP{sess} = RallPP;
end

%% B 

for sess = 4:13
    disp(sess)
    SessionNP = SessionNP_all{sess};
    ntrials = size(SessionNP,1);
    clear spikeind_pfc spikeind_ca1
    spikeind_pfc = spikeind_PFC{sess};
    spikeind_ca1 = spikeind_CA1{sess};
    
    lfpCA1 = lfpCA1_all{sess};
    lfpPFC = lfpPFC_all{sess};
    nCells_pfc = length(unique(spikeind_pfc));
    nCells_ca1 = length(unique(spikeind_ca1));
    phase_ca1 = zeros(50,length(lfpCA1));
    phase_pfc = zeros(50,length(lfpPFC));
    whlrld = whlrld_all{sess};
%     bins = [1, 1.2, 1.4, 1.6, 1.8, 2];
    bins = [1.1 1.4 1.6 1.9];
    nbins = length(bins)-1;
    epoch = cell(nbins,1);
    
    % Start the phase vectors
    parfor ii = 1:50
        phase_ca1(ii,:) = angle(hilbert(eegfilt(lfpCA1,srate_lfp,ii,ii+1)));
        phase_pfc(ii,:) = angle(hilbert(eegfilt(lfpPFC,srate_lfp,ii,ii+1)));
    end
    
    clear lfpPFC lfpCA1 whlrl_speed_all
    
    % Get spike timings in seconds for PFC pyr cells
    for ncell = 1:nCells_pfc
        cell_spks_pfc{ncell} = vertcat(spkts_pfc_t{sess}{:,ncell});
    end
    
    for ncell = 1:nCells_ca1
        cell_spks_ca1{ncell} = vertcat(spkts_ca1_t{sess}{:,ncell});
    end
    
    clear pall kappall RallHH RallHP RallPH RallPP
    %     Rall = nan(length(pyr_cells),50,5);
    
    % Gets all time spent in each bin
    for bin = 1:nbins
        ind = find(whlrld(:,7)>bins(bin) & whlrld(:,7)<bins(bin+1));
        epoch{bin} = round((ind/srate_bhv)*srate_lfp);
        epoch{bin} = upsample(epoch{bin},32);
        
        for k = find(epoch{bin}>0)
            for jj = 1:31
                epoch{bin}(k+jj) = epoch{bin}(k)+jj;
            end
        end
    end
    
    % Compute PLV for each bin
    for j = 1:nCells_ca1
        spks_indices=round(cell_spks_ca1{1,j}*srate_lfp);
        
        if(isempty(spks_indices))
            continue
        end
        
        for bin = 1:nbins
            
            if(length(intersect(spks_indices,epoch{bin}))<10)
                continue
            end
            
            for ii = 1:50
                
                [meanangle,meanvectorlength,angleSD,CI,kappa]= ...
                    anglemean(phase_ca1...
                    (ii,intersect(spks_indices,epoch{bin})));
                
                [meanangle2,meanvectorlength2,angleSD2,CI2,kappa2]= ...
                    anglemean(phase_pfc...
                    (ii,intersect(spks_indices,epoch{bin})));
                
                RallHH(j,ii,bin)=meanvectorlength;
                RallPH(j,ii,bin)=meanvectorlength2;
            end
        end
    end
    Rsess_all_HH{sess} = RallHH;
    Rsess_all_PH{sess} = RallPH;
    
    for j = 1:nCells_pfc
        spks_indices=round(cell_spks_pfc{1,j}*srate_lfp);
        
        if(isempty(spks_indices))
            continue
        end
        
        for bin = 1:nbins
            
            if(length(intersect(spks_indices,epoch{bin}))<10)
                continue
            end
            
            for ii = 1:50
                [meanangle,meanvectorlength,angleSD,CI,kappa]= ...
                    anglemean(phase_ca1...
                    (ii,intersect(spks_indices,epoch{bin})));
                
                [meanangle2,meanvectorlength2,angleSD2,CI2,kappa2]= ...
                    anglemean(phase_pfc...
                    (ii,intersect(spks_indices,epoch{bin})));
                
                RallHP(j,ii,bin)=meanvectorlength;
                RallPP(j,ii,bin)=meanvectorlength2;
            end
        end
    end
    Rsess_all_HP{sess} = RallHP;
    Rsess_all_PP{sess} = RallPP;
end

disp('end')

