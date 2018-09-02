

% Script testing the signal extraction performances of different methods
% Analyzing the Waveshape of Brain Oscillations with Bicoherence:
% Simulations of
% - Section 2.1.3 
% - Figure 4
%

close all
clear
set(0,'defaultfigurecolor',[1 1 1])


% To run this code on other methods as well, they need to be downloaded
% and added to the MATLAB path
% SSD  : https://github.com/svendaehne/matlab_SPoC/SSD/
% TDSep: We cannot distribute the source code and it does not seem to be
%        available to the online community at the moment. 
%        You can contact Andreas Ziehe and ask him to provide you with a 
%        MATLAB implementation:
%        Email:Andreas.Ziehe@tu-berlin.de         
%
%methods     = {'tdsep','ssd','ssdf','hpmax'};
methods     = {'hpmax'};

% simulation settings
nrRepetitions       = 30;
fs                  = 1000; 
T                   = 300; 
mf                  = 50;
t                   = 0:(1/fs):T-(1/fs);
nrSensors           = 20;
nrNoise             = 20;
% Values of InvFreqPows vs. noise color:
% 'brown'    =  2
% 'pink'     =  1
% 'white'    =  0
% 'blue'     = -1
% Entries of InvFreqPows are not limited to the above examples
InvFreqPows  = [2,0,-1];
maxfreq      = 51;
alpha_range  = 10:12; 
snr_values   = logspace(-6,2,9);
nrSNR        = length(snr_values);
nrNoiseTypes = length(InvFreqPows);
nrMethods    = length(methods);
% settings for spectral computations
segleng  = 1000;
segshift = 1000;
epleng   = 1000;
% all performances are stored in this array
results  = zeros(nrNoiseTypes,nrMethods,nrSNR,nrRepetitions);
% settings for target oscillation
nr_harm    = 80;
phi0       = 3*pi;
% loop over noise types
for noise_type_idx = 1:nrNoiseTypes
    InvFreqPow = InvFreqPows(noise_type_idx);
    % loop over repetions
    for rep_idx = 1:nrRepetitions
        % loop over snr values
        for snrIdx = 1:nrSNR
            % create mixing matrix for noise sources
            An         = randn(nrSensors,nrNoise);
            sigToNoise = snr_values(snrIdx);
            % burst type 1: modified sawtooth       
            burst_1    = zeros(1,fs);            
            for harx = 1:nr_harm
                burst_1 = burst_1+(1/(harx^1.5))*sin(harx*10*2*pi*t(1:fs)+(harx-1)*phi0);
            end
            % burst type 2: normal sawtooth              
            burst_2 = sawtooth(2*pi*10*t(51:1050));           
            burst   = [burst_1;burst_2];
            data    = [];
            % create target bursting oscillator
            for idx = 1:T*2
                data = [data,zeros(1,100+randi([0,400])),(10+randn*2)*burst(randi([1,2]),1:100*randi([5,10]))];
            end
            target = data(1:T*fs);
            target = (target-mean(target));
            target = target/(std(target,[],2));
            % generate noise sources
            noise_signals   = generateNoise(nrNoise,length(target),InvFreqPow); 
            % mix noise sources
            mix_noise       = An*noise_signals;
            % check snr
            variance_ratio  = mean(mean(target.^2)./mean(mix_noise.^2,2));
            % adjust to target snr
            adjust_factor   = sqrt(sigToNoise/variance_ratio);
            adjusted_signal = repmat(target,nrSensors,1).*repmat(adjust_factor,nrSensors,length(target));
            % combine signal and mixed noise sources
            M = mix_noise+adjusted_signal;

            %% Demix signal 
            for methodIdx = 1:nrMethods 
                method = methods{methodIdx};
                
                switch method
                    case 'tdsep'
                        thresh          = .0001;
                        [C,D]           = tdsep3(M,[0:100],thresh);
                        W_method        = inv(C);
                        A_method        = C;
                        db              = (W_method*M);
                       
                    case 'ssdf'
                        freq_ssd2                    = [9 11; 7 13; 8 12];
                        data_raw                     = M';
                        Fs                           = fs;
                        [W2,Aret2,lambda,C_s,X_ssdf] = SSD(data_raw, freq_ssd2, Fs, [], []);                  
                        W_method                     = [W2'];
                        A_method                     = Aret2;                         
                        db                           = X_ssdf';

                    case 'ssd'
                        freq_ssd2                        = [9 11; 7 13; 8 12];
                        data_raw                         = M';
                        Fs                               = fs;
                        [W2, Aret2, lambda, C_s, X_ssd2] = SSD(data_raw, freq_ssd2, Fs, [], []);
                        trans_data2                      = (data_raw*W2);
                        W_method                         = [W2'];
                        A_method                         = Aret2;
                        db                               = [trans_data2]';

                     case 'hpmax'
                        alpha_range                            = [11];
                        topMany                                = 1;
                        df                                     = 8;
                        dStop                                  = 1;
                        [data_back_HPMax,top_filters,top_topos] = get_HPMAX_rev(M,segleng,segshift,epleng,alpha_range,mf,topMany,dStop,df);
                        W_method                                = top_filters{1};
                        A_method                                = top_topos{1};
                        db                                      = data_back_HPMax{1};
                end


                %% Best Selection                 
                results(noise_type_idx,methodIdx,snrIdx,rep_idx) = -inf;

                for ix = 1:size(A_method,2)

                    db_meth        = db(ix,:);
                    db_meth        = (db_meth-mean(db_meth));
                    db_meth        = db_meth/std(db_meth,[],2);
                    cdm1           = corrcoef(db_meth(1000:end),target(1000:end));
                    cdm1           = cdm1(1,2);
                    cdm2           = corrcoef(-db_meth(1000:end),target(1000:end));
                    cdm2           = cdm2(1,2);             
                    [corr_max,bix] = max([cdm1,cdm2]);
                    if(bix == 2)
                        db_meth = -db_meth;
                    end

                    if corr_max > results(noise_type_idx,methodIdx,snrIdx,rep_idx)
                        results(noise_type_idx,methodIdx,snrIdx,rep_idx) = corr_max;                    
                    end           
                end
            end
        end
    end
end


%% setting-specific plotting 
figure
set(gcf,'Position',[164,352,1117,324]);
rep_dim=3;
titles={'Brown Noise','White Noise','Blue Noise'};
if (nrMethods == 1)
    rep_dim=2;
end
for noise_type_idx = 1: nrNoiseTypes
    
    noise_type_results = squeeze(results(noise_type_idx,:,:,:));
    noise_type_means   = squeeze(mean(noise_type_results,rep_dim));
    noise_type_stds    = squeeze(std(noise_type_results,[],rep_dim));
    
    subplot(1,nrNoiseTypes,noise_type_idx)
    xsnr_values=repmat(snr_values,nrMethods,1);
    errorbar(log10(xsnr_values'),noise_type_means',noise_type_stds','-o')
    title(titles{noise_type_idx})
    xlabel('Log10 (SNR)')
    ylabel({'Correlation with','Target'})
    legend(methods)
    ylim([-0.2,1.2])
    xlim([-7,3])
end
    
   
