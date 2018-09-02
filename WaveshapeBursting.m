% Script illustrating the Robustness of biphase estimates
% for different levels of pink and brown noise
% Analyzing the Waveshape of Brain Oscillations with Bicoherence:
% Simulations of
% - Section 1.1.4
% - Figure 2
%
%
% Copyright (c) [2018] [Sarah Bartz-Schaechtelin]


close all
clear
set(0,'defaultfigurecolor',[1 1 1])
s             = RandStream('mrg32k3a');
RandStream.setGlobalStream(s)
s.Substream   = 2;
colors        = get(0,'DefaultAxesColorOrder');
nrRepetitions = 50;
fs            = 1000;
segleng       = fs;
T             = 300;
mf            = 50;
xderlims      = [-0.05,0.05;-0.1,0.1;-0.05,0.05;-0.1,0.1;-0.1,0.1];
%cm=[255,255,255]/255;

% Values of InvFreqPows vs. noise color:
% 'brown'    =  2
% 'pink'     =  1
% 'white'    =  0
% 'blue'     = -1
% Entries of InvFreqPows are not limited to the above examples
InvFreqPows  = [2,1];
snr_values   = logspace(-2,0,5);
t            = [0:(1/fs):T-(1/fs)];
f0           = 10;
phis         = [0,pi/2,pi,3*pi/2,-3*pi/4]+1.5*pi;
nrSNR        = length(snr_values);
nrNoiseTypes = length(InvFreqPows);

nrHarmonics        = 200;
all_harmonics      = zeros(length(phis),nrHarmonics,length(t));
periodic_signals   = zeros(length(phis),length(t));

% burst
for phix=1:length(phis)
    phi0 = phis(phix);
    for h = 1:nrHarmonics
        all_harmonics(phix,h,:) = (1/(h))*cos(2*pi*h*f0*t+(h-1)*phi0);
    end
    c_sig                    = sum(squeeze(all_harmonics(phix,:,:)),1);
    ms                       = max(abs(c_sig));
    periodic_signals(phix,:) = c_sig/ms;
end

burst      = periodic_signals(:,1:1000);
burst(3,:) = periodic_signals(3,31:1030);
burst(2,:) = periodic_signals(2,36:1035);
burst(1,:) = periodic_signals(1,76:1075);

data = [];
for k = 1:T*2
    data = [data,zeros(length(phis),100+randi([0,400])),(10+randn*2).*burst(:,1:100*randi([length(phis),10]))];
end
psm = data(:,1:T*fs);
periodic_signals = (psm-mean(psm,2));
periodic_signals = periodic_signals./(std(psm,[],2));

biphase_save = zeros(nrNoiseTypes,nrRepetitions,length(phis),nrSNR);
bimag_save   = zeros(nrNoiseTypes,nrRepetitions,length(phis),nrSNR);
signal_save  = zeros(nrNoiseTypes,length(phis),nrSNR,size(periodic_signals,2));

for sigix=1:length(phis)
    target=periodic_signals(sigix,:);

    for noise_type_idx = 1:nrNoiseTypes
        InvFreqPow = InvFreqPows(noise_type_idx);

        for rix=1:nrRepetitions
            %generate noise for each SNR
            noise_signals = generateNoise(nrSNR,length(target),InvFreqPow);
            power_ratio   = mean(target.^2)./mean(noise_signals.^2,2);
            asnr          = sqrt(snr_values'./power_ratio);
            asignal       = repmat(target,size(snr_values,2),1).*repmat(asnr,1,size(target,2));
            ns_signals    = asignal+noise_signals;


            if(rix==1)
                signal_save(noise_type_idx,sigix,:,:) = ns_signals;
            end

            [bs,bsnr]                 = compute_bispectrum(ns_signals',fs,fs/2,fs,50);
            bicoh                     = bs./bsnr;
            bimag_save(noise_type_idx,rix,sigix,:)   = abs(bicoh(:,11,11))';
            biphase_save(noise_type_idx,rix,sigix,:) = (wrapTo2Pi(angle(bicoh(:,11,11)).'));
        end
    end
end

correction        = [0.5*pi,0,1.5*pi,pi,1.25*pi];
alter             = [-2*pi,-2*pi,2*pi,0,+2*pi];

for noise_type_idx = 1:nrNoiseTypes
    noise_biphase_save = squeeze(biphase_save(noise_type_idx,:,:,:));
    biphase_corrected  = noise_biphase_save;

    for sigix = 1:length(phis)

        bps                          = squeeze(noise_biphase_save(:,sigix,:));
        val1                         = bps+alter(sigix);
        val2                         = bps;
        d1                           = abs(val1-correction(sigix));
        d2                           = abs(val2-correction(sigix));
        ixi                          = d2>d1;
        bps(ixi)                     = val1(ixi);
        biphase_corrected(:,sigix,:) = bps;
    end
    biphase_save(noise_type_idx,:,:,:) = biphase_corrected;
end

%% PLOTTING (setting specific)
offsets    = [0,2,4,6,8];
snr_colors = [145,209,229;119,171,191;80,114,130]/255;
show_idx   = [1,3,5];
figure
spIdx      = [1,3];
set(0,'DefaultAxesColorOrder',snr_colors);

for noise_type_idx = 1:nrNoiseTypes
    example_signals = squeeze(signal_save(noise_type_idx,1,:,1:5000));
    es              = example_signals;
    es              = es-mean(es,2);
    es              = es./max(abs(es),[],2);

    subplot(2,nrNoiseTypes,spIdx(noise_type_idx))
    plot(es(show_idx,:)'+repmat(offsets(show_idx)',1,size(es,2))')
    yticks([])
    box off
    if (noise_type_idx == 1)
        l1 = legend('1','0.1','0.01');
    end
end

title(l1,'Signal-To-Noise Ratio (SNR)')
set(l1,'Position',[0.6277    0.6179    0.2312    0.2620])
set(groot,'defaultAxesColorOrder','remove')

subplot(2,2,4)
noise_colors = [226,131,214;115,99,87]/255;
for noise_type_idx = 1:nrNoiseTypes
    biphase_corrected = squeeze(biphase_save(noise_type_idx,:,1,:));
    correction        = 0.5*pi;
    bias              = mean(biphase_corrected-correction);
    variance          = var(biphase_corrected);
    mse               = mean((biphase_corrected-correction).^2);
    rmse              = sqrt(mse);
    plot(log10(snr_values),rmse','-o','Color',noise_colors(noise_type_idx,:),'LineWidth',1)
    ylim([-0.5,2])
    hold on
end
l2 = legend('Brown Noise','Pink Noise');

subplot(2,2,1)
ylabel({'Brown Noise','','normalized amplitude'})
xticklabels([0,1,2,3,4,5])
xlabel('Time (s)')

subplot(2,2,3)
ylabel({'Pink Noise','','normalized amplitude'})
xticklabels([0,1,2,3,4,5])
xlabel('Time (s)')

subplot(2,2,4)
title({'Root-Mean-Square Error (RMSE) of ';'Biphase Estimates'})
xlabel('Log10 (SNR)')
ylabel('RMSE (rad)');
