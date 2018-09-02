function [coloredNoise] = generateNoise(nrChan,nrSamples,InvFrequencyPower)
%GENERATENOISE Summary of this function goes here
%   Detailed explanation goes here

whiteNoise   = randn(nrChan,nrSamples+mod(nrSamples,2));
% Fourier Transform 
NoiseFT      = fft(whiteNoise').';
NoiseFT      = NoiseFT(:,1:(nrSamples+mod(nrSamples,2))/2+1);  

% create desired noise type
fIndex       = repmat([1:size(NoiseFT,2)],nrChan,1);
NoiseFT      = NoiseFT./sqrt(fIndex.^InvFrequencyPower);
NoiseFT      = [NoiseFT,conj(NoiseFT(:,end-1:-1:2))];

% Inverse Fourier Transform
coloredNoise = real(ifft(NoiseFT.'))';
coloredNoise = coloredNoise(:,1:nrSamples);
coloredNoise = coloredNoise-mean(coloredNoise,2);
coloredNoise = coloredNoise./std(coloredNoise,1,2);
end
