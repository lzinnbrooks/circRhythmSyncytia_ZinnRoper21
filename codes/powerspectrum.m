%Compute the power spectrum for a signal using fast Fourier transform
function [f,power] = powerspectrum(solx,solt)

T=solt(end)-solt(1); % length of time of simulation
solx=solx-mean(solx); %shift to remove constant frequency component

N=length(solx);
Fs=N/T; %sampling frequency

NFFT = 2^nextpow2(N); % Next power of 2 from length of y
Y = fft(solx,NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2+1);
psd=2*abs(Y(1:NFFT/2+1));
power=psd.^2/2;
end
