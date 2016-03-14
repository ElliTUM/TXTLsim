function [fftE,dExpdT,dSimdT,tExpFFT,tSimFFT] = fftError(expDataInFFT,simDataInFFT,expTimeInFFT,simTimeInFFT,cutFFT)
%FFTERROR This function computes the error between the FFT of a
%experimental and simulated data set.

%derive simulated data
dSimdT=[0];
for i=1:length(simDataInFFT)-1
    dSimdT(i,1)=(simDataInFFT(i+1,1)-simDataInFFT(i,1))/(simTimeInFFT(i+1,1)-simTimeInFFT(i,1));
end

% fast fourier transform expData
x = fft(expDataInFFT);
%figure(1); plot(abs(X)); title('Amplitudes as a function of frequency');

% create a filter mask
N = length(expDataInFFT);
mask = zeros(1, N);
f = 0:N/2;
sigmaf = 50;
mask(1:N/2+1) = exp(-(f/(2*sigmaf)).^2);
mask(N:-1:N/2+2) = mask(2:N/2);
% subplot(3,1,1); plot(mask); title('Low-pass mask');
mask=mask';
% Multiply the FT by the mask
xFilt = x .* mask;
% subplot(3,1,2); plot(abs(Xfilt)); title('Weighted transform');

% Transform back
filteredExpData = real(ifft(xFilt)); %somehow there is a small imaginary part left
%subplot(3,1,3); plot(xfilt); title('Low-pass filtered signal');

dExpdT=[0];
for i=1:length(filteredExpData)-1
    dExpdT(i,1)=(filteredExpData(i+1,1)-filteredExpData(i,1))/ ...
        (expTimeInFFT(i+1,1)-expTimeInFFT(i,1));
end


dExpdT=dExpdT(cutFFT:end-cutFFT,1); 
dSimdT=dSimdT(cutFFT:end-cutFFT,1);
tExpFFT=expTimeInFFT(cutFFT:end-cutFFT-1,1); 
tSimFFT=simTimeInFFT(cutFFT:end-cutFFT-1,1); %should be the same as tExpFFT

fftE = sum((dSimdT-dExpdT).^2); %least squared error of derivatives
end

