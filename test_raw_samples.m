
function test_raw_samples(fbin,sampfreq);
% Test Front end

% open the input file with raw samples

% read a snapshot of sample with appropriated format
fseek(fbin,ceil(sampfreq*0),-1);
[data, cnt]= fread(fbin,ceil(sampfreq*0.1),'schar');

data = data(1:1:end);
%data = data -128;

figure(1);
%plot(data(1:ceil(sampfreq*0.005)),'ko-'), grid on;
plot(data(1:1000),'ko-'), grid on;
axis tight;aa=axis;axis([aa(1) aa(2) aa(3)-1 aa(4)+4])
title('Time Domain of the first 1000 samples');
xlabel('Num. of samples');ylabel('Amplitude');


x = linspace(-31,32,64);
figure(2),hist(data(1:ceil(sampfreq*0.02)),x);

% plot signal spectrum
numptsfft = ceil(sampfreq*0.09);%1048576;    
NFFT = 4096; 

[Pxx_L1,F_L1] = psd(data(100:numptsfft)-mean(data(100:numptsfft)),NFFT,sampfreq);
figure(3)
plot(F_L1/1e6,10*log10(abs(Pxx_L1)),'k'),grid on;
axis tight;aa=axis;axis([aa(1) aa(2) aa(3)-1 aa(4)+4])
title('L1 Spectrum');
xlabel('freq (MHz)');ylabel('Mag (dB)');








