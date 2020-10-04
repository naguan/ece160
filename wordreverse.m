function [outputAud] = wordreverse(aud,fs)

%[aud, fs] = audioread('spoken-01.wav');
aud = audioread(aud);
cutoff = fs/4;
aud = lowpass(aud, cutoff, fs);
blkLength = floor(0.05*fs);
hopLength = floor(0.02*fs);
numBlk = floor(length(aud)/hopLength)-ceil(blkLength/hopLength);
time = ((0:numBlk-1)*hopLength+(blkLength*(hopLength/blkLength))/fs);
energy = zeros(numBlk, 1);
logEnergy = zeros(numBlk, 1);
flux = zeros(numBlk,1);
prevDFT = 0;
for n = 1:numBlk
    start = (n-1)*hopLength + 1;
    finish = min(length(aud), start + blkLength - 1 );
    hannwin = hann(finish - start + 1);
    window = aud(start:finish).^2;
    energy(n) = mean(window.*hannwin);
    logEnergy(n) = (log(energy(n)+0.00001));
    currDFT = fft(aud(start:finish));
    exponent = 0.001;
    differenceDFT = sum(abs(currDFT))-sum(abs(prevDFT));
    summation = sum(abs(differenceDFT).^exponent);
    flux(n) = summation.^(1/exponent);   
    prevDFT = currDFT;
end

tt = 0:1/fs:(length(aud)-1)*1/fs;
outputAud = [];
space = zeros(1000,1);
began = 0;
started = 0;
finished = 0;


for i = 1:numBlk
    if logEnergy(i) < -11 && began == 1
        began = 0;
        finished = i;
        wordSegStart = (started - 1)*hopLength + 1;
        wordSegEnd = wordSegStart + blkLength + (finished - started)*hopLength;
        outputAud = [aud(wordSegStart:wordSegEnd);outputAud];
        outputAud = [space; outputAud];
    end
    if logEnergy(i) > -7 && began == 0
        began = 1;
        started = i;
    end
end

soundsc(aud, fs);
pause(3);
soundsc(outputAud, fs);

subplot(6,1,1);
plot(tt, aud);
title("Audio");

subplot(6,1,2);
plot(time, energy); 
title("Avg Energy");

subplot(6,1,3);
plot(time,logEnergy);
title("Log Energy");

subplot(6,1,4);
plot(time, flux);
title("Spectral Flux");

subplot(6,1,5);
plot(outputAud);
title("Output Audio");

end
