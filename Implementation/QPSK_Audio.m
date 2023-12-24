clc
clear

%% Modulation parameters

% Sampler frequency
fsampler = 20000; % 20 kHz
% Quantizing bits
b = 6;
% Bit period
Tb = 0.0001; % 100 microseconds
% Carrier frequency
fc = 40000; % 40 kHz
% Read audio from file
[audio, faudio] =  audioread("audio4.mp3");
% Simulation sampling frequency
fs = 8 * fc;

%% Pulse code modulation

% Preparing audio for modulation
audio = transpose(audio);
audio = audio(1, :);

% Creating message signal for audio
message = audio;
txlimit = (length(message) - 1) / faudio;
t1 = 0 : 1 / faudio : txlimit;

% Sampling
sampled = resample(message, fsampler, faudio);
t2 = 0 : 1 / fsampler : (length(sampled) - 1) / fsampler;

% Quantizing
quantized = floor((0.99 * normalize(sampled, "range")) * 2 ^ b);

% Encoding
encoded = split(strjoin(string(dec2bin(quantized, b)), ""), "");
encoded = encoded(2 : length(encoded) - 1); % remove empty strings
encoded = transpose(str2double(encoded));

% Plot audio, sampled audio, and quantized audio
figure(1)
subplot(3, 1, 1)
plot(t1, message)
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Message signal")
ylabel("m(t)")
subplot(3, 1, 2)
stem(t2, sampled, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Sampled signal")
subplot(3, 1, 3)
stem(t2, quantized, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Quantized signal")
xlabel("time (s)")

% Clearing memory for better performance
clear sampled
clear quantized

%% QuadriPhase shift keying

% Symbolizing encoded signal
symbolized = reshape(encoded, [2, length(encoded) / 2]);
symbolized = string(symbolized);
symbolized = symbolized(1, :) + symbolized(2, :);
symbolized = bin2dec(symbolized);

% Changing the map causes distortion
map = [2, 3, 1, 4];
% map = [3, 1, 2, 4];
symbolized = map(symbolized + 1);

T = 2 * Tb;
t3 = 0 : 1 / fs : length(symbolized) * T;
t3 = t3(1:length(t3)-1);

% Serializing symbolized signal (non-return to zero)
serial = round(resample(symbolized, fs, ceil(1 / T), 0));

% Modulating the serialized signal
carrier = cos(2*pi*fc*t3);
% Note: filter bitstream first
modulated = cos(2*pi*fc*t3 + (2*serial - 1) * pi/4);

% Plotting serialized signal and modulated signal
figure(2)
subplot(3, 1, 1)
plot(t3, serial)
title("Serialized signal")
subplot(3, 1, 2)
plot(t3, carrier)
title("Carrier signal")
ylabel("c(t)")
subplot(3, 1, 3)
plot(t3, modulated)
title("Modulated signal")
ylabel("s(t)")
xlabel("time (s)")

% Clearing memory for better performance
clear encoded
clear serial

%% Demodulation

% Coherent detector
detected1 = modulated .* carrier;
detected2 = modulated .* imag(hilbert(carrier));

% Low pass filter and normalizer
f0 = 1 / T;
[bl, al] = butter(1, f0 / (fs / 2), "low");
filtered1 = filter(bl, al, detected1);
filtered1 = normalize(filtered1, "range");
filtered2 = filter(bl, al, detected2);
filtered2 = normalize(filtered2, "range");

% Sampling
sampled1 = filtered1(fs * T / 4 : fs * T : length(filtered1));
sampled2 = filtered2(fs * T / 4 : fs * T : length(filtered2));

% Decision making
decoded1 = sampled1;
decoded1(decoded1 < 0.5) = 0;
decoded1(decoded1 >= 0.5) = 1;
decoded2 = sampled2;
decoded2(decoded2 < 0.5) = 0;
decoded2(decoded2 >= 0.5) = 1;

% Plotting modulated signal and filtered output of coherent detector
figure(3)
subplot(3, 1, 1)
plot(t3, modulated);
title("Modulated signal")
subplot(3, 1, 2)
ylabel("s(t)")
plot(t3, filtered1);
title("First sampled signal")
subplot(3, 1, 3)
plot(t3, filtered2);
title("Second sampled signal")
xlabel("time (s)")

% Clearing memory for better performance
clear carrier
clear modulated
clear detected1
clear detected2
clear filtered1
clear filtered2
clear sampled1
clear sampled2

% Combining and deserializing decoded signals
deserialized = string(decoded1) + string(decoded2);
deserialized = [deserialized{:}];
deserialized = mat2cell(deserialized, 1, b * ones(1, length(deserialized) / b));
deserialized = transpose(bin2dec(deserialized));

% Dequantizing deserialized signal to obtain original message
output = (deserialized - (2 ^ (b-1) - 0.5)) / (2 ^ (b-1));

% Plotting decoded signal and demodulated signal
figure(4)
subplot(2, 1, 1)
stem(t2, deserialized, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Combined deserialized signal")
subplot(2, 1, 2)
plot(t2, output, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Demodulated signal")
ylabel("m'(t)")
xlabel("time (s)")

player = audioplayer(output, fsampler);
player.play