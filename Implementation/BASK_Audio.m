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
[audio, faudio] =  audioread("audio3.mp3");
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
encoded = encoded(2 : length(encoded) - 1);
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

%% Amplitude shift keying

t3 = 0 : 1 / fs : length(encoded) * Tb;
t3 = t3(1:length(t3)-1);

% Serializing encoded signal (on-off signaling)
serial = resample(encoded, fs, ceil(1 / Tb));
serial(serial < 0.5) = 0;
serial(serial >= 0.5) = 1;

% Modulating the serialized signal
carrier = cos(2*pi*fc*t3);
modulated = serial .* carrier;

% Plotting serialized signal and modulated signal
figure(2)
subplot(2, 1, 1)
plot(t3, serial)
title("Serialized signal")
subplot(2, 1, 2)
plot(t3, modulated)
title("Modulated signal")
ylabel("s(t)")
xlabel("time (s)")

% Clearing memory for better performance
clear encoded
clear carrier
clear serial

%% Demodulation

% Envelope detector
x1 = abs(modulated);

% Low pass filter and normalizer
f02 = 1 / Tb;
[bl, al] = butter(1, f02 / (fs / 2), "low");
x1 = filter(bl, al, x1);
x1 = normalize(x1, "range");

% Sampling
x2 = x1(fs * Tb / 4 : fs * Tb : length(x1));

% Decision making
decoded = x2;
decoded(decoded < 0.5) = 0;
decoded(decoded >= 0.5) = 1;

% Plotting modulated signal and output of envelope detector and sampler
figure(3)
subplot(2, 1, 1)
plot(t3, modulated);
title("Modulated signal")
ylabel("s(t)")
subplot(2, 1, 2)
plot(t3, x1);
title("Detected signal")
xlabel("time (s)")

% Clearing memory for better performance
clear modulated
clear x1
clear x2

% Deserializing
deserialized = string(decoded);
deserialized = [deserialized{:}];
deserialized = mat2cell(deserialized, 1, b * ones(1, length(deserialized) / b));

% Deserializing decoded signal
deserialized = transpose(bin2dec(deserialized));

% Dequantizing deserialized signal to obtain original message
output = (deserialized - (2 ^ (b-1) - 0.5)) / (2 ^ (b-1));

% Plotting decoded signal, demodulated signal, and original message
figure(4)
subplot(2, 1, 1)
stem(t2, deserialized, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Deserialized signal")
subplot(2, 1, 2)
plot(t2, output, 'Marker', 'None')
xlim([-0.01 * txlimit, 1.01 * txlimit])
title("Demodulated signal")
ylabel("m'(t)")
xlabel("time (s)")

player = audioplayer(output, fsampler);
player.play
