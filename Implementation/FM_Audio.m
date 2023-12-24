clc
clear

%% Modulation parameters

% Carrier amplitude
Ac = 1;
% Carrier frequency
fc = 40000; % 40 kHz
% Reading audio file
[audio, faudio] = audioread("audio1.mp3");
% Frequency sensitivity
kf = 100; % 100 Hz/V
% Simulation sampling frequency
fs = 8 * fc;

%% Modulation

% Preparing audio for modulation
audio = transpose(audio);
audio = audio(1, :);

% Creating message signal for audio
message = resample(audio, fs, faudio);
t = 0 : 1 / fs : (length(message) - 1) / fs;

% Frequency modulation
carrier = cos(2*pi*fc*t);
integralm = cumsum(message) / fs;
modulated = Ac * cos(2*pi*fc*t + 2*pi*kf*integralm);

% Plotting message, carrier, and modualted signals
figure(1)
subplot(3, 1, 1);
plot(t, message);
title("Message signal")
ylabel("m(t)")
subplot(3, 1, 2);
plot(t, carrier);
title("Carrier signal")
ylabel("c(t)")
subplot(3, 1, 3)
plot(t, modulated);
title("Modulated signal")
ylabel("s(t)")
xlabel("time (s)")

%% Demodulation

% Obtaining instantaneous frequency
x = hilbert(modulated);
% unwrapping angle for differentiation
phase = unwrap(angle(x));
freq = [0, diff(phase)];

% Filtering
f01 = 100;
f02 = 20000;
filtered = bandpass(freq, [f01, f02], fs);

% Amplifying and normalizing
normalized = fc/kf * filtered;
normalized = detrend(normalized);
normalized(normalized > 1) = 1;
normalized(normalized < -1) = -1;

% Removing DC component
output = normalized;

% Plotting modulated and demodulated signals
figure(2)
subplot(2, 1, 1)
plot(t, modulated);
title("Modulated signal")
ylabel("s(t)")
subplot(2, 1, 2)
plot(t, output);
title("Demodulated signal")
ylabel("m'(t)")
xlabel("time (s)")

% Retreiving audio with original sampling frequency
output = resample(output, faudio, fs);
player = audioplayer(output, faudio);
player.play;

