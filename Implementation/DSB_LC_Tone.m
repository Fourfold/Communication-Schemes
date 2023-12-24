clc
clear

%% Modulation Parameters

% Carrier amplitude
Ac = 1;
% Carrier frequency
fc = 40000; % 40 kHz
% Message tone amplitude
Am = 1;
% Message tone frequency
fm = 1000;
% Modulation index
ka = 0.75 / Am;
% Simulation sampling frequency
fs = 8 * fc;

%% Modulation

% Creating message signal for tone
t = 0 : 1 / fs : 0.005;
message = Am * cos(2*pi*fm*t);

% Amplitude Modulation
carrier = Ac * cos(2*pi*fc*t);
modulated = (1 + ka * message) .* carrier;

% Plotting message, carrier, and modulated signals
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
xlabel ("time (s)")

%% Demodulation (Envelope Detector)

% Rectification
x1 = abs(modulated);

% Low pass filter
f0 = 20000;
[bl, al] = butter(5, f0 / (fs / 2), "low");
x2 = filter(bl, al, x1);

% Remove DC component
x3 = detrend(x2);

% Restore amplitude
output = x3 * 2;

% Plotting modulated, rectified, and demodulated signals
figure(2)
subplot(3, 1, 1)
plot(t, modulated)
title("Modulated signal")
ylabel("s(t)")
subplot(3, 1, 2)
plot(t, x1)
title("Rectified signal")
subplot(3, 1, 3)
plot(t, output)
title("Demodulated signal")
ylabel("m'(t)")
xlabel("time (s)")

