clc
clear

%% Modulation parameters

% Carrier amplitude
Ac = 1;
% Carrier frequency
fc = 40000; % 40 kHz
% Message tone amplitude
Am = 1;
% Message tone frequency
fm = 1000;
% Frequency sensitivity
kf = 10000; % 100 Hz/V
% Simulation sampling frequency
fs = 8 * fc;

%% Modulation

% Creating message signal for tone
t = 0 : 1 / fs : 0.005;
message = Am * cos(2*pi*fm*t);

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

% Amplifying and normalizing
normalized = fc/kf * freq;
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

