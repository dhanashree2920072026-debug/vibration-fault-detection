// Machine Fault Detection Using Vibration Signal Analysis and Fourier Transform
// Scilab script that matches the project website (same signals, same results)
// Compares a healthy machine with four fault types:
// Unbalance, Misalignment, Bearing defect, Mechanical looseness

clc; clear; xdel(winsid());

// ---------------- 1. Parameters ----------------
fs = 1024;              // sampling frequency (Hz)
N  = 1024;              // number of samples (1 second of data)
T  = N / fs;            // duration (s)
t  = (0:N-1) / fs;      // time vector
fr = 20;                // shaft frequency, 1X (Hz)

s  = 1;                 // fault severity: 0.5 = mild, 1 = moderate, 1.5 = severe

names = ["Healthy", "Unbalance", "Misalignment", "Bearing defect", "Mech. looseness"];
short = ["Healthy"; "Unbal."; "Misal."; "Bearing"; "Loose"];

// ---------------- 2. Generate the signals ----------------
x = zeros(5, N);

// 1) Healthy: only the shaft frequency
x(1,:) = sin(2*%pi*fr*t);

// 2) Unbalance: the 1X peak grows
x(2,:) = (1 + 1.2*s) * sin(2*%pi*fr*t);

// 3) Misalignment: extra 2X and 3X components
x(3,:) = sin(2*%pi*fr*t) + 0.8*s*sin(2*%pi*2*fr*t) + 0.4*s*sin(2*%pi*3*fr*t);

// 4) Bearing defect: bursts at 70 Hz fault rate ringing at a 150 Hz resonance
x(4,:) = sin(2*%pi*fr*t) + 1.2*s*0.5*(1 + cos(2*%pi*70*t)) .* sin(2*%pi*150*t);

// 5) Mechanical looseness: 0.5X sub-harmonic plus 2X, 3X, 4X harmonics
x(5,:) = sin(2*%pi*fr*t) + s*( 0.4*sin(2*%pi*(fr/2)*t) + 0.6*sin(2*%pi*2*fr*t) ...
         + 0.5*sin(2*%pi*3*fr*t) + 0.4*sin(2*%pi*4*fr*t) );

// ---------------- 3. Energy, power, RMS and FFT ----------------
E = zeros(1,5); P = zeros(1,5); R = zeros(1,5); nump = zeros(1,5);
f   = (0:N/2-1) * fs / N;           // frequency axis (Hz)
mag = zeros(5, N/2);                // magnitude spectra
peaklist = list();

for i = 1:5
    E(i) = sum(x(i,:).^2);          // Energy  E = sum of x[n]^2
    P(i) = E(i) / N;                // Average power  P = E/N
    R(i) = sqrt(P(i));              // RMS = sqrt(P)

    X = fft(x(i,:));
    m = 2*abs(X(1:N/2)) / N;
    m(1) = m(1) / 2;
    mag(i,:) = m;

    // significant peaks: local maxima above 0.1 amplitude
    pk = [];
    for k = 2:N/2-1
        if m(k) > 0.1 & m(k) >= m(k-1) & m(k) >= m(k+1) then
            pk = [pk, f(k)];
        end
    end
    peaklist(i) = pk;
    nump(i) = length(pk);
end

// Parseval's theorem check (power from the spectrum should equal P)
Pfft = sum(abs(fft(x(1,:))).^2) / N^2;
mprintf("\nParseval check (healthy): time-domain P = %.4f, spectrum P = %.4f\n", P(1), Pfft);

// ---------------- 4. Decision ----------------
// ABNORMAL if more than one significant peak OR power > 1.2 x healthy power
status = emptystr(5,1);
for i = 1:5
    if (nump(i) > 1) | (P(i) > 1.2*P(1)) then
        status(i) = "ABNORMAL";
    else
        status(i) = "NORMAL";
    end
end

// ---------------- 5. Results table ----------------
mprintf("\nSeverity s = %.1f\n", s);
mprintf("%-18s %8s %10s %8s %10s %6s  %-9s %s\n", "Condition", "RMS", "Energy", "Power", "x healthy", "Peaks", "Status", "Peak frequencies (Hz)");
for i = 1:5
    mprintf("%-18s %8.3f %10.1f %8.3f %9.2fx %6d  %-9s %s\n", names(i), R(i), E(i), P(i), P(i)/P(1), nump(i), status(i), strcat(string(peaklist(i)), ", "));
end

// ---------------- 6. Plots: time domain and FFT ----------------
scf(0);
for i = 1:5
    subplot(5, 2, 2*i-1);
    if i == 1 then
        plot(t(1:205), x(i,1:205), "g");    // healthy in green
    else
        plot(t(1:205), x(i,1:205), "r");    // faults in red
    end
    a = gca(); a.data_bounds = [0 -4; 0.2 4];
    xtitle(names(i) + " - time domain", "Time (s)", "x(t)");

    subplot(5, 2, 2*i);
    plot2d3(f(1:251), mag(i,1:251));
    a = gca(); a.data_bounds = [0 0; 250 4];
    xtitle(names(i) + " - magnitude spectrum", "Frequency (Hz)", "|X(f)|");
end

// ---------------- 7. Comparison bar charts ----------------
scf(1);
vals = list(E, P, R);
ttl  = ["Energy  E = sum(x^2)", "Average power  P = E/N", "RMS = sqrt(P)"];
for j = 1:3
    subplot(1, 3, j);
    bar(1:5, vals(j));
    a = gca();
    a.x_ticks = tlist(["ticks","locations","labels"], (1:5)', short);
    xtitle(ttl(j));
end
