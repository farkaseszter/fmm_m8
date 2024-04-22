% ========================== Import data =========================== %
clc
clearvars
close all
format long

%% USER INPUTS

filename = "\pData_ref90_p";   % Set a consistent name for the files please
path_int = "E:\Dokumentumok\MSc_1\FMM\measurements\04_16";   % Path for the folder where you have your data
N = 90000;   % Number of datapoints
fs = 1000;   % Sampling frequency
maxPhi = 180;   % Largest value for phi
dPhi = 15;   % The size of the angle increments between 2 measurements
numAvg = 200;   % Number of averages in PSD
plot_data = 60;   % Which data do you want to plot? Give value in degrees
% IF YOU CHANGE THESE INPUTS AFTER RUNNING THE CODE AT LEAST ONCE, DELETE THE EXISTING DATA.MAT
% FILE BEFORE RERUNNING!

%% Defining variables
N_meas = maxPhi/dPhi+1;   % Number of measurements
% frequency
df = fs/N;
f = 0:df:(fs/2-df);
% time
dt = 1/fs;
tLength = N*dt;
t = 0:dt:(tLength-dt);

%% File reader

N_lines = N + 2;
N_chan = 2;
linesToSkip = 2;
fileFormat = '%s';

if exist("data.mat","file")==2
    load data.mat
else
    press = zeros(N, N_chan, N_meas);
    path = string(N_meas);
    for i = 1:N_meas
     path(i) = path_int + filename + num2str(maxPhi - ((i-1)*dPhi) + ".dat");
     fileID = fopen(path(i), "r");
     for j = 1:N_lines
         line = fgetl(fileID);
         if j <= linesToSkip
             continue;
         else
             p1 = str2double(extractBefore(line, '	'));
             p2 = str2double(extractAfter(line, '	'));
             press(j-2, 1, i) = p1;
             press(j-2, 2, i) = p2;
         end
     end
    end
    save data.mat
end

%% PSD calculation

pFFT = zeros(N_lines-2, N_chan, N_meas);
pFFT = fft(press);

pPSD = zeros(N_lines-2, N_chan, N_meas);

for i = 1:N_meas
        
        pPSD(end, 1, i) = (1/(fs*N))*abs(pFFT(end, 1, i)).^2;
        pPSD(end, 2, i) = (1/(fs*N))*abs(pFFT(end, 2, i)).^2;
        pPSD(2:(end-1), 1, i) = (2/(fs*N))*abs(pFFT(2:(end-1), 1, i)).^2;
        pPSD(2:(end-1), 2, i) = (2/(fs*N))*abs(pFFT(2:(end-1), 2, i)).^2;
end

%Truncating PSD at Nyquist frequency
pPSD = resize(pPSD, length(pPSD)/2);
fAper = resize(f, length(f)/2);


%% averaging FFT output

pPSDAvg = zeros(length(pPSD)/numAvg-1, N_chan, N_meas);
fAvg = 0:df*numAvg*2:fs-2*(df*numAvg*2);
fAvg = resize(fAvg,fs/2);

for i = 1:N_meas
    for j = 2:length(pPSDAvg)
        sumAvg1 = sum( pPSD( (j-1)*numAvg:j*numAvg, 1, i ) );
        sumAvg2 = sum( pPSD( (j-1)*numAvg:j*numAvg, 2, i ) );
        avg1 = sumAvg1/numAvg;
        avg2 = sumAvg2/numAvg;
        pPSDAvg(j-1, 1, i) = avg1;
        pPSDAvg(j-1, 2, i) = avg2;
    end
end
pPSDAvg = resize(pPSDAvg,fs/2);


%% Determining vortex shedding frequency

f_vs_all = zeros(N_chan, N_meas);
for i = 1:N_meas
    for j = 1:N_chan
        f_vs_all(j,i) = find(pPSDAvg(:,j,i)==max(pPSDAvg(2:end,j,i)))*df*numAvg*2;
    end
end
f_vs = mode(f_vs_all(2,:))

%% Calculating phase shift

phases = zeros(N_meas,1);
for i = 1:N_meas
    a_sum=0;
    b_sum=0;
    for k=0:(N-1) 
        time = k*dt;
        a_sum = a_sum + press(k+1,1,i)*cos(f_vs*time);%perform fourier correlation with ONLY THE RELEVANT FREQUENCY
        b_sum = b_sum + press(k+1,1,i)*sin(f_vs*time);
    end
    a = a_sum*2/N;
    b = b_sum*2/N;
    phase = atan2(a,b);  %phase is arctangent of fourier coefficients
    
    
    a_sum=0;  
    b_sum=0;
    for k=0:(N-1)
        time = k*dt;
        a_sum = a_sum + press(k+1,2,i)*cos(f_vs*time);%perform fourier correlation with ONLY THE RELEVANT FREQUENCY
        b_sum = b_sum + press(k+1,2,i)*sin(f_vs*time);
    end
    a = a_sum*2/N;
    b = b_sum*2/N;
    phase1 = atan2(a,b);
    
    phases(i)=(phase-phase1);
end
phases
phase = phases(fix(N_meas/2)+1)


%% PLOTS
%time series

num = 1 + 1/30*(180 - plot_data);
fig1 = figure("Name","Plots","WindowState","maximized");
subplot(2, 2, 1)
plot(t, press(:, 1, num))
xlabel('t [s]');
ylabel ('p_1 [Pa]')
title("Time series of " + num2str(plot_data) + "° pressure signals");
subplot(2, 2, 3)
plot(t, press(:, 2, num))
xlabel('t [s]');
ylabel ('p_2 [Pa]');
title("Time series of  270° pressure signals")

% PSD
subplot(2, 2, 2)
plot(fAvg, pPSDAvg(:, 1, num))
xlabel('f [Hz]');
ylabel ('avg(PSD(p_1))');
xlim([0,fs/2])
title("PSD of " + num2str(plot_data) + "° pressure signals")
subplot(2, 2, 4)
plot(fAvg, pPSDAvg(:, 2, num))
xlabel('f [Hz]');
ylabel ('avg(PSD(p_1))');
xlim([0,fs/2])
title("PSD of 270° pressure signals")