%% Exercise 3: Convolution to estimate voltage response to a spike train
clear all, close all, clc
% A convolution is a mathematical operation on two functions 
% that expresses the amount of overlap of one function as 
% it is shifted over the other. 
% We will use convolution to estimate the voltage response of 
% a neuron to an incoming spike train

% Generate 3 seconds of a Poisson spike train with firing rate 20Hz.
targetRate = 20; % Measured in Hz.
N = 10000; % Measured in ms.
firingRate = 0;
p = 0;

while 1
    if abs( firingRate - targetRate ) < 0.1
        break
    end

    spikeTrain = rand(1, N) > (1 - p);
    firingRate = sum(spikeTrain == 1) / (N/1000);
    p = p + 0.00001;
    
    if p > 1
        disp('did not find a suitable p')
        break
    end 
end

% Construct a kernel, the response of the neuron to one spike at t=0. 
% We assume the neuron is linear, that is, 
% the response to multiple spikes is the sum of the responses to each
% individual spike. For the kernel, use an exponential with mean mu of 5ms. 
% Calculate the kernel for values between -50 and 50ms
mu = 5;
k = exppdf(-50:50, mu); 

% Find spike times first
ts = 1 : N; % time vector in ms
spkTimes = ts(spikeTrain);

% convolution to estimate the voltage response of the cell to the spike train
spikeTrain = double(spikeTrain);
estVoltage = conv(spikeTrain, k, 'same');

% plot
figure
subplot(2,1,1)
plot(-50:50, k)
box off
title('kernel')
xlabel('tau (ms)')

subplot(2,1,2)
plot(ts, estVoltage)
hold on
yLevel = 0.3;
plot(spkTimes, yLevel, 'k', 'Marker', '*', 'MarkerSize', 4)
box off
ylabel('voltage (a.u.)')
xlabel('time (ms)')
xlim([0 1000])
ylim([0 0.4])
legend('voltage', 'spike times')

%% Exercise 4: Convolution to detect edges in images
clear all, close all, clc
% In an image, edges are where the image is different from its neighbors. 
% Convolution in two dimensions is often used for edge detection in 
% image processing. The output of this operation is very similar to 
% the response of cells in primary visual cortex, 
% which respond selectively to oriented edges. 
% The following kernel will be zero in regions of the image where 
% neighboring pixels have the same value, and nonzero for edges

% Load image.(from http://www.bbc.co.uk/nature/life/Octopus)
load('octopus.mat')

figure
imagesc(octopus);
colormap gray

% Make the kernel.
k = [...
    0, 0, 0;...
    0, 1.125, 0 ;...
    0, 0, 0 ] - 0.125 * ones( 3, 3 );
% k = double( k );
octopusEdges = conv2( octopus, k, 'same' );

figure
imagesc(octopusEdges);
colormap gray

% Plot the absolute value of the convolution 
% (this will show both positive and negative edges as lighter than the background).
edgeImage = abs(octopusEdges);

figure
imagesc(edgeImage);
colormap gray
caxis([-20 100])

%% Exercise 5: Correlation to analyze premotor neural data
clear all, close all, clc
% We’ll use cross-correlation to analyze 
% the timing relation between the song of a juvenile zebra finch 
% and neural data recorded while he was singing.

load('MackeviciusData.mat')

T = length(song) / fs;
dt = 1 / fs;
k = length(song);
ts = dt:dt:T;

% In two panels, plot the song data and the neural data, 
% with the correct time axes.
figure
subplot(2,1,1)
plot(ts, song)
title('Bird song.')
ylabel('Amplitude (mV)')

subplot(2,1,2)
plot(ts, units)
title('Units firing.')
ylabel('Amplitude (\muV)')

% Play the song if you wish.
% sound( song, fs )

% You’ll notice that both the sound data and the neural trace 
% have a lot of fluctuations. We want to ignore the fast fluctuations, 
% and instead focus on whether broad bursts in neural activity precede
% song syllables (broad bursts of sound). 
% Therefore, we will compute the log power of each signal,
% so we're left with signals showing the broad changes 
% in power of each signal

% Smooth song and spike trains.
gWin = gausswin(ceil(fs * 0.026));
smoothUnits = 10 * log10(conv(units .^ 2, gWin, 'same'));
smoothSong = 10 * log10(conv(song.^2, gWin, 'same' ));

figure
subplot(2,1,1)
plot(ts, smoothSong)
ylim([-30 10])
title('Smoothed bird song power.')
ylabel('Power (dB)')

subplot(2,1,2)
plot(ts, smoothUnits)
ylim([0 7])
title('Units firing power.')
ylabel('Power (dB)')
xlabel('time (s)')

% Compute the cross-correlation between the two signals. Limit the lags to
% [-0.5, 0.5] seconds.
maxLag = 0.5;
lagSamp = round( maxLag * fs );
[songUnits, tauSamp] = xcorr(smoothUnits, smoothSong, lagSamp);

% Generate tau lags for plotting.
% tau = ( 1 : length( songUnits ) ) .* dt - 0.5 * length( songUnits ) * dt;
tau = tauSamp ./ fs;

figure
plot(tau, songUnits)
hold on
xlim([-0.5 0.5])
xlabel('\tau (s)')
ylabel('Corr. Value (a.u.)')

% Find time at which units fire in relation to song production.
[pks, locs] = findpeaks(songUnits);
[~, maxIdx] = max(pks);
maxTau = tau(locs(maxIdx));
fprintf('Unit firing precedes song production by %5.2g s.', -maxTau);

% Mark max tau.
plot( maxTau, pks( maxIdx ), 'ro',...
    'MarkerFaceColor', 'r',...
    'MarkerSize', 6 )
