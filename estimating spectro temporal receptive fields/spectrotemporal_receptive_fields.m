%% Exercise 7: Estimation of spectro-temporal receptive fields with spike-triggered averaging
clear all , close all, clc

% we will demonstrate how to estimate the spectro-temporal receptive field
% (STRF) using spike-triggered averaging (STA). 
% First, we will generate a spike train from a model neuron 
% with a spectro-temporal kernel that we will provide. 
% Whenever the stimulus is sufficiently correlated with the kernel, 
% the neuron will fire. Next, we will compute the STA 
% (the average stimulus that precedes a spike). 
% We will compare the STA with the kernel used to generate the spike train,
% note that they look similar, and comment on the limitations that prevent
% us from perfectly estimating the STRF.

% We begin by constructing a stimulus to present to the neuron. 
% We know the neuron is selective to time-varying sound, 
% but we want to know which sound best excites the neuron. 
% Therefore, we want a stimulus that will 
% try many different sequences of sounds with equal probability. 
% The MATLAB file generatestimulus.m generates such
% a stimulus in the form of a 2D matrix where the rows 
% represent 50 logarithmically spaced tone frequencies 
% and the columns represent time bins of width 1ms. 
% Entry (i, j) of the matrix is the amplitude of tone i at time step j. 
% The stimulus is constructed so that tones turn on and last for 30 ms, 
% with a slight ramping at the onset and offset 
% to make transitions sound less abrupt. 
% Each tone has a 10% probability of turning on in a given 30 ms window.

stimulus = generatestimulus;

freq = logspace(2, 4, 50);% 50 log-spaced freqs 100 - 10,000 Hz
dt = 0.01; % sampling interval in seconds

% Listen to the first 3 seconds of the stimulus.
% playstim( stimulus( :, 1 : ( 3 / dt ) ), freq, dt);

% Next we will construct an STRF for our model neuron. 
% The file generatekernel.m returns an array of 100 time bins by 50 frequency bins, 
% containing the STRF. The STRF is the sum of bivariate Gaussian distributions.

kernel = generatekernel;
% playstim( Kernel, freq, 0.001 )

% Now we are ready to simulate how the neuron responds to the stimulus. 
% To do this, you will slide the kernel across the stimulus, 
% and at each time bin t calculate the integral (sum) 
% of the element-wise product between the kernel and the stimulus
% that occurred between time (t - 100) and time (t – 1). 
% This is our estimate of how strongly the neuron is excited 
% by the stimulus at each time. 
% Our simulated neuron should spike whenever this excitatory drive 
% exceeds some threshold.

estResponse = zeros(1, size(stimulus, 2));
for i = size(kernel, 2) + 1:size(stimulus, 2)
    response = kernel .* stimulus(:, i - size(kernel, 2):i - 1);
    estResponse(i) = sum(response(:));
end

% calculate the excitatory drive to the neuron at each time in the stimulus.
% You can start 100ms into the stimulus (the STRF is 100 ms long). 
% code should record a spike whenever excitatory drive 
% to the neuron exceeds a threshold. 
% Choose a threshold such that the neuron spikes 
% approximately 10,000 times over the entire stimulus.

sortedEstResp = flip(sort(estResponse));
thresh = sortedEstResp(1e4);
spikes =  estResponse >= thresh;
spikes = double(spikes);

% Plot variables.
figure
subplot(3,1,1)
imagesc(stimulus(:, 1:(5 / dt)))

subplot(3,1,2)
plot(estResponse(1 : (5 / dt)))

subplot(3,1,3)
plot(spikes(1 : (5 / dt)))

% Find STA
%  calculate the average 150 ms of stimulus that precedes each spike
sta = zeros(50, 150);
for i = find(spikes == 1)
    sta = sta + (1 / sum(spikes)) .* stimulus(:, i-150:i-1);
    
end
    
figure
surf(1:150, logspace(2, 4, 50), sta, 'edgecolor', 'none');
axis tight;
view(0, 90);
shg
set(gca, 'yscale', 'log')
xlabel('time (ms)')
ylabel('Frequency (Hz)')
title('STA')
colorbar;

figure
imagesc(kernel)
title('Kernel (STRF)')

% playstim( sta, freq, 0.001 )