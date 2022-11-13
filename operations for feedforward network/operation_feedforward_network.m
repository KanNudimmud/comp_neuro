%% Exercise 1: Matrix Operations for a Feedforward Network
% We will construct a 2 layer feedforward linear network, 
% and use matrix operations to calculate its
% outputs, given its inputs and weights.
% We’ll call the output neurons y1, …, yM and input neurons
% x1, …, xN. Wij is the connection strength (weight) 
% onto neuron yi from neuron xj. We refer to W as the weight matrix
close all, clear all, clc

M = 10;
N = 50;

% Generate a weight matrix. 
% Assume the weights are random and uniform between -1 and 1
W = rand( M, N );
W = 2.*( W - 0.5 );

% Generate a 50-dimensional pattern of inputs consisting of Gaussian entries
x = randn( 1, 50 )';

% Calculate the network output
y = W * x;

figure
subplot(3,1,1 )
plot(x)
box off
xlabel('sample number')
ylabel('input values (a.u.)')

subplot(3,1,2)
pcolor(W)
title('Weights matrix')

colorbar
subplot(3,1,3)
plot(y)
box off
xlabel('sample number')
ylabel('output values (a.u.)')

%% Exercise 2: Logical operations, for-loops, and plotting (random walk)
% We will construct a biased random walk with a reflecting barrier and reset.
% You can think of this as a simple model of the voltage of a neuron. 
% The voltage increases, with some fluctuations,until it reaches a threshold,
% at which point it spikes then resets to a resting value, and starts the
% process again.

% Specifically, the voltage of the cell is given by: V(t+1) = V(t) +dV(t), 
% where dV(t) is 1 with probability p, and -1 with probability 1 - p. 
% We want the system to reset once it hits a maximum ceiling, 
% so we will add the condition: V(t+1) = Vreset if Y(t) > Vthres.

clear all, close all, clc

% simulate a neuron with initial voltage of -65mV, threshold of -45mV, 
% and reset voltage of -70mV. 
% Choose p so that your neuron has an average firing rate of approximately 10Hz, 
% assuming each time step corresponds to 1ms. 
% Calculate the voltage values for 1 second

T = 10000; % measured in ms
Vreset = -70; % measured in mV
Vthresh = -45; % measured in mV
V0 = -65; % measured in mV

% find p such that the firing rate is approximately 10Hz
firingRate = 0;
p = 0;
while 1
    if abs( firingRate - 10 ) < 0.1
        break
    end
    
    V = generatevoltage( p, T, Vreset, Vthresh, V0 );
    firingRate = sum( V == Vthresh ) / ( T / 1000 );
    p = p + 0.0001;
    
    if p > 1
        disp( 'did not find a suitable p' )
        break  
    end 
end

ts = 1 : T; % time vector in ms

figure
plot(ts,V)
box off
xlabel('time (ms)')
ylabel('voltage (mV)')
xlim([0 1000])
ylim([-80 -40])
% title(['probability that dV = 1 is ', num2str(p)])