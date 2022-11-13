function V = generatevoltage( p, T, Vreset, Vthresh, V0 )
% takes as input the probability of going up, 
% p, the number of time steps to simulate, T, the reset voltage, Vreset, 
% the spike threshold, Vthres, and the initial voltage, V0.

% create the vector V = zeros(1,T); that will hold voltage values at each
% step of the process 
V = zeros( 1, T );
V( 1 ) = V0;

for t = 1 : T - 1
    if V(t) == Vthresh
        V(t + 1) = Vreset; 
    else
        dV = (( rand(1)<p) - 0.5) * 2;
        V(t+1) = V(t) + dV;  
    end
    if V(t + 1) > Vthresh
        V(t + 1) = Vthresh;   
    end
end