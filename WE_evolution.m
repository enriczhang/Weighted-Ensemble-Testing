function [xs,obs_avg] = WE_evolution(xs,ws,obs_avg,dt,deltaT,bta,F,N)

%%%%%%%%%%%%%%%%%%% WE evolution function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function evolves the WE particles 
%for deltaT time steps each of size dt

%INPUTS:
%xs = dxN particle matrix, dt = time step 
%deltaT = number of time steps per evolution step
%bta = inverse temperature
%F = force, d = dimension, N = number of particles

%OUTPUTS:
%xs = dxN updated particle matrix

%NOTES:
%first d rows of xs are the particle positions in d dimensions
%where each column of xs represents a single particle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%evolve the WE particle vector for DeltaT time steps
for t=1:deltaT
    
    %evolve WE particles for 1 integrator time step of size dt
    xs = xs + F(xs)*dt + sqrt(2*dt/bta)*normrnd(0,1,[2,N]);
    
    %recycle particles that reach the sink
    sink = (xs(2,:)>0.9 & xs(1,:)>0.5 & xs(1,:)<0.6);
    xs = [0.1;0.5].*sink + xs.*(1-sink);
    
    %restrict particles to unit square
    tol = 10^(-5);
    xs = min(max(xs,tol),1-tol);
    
    %update observable average
    obs_avg = obs_avg + sum(ws(sink));
    
end
