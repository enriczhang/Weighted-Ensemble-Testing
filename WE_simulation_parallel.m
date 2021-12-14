%%%%%%%%%%%%%%%%%%%% main WE simulation program %%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function runs a WE simulation of the 2D entropic switch problem

%NOTES:
%this runs the "entropic switch" problem on a 3-well potential energy
%depending on temperature, particles either transition directly from 
%the left well, A, to the right well, B, or they transition via state C

%at bta = 1.67, the particles tend to transition directly from A to B
%at bta = 6.67, the particles tend to transition from A to B through C

%runtime on a Surface Pro 7 is ~15 minutes for T = 10^5 select/mutate steps 
%each of length deltaT = 10 time steps, with N = 1000 particles and 10 bins

%this run time is sufficient to equilibrate at bta = 6.67 
%assuming that the integrator time step is dt = 0.01

%particle vector xs is augmented with a row labeling the last visited state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close previous mt session
delete(gcp('nocreate'))

%close all figures
close all

%initialize parallel pool
parpool(28)

%initialize WE simulation
%WE_type = input('input type: adaptive (a) or uniform (u): ','s');
%%% remove console input for remote runs

%%% List of available types
% u = uniform binning, uniform allocation (does not have difference for
%   levels only between start and end or levels "around" start)
% c = uniform binning, uniform allocation, levels between end and start
% h = heirarchal binning, wv allocation
% a = kmeans binning, wv allocation
% m = adaptive binning over muv to get uniform binning over h
% g = static binning over muv to get uniform binning over h
% s = same as g but with different creation of x,y pairs
% p = Aristoff suggestion on binning in h
% t = kmeans h density (static polling in init)
% v = kmeans on muv density (static polling in init) w/ reverting back to h

for new_type=['u','a','g','t','v','m']
for max_bins=[6,3]
for init_pop=[50]

rerun = 1; %parameter to reload init files
WE_type = new_type;
%max_bins = 6;
% init_pop = 500; %can change for new populations
if WE_type == 'f'
    Bin_Floor =  5;
else
    Bin_Floor = 1;
end

loadname = sprintf('initialization_data_%s_%d_%d.mat',WE_type,Bin_Floor,init_pop);
if rerun == 0
    try
        load(loadname)
    catch
        WE_live_init(WE_type,Bin_Floor,init_pop,max_bins);
        load(loadname)
    end
else
    WE_live_init(WE_type,Bin_Floor,init_pop,max_bins);
    load(loadname)
end

trials = 500;
T = 10^5;	%%% ensure convergence
data = zeros(trials,1);
type = WE_type;

%loop over independent trials
parfor trial=1:trials %parfor
    xs = [0.1;0.5] + [normrnd(0,0.02,[2 N])];
    ws = ones(1,N)/N;
    obs_avg = 0;
%main WE simulation loop
for t=1:T
    
    %if(t/1000) == floor(t/1000)
    %    t/T
    %end

    %plot WE particles
    %WE_plot      %comment this out for long simulations

    %define WE bins, u, and particle allocation, Nu
    [wus,copies] = WE_parameters_(xs,ws);
    
    %perform WE selection or resampling step
    [xs,ws] = WE_selection_(xs,wus,copies);    
    
    %perform WE evolution or mutation step
    [xs,obs_avg] = WE_evolution_(xs,ws,obs_avg);   

end

data(trial) = obs_avg/(T*deltaT);

end

%display WE observable average
%disp('WE observable average = ...')
WE_average = mean(data);

%disp('WE observable error = ...')
WE_std = std(data)/sqrt(trials);

%disp('exact observable average = ...')
exact_average;

savename = sprintf('RERUN_WE_data_%s_%d_%d_%dbins_%dtrials.mat',WE_type,Bin_Floor,init_pop,max_bins,trials);
save(savename)
end %max_bins
end %population loop
end %WE_type loop
% end order is moved, type inside pop inside bins
