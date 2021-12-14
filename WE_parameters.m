function [wus,copies] = ... 
    WE_parameters(xs,ws,N,h_,v_,RMSD_,levels,type,tol,bin_floor,num_bins)

%%%%%%%%%%%%%%%%%%%%% WE parameters function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function computes the parameters needed for the WE selection step

%INPUTS:
%xs = (d+1)xN particle matrix, ws = 1xN weight vector
%levels = level sets of reaction coordinate used to define fixed WE bins 
%tol = weight tolerance -- sets lower bound on weights

%OUTPUTS:
%wus = 1xM vector of uniform weights in bins after the WE selection step
%copies = 1xN vector defining the number of copies of each particle

%NOTES:
%the WE bins are chosen based on level sets of the reaction coordinate xi
%and the WE allocation is uniform among bins whose weight is > tol
%the parameter tol prevents splitting of particles with very small weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define bin occupancy matrix u
%(i,j)th entry of u equals 1 if bin i contains particle j, and 0 else
if type == 'u' || type == 'c'
    u = ([-inf;levels]<=RMSD_(xs))&(RMSD_(xs)<[levels;inf]);
    u = u(any(u,2),:);
elseif type == 'h'
    bins = clusterdata(h_(xs)',num_bins);
    u = (bins == linspace(1,num_bins,num_bins))';
elseif type == 'm'
    % Compare this wv to muv
    % maybe in both
    % set up initialization to handle this
    h_data = h_(xs);
    [sorted_h, sort_h_ind] = sort(h_data); %~ ignores the output like _ in python
    %sort everything by h to make it easier on indexing
    xs = xs(:,sort_h_ind);
    ws = ws(sort_h_ind);
    val_sort = v_(xs); %% called twice in this implementation
    sort_wv = ws.*val_sort;
    cs_wv = cumsum(sort_wv);
    
    cs_muv_bins = linspace(cs_muv(1),cs_muv(end),num_bins+1);
    cs_muv_bins = cs_muv_bins(2:end-1);
    
    %convert these to bins in h
    h_bins = zeros(1,num_bins-1);
    for i=1:length(h_bins)
        h_bins(i) = sorted_h(find(cs_muv > cs_muv_bins(i),1));
    end
    levels = h_bins';

    u = ([-inf;new_lvls]<=cs_wv)&(cs_wv<[new_lvls;inf]);
    u = u(any(u,2),:); % should be unnecessary by design
elseif type == 'g' || type == 's'
    opt_all = h_(xs);
    u = ([-inf;levels]<=opt_all)&(opt_all<[levels;inf]);
    u = u(any(u,2),:); % should be unnecessary by design
elseif type == 'p'
    H = h_(xs);
    H = H + h_min;
    H = H / h_max;
    H = H.^(1/root);
    u = ([-inf;levels]<=H)&(H<[levels;inf]);
    u = u(any(u,2),:);
elseif type ==  't'
    u = ([-inf;levels]<=h_(xs))&(h_(xs)<[levels;inf]);
    u = u(any(u,2),:);
elseif type == 'v'
    u = ([-inf;levels]<=h_(xs))&(h_(xs)<[levels;inf]);
    u = u(any(u,2),:);
else
    %levels = linspace(-0.01,1,numbins)';    
    %u = ([-inf;levels]<=h_(xs))&(h_(xs)<[levels;inf]);
    %u = u(any(u,2),:);
    
    bins = kmeans(h_(xs)',num_bins);
    u = (bins == linspace(1,num_bins,num_bins))';
end

%define bin weight matrix Pu
%(i,j)th entry of Pu equals weight of particle j in bin i
Pu = ws.*u;   

%defines total bin weight vector wus
wus = sum(Pu,2);

%define WE allocation distribution, distr
if type == 'u'
    %here, distr is the uniform distribution in bins with total weight > tol
    distr = ((wus > tol)/sum(wus > tol))';  
else
    %here, distr is the optimal choice, weights*v
    vals = v_(xs).*Pu;
    distr = sum(vals,2);
    distr = distr'/sum(distr);
end

%define integer particle allocation, Nu, approximately following distr
% Nu = 1 + WE_resample(N-size(u,1),distr)';
Nu = bin_floor + WE_resample(N-bin_floor*size(u,1),distr)';

%normalize the bin weight matrix Pu
Pu = Pu./wus;

%compute the uniform weights in each bin, wus, after selection
wus = wus./Nu;

%compute the matrix of the number of copies of each particle in each bin
%(i,j)th entry of copies = number of copies of particle j in bin i
copies = WE_resample(Nu,Pu);

