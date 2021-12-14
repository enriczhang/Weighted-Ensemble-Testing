function [blank] = WE_live_init(IType, BFloor, IPop, MBins)

%%% List of available types
% u = uniform binning, uniform allocation (does not have difference for
%   levels only between start and end or levels "around" start
% h = heirarchal binning, wv allocation
% a = kmeans binning, wv allocation
% f = 5_min with kmeans
% m = adaptive binning over muv to get uniform binning over h
% g = static binning over muv to get uniform binning over h
% s = same as g but with different creation of x,y pairs

bta = 15;      %6.67;       %bta = inverse temperature
dt = 0.001;      %dt = time step
deltaT = 10;      %deltaT = number of time steps between resamplings
T = 10^4;        %T = total number of time steps
N = IPop;         %N = number of particles
d = 2;           %d = dimension of space of particles
tol = 10^(-20);  %tol = tolerance for weights -- sets lower bound
obs_avg = 0;     %obs_avg = WE observable average
type = IType;      %u = uniform binning+alloc, a = h-adaptive binning+v-alloc
num_bins = MBins;	 %can become input variable
if type == 'c'
    levels=(linspace(0,1,6)*sqrt((0.1-.5)^2+(0.5-1)^2))'; %levels for uniform binning between
else
    levels = [0; 1/4; 2/4; 3/4; 1; 5/4];   %levels for uniform binning
end

bin_floor = BFloor %sets the minimum allocation per occupied bin

%define WE potential, V, and force, F = -\nabla V, RMSD reaction coordinate, 
%discrepancy h, and variance function v

%define potential, V, and force, V
a = 50.5;
b = 49.5;
c = 10^5;
d = 51;      %fix this: it is used for dimension below.
e = 49;
del = 0.2;
V = @(x,y) exp(-(a*(x-0.25).^2+a*(y-0.75).^2+2*b*(x-0.25).*(y-0.75))) ... 
                + exp(-c*(x.^2.*(1-x).^2.*y.^2.*(1-y).^2)) ... 
                + 0.5*exp(-(d*x.^2+d*y.^2-2*e*x.*y)); ... 

F = @(x) -[- exp(-c*x(1,:).^2.*x(2,:).^2.*(x(1,:) - 1).^2.*(x(2,:) - 1).^2).*(2*c*x(1,:).*x(2,:).^2.*(x(1,:) - 1).^2.*(x(2,:) - 1).^2 + c*x(1,:).^2.*x(2,:).^2.*(2*x(1,:) - 2).*(x(2,:) - 1).^2) - exp(- a*(x(1,:) - 1/4).^2 - a*(x(2,:) - 3/4).^2 - 2*b*(x(1,:) - 1/4).*(x(2,:) - 3/4)).*(2*b*(x(2,:) - 3/4) + a*(2*x(1,:) - 1/2)) - (exp(- d*x(1,:).^2 + 2*e*x(1,:).*x(2,:) - d*x(2,:).^2).*(2*d*x(1,:) - 2*e*x(2,:)))/2; ...
    - exp(-c*x(1,:).^2.*x(2,:).^2.*(x(1,:) - 1).^2.*(x(2,:) - 1).^2).*(2*c*x(1,:).^2.*x(2,:).*(x(1,:) - 1).^2.*(x(2,:) - 1).^2 + c*x(1,:).^2.*x(2,:).^2.*(2*x(2,:) - 2).*(x(1,:) - 1).^2) - exp(- a*(x(1,:) - 1/4).^2 - a*(x(2,:) - 3/4).^2 - 2*b*(x(1,:) - 1/4).*(x(2,:) - 3/4)).*(2*b*(x(1,:) - 1/4) + a*(2*x(2,:) - 3/2)) - (exp(- d*x(1,:).^2 + 2*e*x(1,:).*x(2,:) - d*x(2,:).^2).*(2*d*x(2,:) - 2*e*x(1,:)))/2];

%construct coarse matrix, P
sig = sqrt(2*dt/bta);
n = 100;
P = zeros(n,n,n,n);
V_ = zeros(n,n);
f = zeros(n,n);
RMSD = zeros(n,n);
for i=1:n
    i/n
    for j=1:n
        V_(i,j) = V((i-0.5)/n,(j-0.5)/n);
        RMSD(i,j) = norm([0.55, 0.95]-[(i-0.5)/n, (j-0.5)/n]);
        if i > n/2 && i <= n/2 + 10 && j>n-10
            P(i,j,10,floor(n/2)) = 1;
            f(i,j) = 1;
        else
            f(i,j) = 0;
            for k=1:n
                for l=1:n
                    x = (i-0.5)/n;
                    y = (j-0.5)/n;
                    z = (k-0.5)/n;
                    w = (l-0.5)/n;
                    P(i,j,k,l) = ...
                    exp(-norm([z; w] - ([x; y] + F([x;y])*dt))^2/(2*sig^2));
                end
            end
        end
    end
end

%solve for mu, h and v, in vector form
P = reshape(P,n^2,n^2);
P = P./sum(P,2);
f = reshape(f,n^2,1);
[mu,lam] = eigs(P',1);
mu = mu/sum(mu);
h = linsolve(eye(n^2)-P-mu*mu',f-(mu'*f)*ones(n^2,1));
v = sqrt(P*h.^2-(P*h).^2);
v(1) = tol;
v(n^2) = tol;
h(1) = min(h);
h(n^2) = min(h);
RMSD = reshape(RMSD,n^2,1);

%compute ``exact'' observable average
exact_average = 0;
mu2 = reshape(mu,n,n);
for i=n/2+1:n/2+10
    for j=n-9:n
        exact_average = exact_average + mu2(i,j);
    end
end

%plot optimal allocation
% v2 = reshape(v,n,n);
% surf(mu2.*v2)

%define WE functions mu, h, and v in function form
grd = @(xs) max(min(n*floor(n*xs(2,:)) + ceil(n*xs(1,:)),n^2),1);
h_ = @(xs) h(grd(xs))';
v_ = @(xs) v(grd(xs))';
mu_ = @(xs) mu(grd(xs))';
RMSD_ = @(xs) RMSD(grd(xs))';

%fix global variables for WE functions
%WE_parameters_ = @(xs,ws) ...
%    WE_parameters(xs,ws,N,h_,v_,RMSD_,levels,type,tol,bin_floor,num_bins);
%WE_selection_ = @(xs,wus,copies) WE_selection(xs,wus,copies);
%WE_evolution_ = @(xs,ws,obs_avg) ...
%    WE_evolution(xs,ws,obs_avg,dt,deltaT,bta,F,N);

%define initial particles in A: last in A labeled 0, last in B labeled 1
xs = [0.1;0.5] + [normrnd(0,0.02,[2 N])];
%xs = rand([2 100*N]);

%define initial uniform weights
ws = ones(1,N)/N;
%ws = mu_(xs);
%ws = ws/sum(ws);

if type == 'g'
    xgrd = linspace(0.01,0.99,100);
    [X,Y] = meshgrid(xgrd,xgrd);
    Z=cat(2,X',Y');
    W=reshape(Z,2,[]);
    H=h_(W);
    [sorted_h, sort_h_ind] = sort(h); %~ ignores the output like _ in python
    %sort everything by h to make it easier on indexing
    W = W(:,sort_h_ind);
    val_sort = v_(W); %% called twice in this implementation
    mu_sort = mu_(W);
    sort_muv = mu_sort.*val_sort;
    cs_muv = cumsum(sort_muv);
    
    %get bins equal in mu*v first
    cs_muv_bins = linspace(cs_muv(1),cs_muv(end),num_bins+1);
    cs_muv_bins = cs_muv_bins(2:end-1);
    
    %convert these to bins in h
    h_bins = zeros(1,num_bins-1);
    for i=1:length(h_bins)
        h_bins(i) = sorted_h(find(cs_muv > cs_muv_bins(i),1));
    end
    levels = h_bins';
elseif type == 's'
    xgrd = linspace(0.01,0.99,100);
    [X,Y] = meshgrid(xgrd,xgrd);
    W=combvec(xgrd,xgrd);
    H=h_(W);
    [sorted_h, sort_h_ind] = sort(h); %~ ignores the output like _ in python
    %sort everything by h to make it easier on indexing
    W = W(:,sort_h_ind);
    val_sort = v_(W); %% called twice in this implementation
    mu_sort = mu_(W);
    sort_muv = mu_sort.*val_sort;
    cs_muv = cumsum(sort_muv);
    
    %get bins equal in mu*v first
    cs_muv_bins = linspace(cs_muv(1),cs_muv(end),num_bins+1);
    cs_muv_bins = cs_muv_bins(2:end-1);
    
    %convert these to bins in h
    h_bins = zeros(1,num_bins-1);
    for i=1:length(h_bins)
        h_bins(i) = sorted_h(find(cs_muv > cs_muv_bins(i),1));
    end
    levels = h_bins';
elseif type == 'p'
    [sorted_h, sort_h_ind] = sort(h);
    h_min = sorted_h(1);
    h_max = sorted_h(end);
    root = 4;
    sorted_h = sorted_h - sorted_h(1); 
    sorted_h = sorted_h / sorted_h(end);
    sorted_h = (sorted_h).^(1/4);
    
    levels = linspace(sorted_h(1), sorted_h(end), num_bins+1);
    levels = (levels(2:end-1))';
elseif type == 't'
    [sorted_h, sort_h_ind] = sort(h);
    bins = kmeans(sorted_h,num_bins);
    
    [sb,si]=sort(bins);
    sorted_h = sorted_h';
    resorted_h = sorted_h(:,si);
    levels = zeros(1,num_bins);
    for i=1:num_bins
        levels(i) = resorted_h(find(sb == i,1));
    end
    levels = sort(levels);
    levels = levels(2:end);
    levels = levels';
elseif type == 'v'
    xgrd = linspace(0.01,0.99,100);
    [X,Y] = meshgrid(xgrd,xgrd);
    W=combvec(xgrd,xgrd);
    H=h_(W);
    [sorted_h, sort_h_ind] = sort(h); %~ ignores the output like _ in python
    %sort everything by h to make it easier on indexing
    W = W(:,sort_h_ind);
    val_sort = v_(W); %% called twice in this implementation
    mu_sort = mu_(W);
    sort_muv = mu_sort.*val_sort;
    cs_muv = cumsum(sort_muv);
    
%     bins = kmeans(sort_muv',num_bins);
    bins = kmeans(cs_muv',num_bins);
    
    [sb,si]=sort(bins);
    sorted_h = sorted_h';
    resorted_h = sorted_h(:,si);
    levels = zeros(1,num_bins);
    for i=1:num_bins
        levels(i) = resorted_h(find(sb == i,1));
    end
    levels = sort(levels);
    levels = levels(2:end);
    levels = levels';
end

%%fix global variables for WE functions
WE_parameters_ = @(xs,ws) ...
    WE_parameters(xs,ws,N,h_,v_,RMSD_,levels,type,tol,bin_floor,num_bins);
WE_selection_ = @(xs,wus,copies) WE_selection(xs,wus,copies);
WE_evolution_ = @(xs,ws,obs_avg) ...
    WE_evolution(xs,ws,obs_avg,dt,deltaT,bta,F,N);

fname = sprintf('initialization_data_%s_%d_%d.mat', type,bin_floor,N);
save(fname)
