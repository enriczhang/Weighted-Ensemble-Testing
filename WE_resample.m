function samples = WE_resample(n,distr)

%%%%%%%%%%%%%%%%%%%%% WE resampling function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function performs generic residual multinomial resampling

%INPUTS: 
%n = column vector giving number of samples 
%distr = matrix giving distributions to sample from

%OUTPUTS:
%samples = matrix giving the number of samples of each type

%NOTES:
%performs multinomial residual resampling on multiple distributions, 
%defined by the rows of matrix distr, with the number of samples 
%of each distribution defined by the column vector n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define residuals and numbers of samples for each distribution
ndistr = n.*distr;           %ideal, noninteger # of samples of each type
ninit = floor(ndistr);       %integer deterministic number of samples 
nresid = n - sum(ninit,2);   %residual sampling distribution

%perform multinomial resampling on the residuals
samples = mnrnd(nresid,(ndistr-ninit)./sum(ndistr-ninit,2)); %residuals
samples(isnan(samples)) = 0;  %set residuals to zero if NaN
samples = ninit+samples;      %add residuals to deterministic samples