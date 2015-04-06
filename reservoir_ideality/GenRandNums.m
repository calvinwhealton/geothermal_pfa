function [values] = GenRandNums(mean, unc, dist, reps)
% function to generate random numbers
% inputs are: 
% mean
% uncertainty (%)
% distribution
% number of replicates

% initializing random number generator
rng(10);

if dist == 1 % uniform distribution
    
    lb= (1 - unc/100)*mean; % lower bound of the uniform
    ub = (1 + unc/100)*mean; % upper bound of the uniform
    
    % random values
    values = rand(reps,1)*(ub- lb) + lb; 
    
elseif dist == 2 % triangular distribution
    
    lb = (1 - unc/100)*mean; % lower bound of the triangular
    ub = (1 + unc/100)*mean; % upper bound of the triangular
    
    % random values
    rv = rand(reps,1);
    % converting from cdf to quantile, first term ensures the correct part
    % of the cdf is used to make the second part valid
    values1 = (rv < 0.5).*(lb + sqrt(rv*(ub-lb)*(mean-lb)));
    values2 = (rv > 0.5).*(ub - sqrt((1-rv)*(ub-lb)*(ub-mean)));
    values  = values1 + values2; 
    
elseif dist == 3 % normal distribution
    
    % assuming that uncertainty maps to 95% confidence interval
    std_dev = (unc/100)*mean/norminv(0.975);
    
    % generating random values
    values = normrnd(mean,std_dev,[reps,1]);
    
elseif dist == 4 % lognormal distribution
    % for lognormal distribution
    % real-space mean is exp(mu + sigma^2/2)
    % real-space variance is (exp(sigma^2)-1)*(exp(2mu + sigma^2))
    % need mu and sigma to generate values
    
    % assume uncertainty is standard deviation as a percentage of the mean
    
    % squared real-space coefficient of variation
    cv2 = (unc/100)^2; 
    
    % solve for sigma^2
    sigma2 = log(1 + cv2);
    
    % solve for mu
    mu = log(mean) - sigma2/2;
    
    % renerating random values
    values = exp(normrnd(mu,sqrt(sigma2),[reps,1]));
    
else % case for not a vaild distribution type
    print('No valid distribution selected.')
end

end

