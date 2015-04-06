function [idealMC] = MonteCarloApprox(mat, ideality)
% input is matrix with k, T, and P (in that order) for the columns
% selection of ideality definition


if ideality == 1 % definition 1 of ideality
    % ideality 1 = k*T/P
    
    idealMC = (mat(:,1).*mat(:,2))./(mat(:,3));
    
else
    print('Not a valid ideality criterion.')
end

end

