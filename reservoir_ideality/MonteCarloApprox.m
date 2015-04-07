function [idealMC] = MonteCarloApprox(mat, ideality, mu, Rb)
% input is matrix with k, T, and P (in that order) for the columns
% selection of ideality definition


if ideality == 1 % definition 1 of ideality
    
    % ideality 1 = (2*pi/mu)k*H*(Po-Pb)/ln(Rb/Ro)
    
    % calculating the pressure difference Po - Pb
    % Po < Pb, which should almost always be the case
    % using absolute values makes distribution foled normal, not normal, but it should be
    % very clsoe because of the difference in Po and Pb is presumably large
    % (hydrostatic versus lithostatic)
    Pdiff = -1*abs(mat(:,3) - mat(:,4));
    
    % using log(a/b) = log(a) - log(b)
    % log(x) in MATLAB is natrual logarithm, log10(x) is base 10 logarithm
    %         const        k         H          Po-Pb    log(Rb)    log(Ro)
    idealMC = (2*pi/mu)*((mat(:,1).*mat(:,2)).*Pdiff)./(log(Rb) - log(mat(:,5)));
  
else
    print('Not a valid ideality criterion.')
end

end

