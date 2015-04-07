% main script to perform uncertainty analysis of reservoir ideality
% coded by Calvin Whealton (caw324@cornell.edu)

% import data
% may need to do this manually
% I could not find a function that worked well
% all data must be in the form of numbers
% for distributions, 1=uniform, 2=triangular,3=normal,4=lognormal
% all columns must be in the same order as the example
% [k, H, Po, Pb, Ro]
formdata = csvread('TestFormationData.csv',1,0);

% undertianty map has column 1 as 1-5 (uncertainty levels)
% columns [2,3,4,5,6] for [k, H, Po, Pb, Ro] respectively
uncermap = csvread('TestUncertaintyLevels.csv',1,0);

% setting number of replicates for Monte Carlo calculation
repsMC = 100000;

% physical constants
mu = 10; % units
Rb = 2; % units

% finding number of formations
forms = length(formdata(:,1));

% initializing matrix to hold all Monte Carlo calcutions for the formations
allFormMC = zeros(repsMC,forms);

%column index of k, H, Po, Pb, and Ro  for mean, uncertainty level, and distribution
ind_k = [2,3,4];
ind_H = [5,6,7];
ind_Po = [8,9,10];
ind_Pb = [11,12,13];
ind_Ro = [14,15,16];

% rows are variables, columns are [mean, uncertainty level, distribution]
inds_mat = [ind_k; ind_H; ind_Po; ind_Pb; ind_Ro]; % converting to matrix

for i = 1:forms % loop over formations
    
    % initializing matrix to hold the random numbers
    % columns will be k, H, Po, Pb, Ro
    rand_nums = zeros(repsMC,length(inds_mat(:,1)));
    
    % generating values for variables, looping to make code cleaner
    for j = 1:length(inds_mat(:,1))
        %                            mean value                 uncertainty level                   distribution
        rand_nums(:,j) = GenRandNums(formdata(i,inds_mat(j,1)),uncermap(formdata(i,inds_mat(j,2)),j+1),formdata(i,inds_mat(j,3)),repsMC);
    
        %hist(rand_nums(:,j))
        %pause
    
    end
    
    allFormMC(:,i) = MonteCarloApprox(rand_nums, 1, mu, Rb);
    
    %hist(allFormMC(:,i))
    %pause
    
end

% can plot histograms
%hist(allFormMC(:,1))
%hist(allFormMC(:,1))

% can calculate other statistics