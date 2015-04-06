% main script to perform uncertainty analysis of reservoir ideality
% coded by Calvin Whealton (caw324@cornell.edu)

% import data
% my need to do this manually
% I could not find a function that worked well
% all data must be in the form of numbers
% for distributions, 1=uniform, 2=triangular,3=normal,4=lognormal
% all columns must be in the same order as the example
formdata = csvread('TestFormationData.csv',1,0);

% undertianty map has column 1 as 1-5 (uncertainty levels)
% columns 2,3,4 for k, T, and P, respectively
uncermap = csvread('TestUncertaintyLevels.csv',1,0);

% setting number of replicates for Monte Carlo calculation
repsMC = 100000;

% finding number of formations
forms = length(formdata(:,1));

% initializing matrix to hold all Monte Carlo calcutions for the formations
allFormMC = zeros(repsMC,forms);

%column index of k, T, and P  for mean, uncertainty level, and distribution
ind_k = [2,3,4];
ind_T = [5,6,7];
ind_P = [8,9,10];

for i = 1:forms % loop over formations
    
    % initializing matrix to hold the random numbers
    % columns will be k, T and P (in that order)
    rand_nums = zeros(repsMC,3);
    
    % generating values for variables
    % for k, column 1
    rand_nums(:,1) = GenRandNums(formdata(i,ind_k(1)),uncermap(formdata(i,ind_k(2)),2),formdata(i,ind_k(3)),repsMC);
    
    % for T, column 2
    rand_nums(:,2) = GenRandNums(formdata(i,ind_T(1)),uncermap(formdata(i,ind_T(2)),3),formdata(i,ind_T(3)),repsMC);
    
    % for P, column 3
    rand_nums(:,3) = GenRandNums(formdata(i,ind_P(1)),uncermap(formdata(i,ind_P(2)),4),formdata(i,ind_P(3)),repsMC);    
    
    % calculating the MC ideality
    allFormMC(:,i) = MonteCarloApprox(rand_nums,1);
    
end

% can plot histograms
hist(allFormMC(:,2))

% can calculate statistics