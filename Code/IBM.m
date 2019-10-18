%IBM Simulate a realisation of the IBM
%   IBM(params,domain,IC,T,Nmax,seed) simulates a realisation from the
%   individual based model presented in the main document.
%
%   INPUTS:
%       params  : [m,p,gm,gp,gb,sig,mu_s]  : (1 x 7)  vector
%       domain  : [L,H]                    : (1 x 2)  vector
%       IC      : [X1,Y1; ...; XN,YN]      : (N0 x 2) matrix
%       T       : (36)                     : Solve for 0 < t < T
%       Nmax    : (5000)                   : Early stop population
%       seed    : (int)                    : rng seed
%
%   OUTPUTS:
%       NT                       = BinnyIBM( ... )
%       [XT,X18]                 = BinnyIBM( ... )
%       [NT,PCT,DT,N18,PC18,D18] = BinnyIBM( ... )
%
%   where:
%       N(18/T)     :  Population at time 18 / T
%      PC(18/T)     :  Pair correlation fcn at time 18 / T
%       D(18/T)     :  Density profile at time 18 / T
%       X(18/T)     :  [x,y] matrix of agent locations at time 18 / T
%
%   Copyright 2019 Alexander P. Browning
%   MEX Function