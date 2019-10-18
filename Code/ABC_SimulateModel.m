function [PriorSamples,Seeds,N,PC,D] = ABC_SimulateModel(Nsamples,PriorSeed,SeedOffset)
%ABC_SimulateModel Simulate IBM model (must be in same directory)
%with samples from prior and save in file with name specified.

% REQUIRED FILES
%   ../Data folder
%   IBM.mex (OS dependent)
%   ABC_Prior.m

    repnames = {'8000_E2', '8000_G2', '10000_H1',  ...
                '10000_H2','8000_F2', '10000_A2', ...
                '12000_H3','12000_F3','12000_D2'};

    Nmax     = 4000;
    domain   = [1440,1900];
    T        = 36;
    
    % Get Initial Conditions
    IC      = cell(1,9);
    for r = 1:9

        X           = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_12h.csv']);
        IC{r}       = X;

    end
    
    % Sample parameters + Cpp seeds
    rng(PriorSeed);
    PriorSamples = ABC_Prior(Nsamples);
	rng('Shuffle');
    
    % Seeds
    Seeds    = reshape((1:(Nsamples * 9)) + SeedOffset,Nsamples,9);

    % Initialise outputs
    N        = zeros(Nsamples,9,2);
    PC       = zeros(Nsamples,9,2,20);
    D        = zeros(Nsamples,9,2,80);

    % Loop through samples (in parallel)
    parfor i = 1:Nsamples

        % Get parameters
        params = PriorSamples(i,:);

        % Loop through densities
        for r = 1:9

            [NT,PCT,DT,N18,PC18,D18] = IBM([params,24,params(end)],domain,IC{r},T,Nmax,Seeds(i,r)); %#ok<PFBNS>

            % Save results
            N(i,r,:)    = [N18,NT];
            PC(i,r,:,:)   = reshape([PC18,PCT]',1,1,2,20);
            D(i,r,:,:)    = reshape([D18,DT]',1,1,2,80);

        end

    end
            
end

