function [Weights,Eps,Models,Seeds,Samples] = ABC_SMC(Nsamples,sigma,thresh_seq,SeedOffset)
%ABC_SMC Perform ABC SMC with IBM model (must be in same 
% directory) and save results in file with name specified.

% REQUIRED FILES
%   ../Data folder
%   IBM.mex (OS dependent)
%   ABC_Prior.m
%   Data_PairCorrelationFcn.m

%% SETUP

    % Suppress warnings
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

    % Models       m  p  gm gp gb
    ModelNZ     = [1, 1, 1, 1, 1;
                   1, 1, 0, 1, 1;
                   1, 1, 1, 1, 0;
                   1, 1, 0, 1, 0;
                   1, 1, 0, 0, 0];
    ModelNZ     = logical(ModelNZ);
               
    extinct     = logical([0,0,0,0,0]);
    
    n_thresh    = length(thresh_seq);
    ss_weights  = [1,1,1];
       
    % SMC kernel
    Ksd         = [0.25,0.006,0.3,0.005,25];
    QV          = diag(Ksd.^2);
    Ksample     = @(ThStar,QVnz)       mvnrnd(ThStar,QV.*QVnz);
    Kdensity    = @(Theta,ThStar,QVnz) mvnpdf(Theta(:,QVnz),ThStar(QVnz),QV(QVnz,QVnz));
    
    % Simulation settings
    Nmax        = 4000;
    domain      = [1440,1900];
    T           = 36;
       
    % Storage
    Samples     = cell(1,n_thresh);
    Weights     = zeros(Nsamples,n_thresh);
    Model       = zeros(Nsamples,n_thresh);
    Eps         = zeros(Nsamples,n_thresh);
    Seeds       = zeros(Nsamples,n_thresh);
    
%% EXPERIMENTAL DATA

    repnames = {'8000_E2', '8000_G2', '10000_H1',  ...
                '10000_H2','8000_F2', '10000_A2', ...
                '12000_H3','12000_F3','12000_D2'};
            
    % Get Initial Conditions and Summary Statistics
    IC          = cell(1,9);
    Nexp18      = zeros(1,9);   % Population (t = 18 h)
    NexpT       = zeros(1,9);   % Population (t = 30 h)
    Dexp18      = zeros(80,9);  % Density Profile (t = 18 h)
    DexpT       = zeros(80,9);  % Density Profile (t = 30 h)
    PCexp18     = zeros(20,9);  % Pair Correlation (t = 18 h)
    PCexpT      = zeros(20,9);  % Pair Correlation (t = 30 h)
    
    for r = 1:9

        X           = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_12h.csv']);
        IC{r}       = X;
        
        X               = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_30h.csv']);
        Nexp18(r)       = size(X,1);
        Dexp18(:,r)     = histcounts(X(:,2),0:1900/80:1900)';
        PCexp18(:,r)    = Data_PairCorrelationFcn(X,400,5,100)';
        
        X               = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_48h.csv']);
        NexpT(r)        = size(X,1);
        DexpT(:,r)      = histcounts(X(:,2),0:1900/80:1900)';
        PCexpT(:,r)     = Data_PairCorrelationFcn(X,400,5,100)';

    end
    
%% INITIALISE
    
    WeightLast  = ones(Nsamples,1);
    ThetaCur    = zeros(Nsamples,5);
    WeightCur   = ones(Nsamples,1);
    ModelCur    = ones(Nsamples,1);
    epsTheta    = zeros(Nsamples,1);
    SeedCur     = zeros(Nsamples,1);
    
    tic;
    % LOOP THROUGH PARTICLES
    thresh_cur  = thresh_seq(1);
    parfor j = 1:Nsamples
        
        % initialise seed, i = 1
        seed = SeedOffset + (j-1) * 1e3;
               
        accepted = false;
        while ~accepted
            
            % Sample model
            Mk0       = SampleModel(extinct);
            
            % Propose (adjust for model)
            Theta0   = ABC_Prior(1) .* ModelNZ(Mk0,:);

            % Initialise Discrepancy
            kappa0   = 0;
            
            % Simulate model (look through data sets)
            for r = 1:9
               
                % Determine Nmax
                Nmax = DetNmax(NexpT(r),thresh_cur - kappa0,ss_weights(1));
                
                % Run Model
                [NT,PCT,DT,N18,PC18,D18] = BinnyIBM([Theta0,sigma,24,sigma],[1440,1900],IC{r},36,Nmax,seed);
                seed = seed + 1;
                
                if NT > Nmax
                    kappa0 = kappa0 + 1e9;
                    break;
                end
            
                % Increase eps
                kappa0  = kappa0 + d(r,NT,PCT,DT,N18,PC18,D18,NexpT,PCexpT,DexpT,Nexp18,PCexp18,Dexp18);
                
                if kappa0 >= thresh_cur
                    break;
                end
  
            end
            
            if kappa0 < thresh_cur
                accepted = true;
            end
            
        end
            
        epsTheta(j)     = kappa0;
        ThetaCur(j,:)   = Theta0;
        ModelCur(j)     = Mk0;
        SeedCur(j)      = seed - 9; % Corresponding to r = 1
            
    end
                
% Update extinct models
extinct = UpdateExtinct(ModelCur);
    
% Store
Weights(:,1) = WeightCur;
Eps(:,1)     = epsTheta;
Models(:,1)  = ModelCur;
Seeds(:,1)   = SeedCur;
Samples{1}   = ThetaCur;
    
%% SMC
for i = 2:n_thresh
   
    WeightLast  = WeightCur;
    ThetaLast   = ThetaCur;
    ModelLast   = ModelCur;
    
    thresh_cur  = thresh_seq(i);
       
    % LOOP THROUGH PARTICLES
    parfor j = 1:Nsamples
        
        % initialise seed, i = 1
        seed = SeedOffset + (i-1)*1e8 + (j-1) * 1e3;
        
        accepted = false;
        while ~accepted
            
            % Sample model
            Mk       = SampleModel(extinct);
                        
            % SAMPLE
            valid = false;
            while ~valid
                
                % RESAMPLE
                index   = SampleParticle(Mk,ModelLast,WeightLast);
                ThStar  = ThetaLast(index,:);
                
                % PERTURB (RESAMPLE OF PRIOR DENSITY ZERO)
                Theta   =  Ksample(ThStar,ModelNZ(Mk,:));
                
                % CHECK PRIOR DENSITY IS NON-ZERO
                valid   = ABC_Prior(Theta);
                
            end
            
            % INITIALISE EPS
            eps     = 0;
            
            % SIMULATE MODEL (LOK THROUGH DATA SETS)
            for r = 1:9
               
                % Determine Nmax
                Nmax = DetNmax(NexpT(r),thresh_cur - eps,ss_weights(1));
                
                % Run Model
                [NT,PCT,DT,N18,PC18,D18] = BinnyIBM([Theta,sigma,24,sigma],[1440,1900],IC{r},36,Nmax,seed);
                seed = seed + 1;
                
                if NT > Nmax
                    eps = eps + 1e9;
                    break;
                end
            
                % Increase eps
                eps  = eps + d(r,NT,PCT,DT,N18,PC18,D18,NexpT,PCexpT,DexpT,Nexp18,PCexp18,Dexp18);
                if eps >= thresh_cur
                    break;
                end
                
            end % end for r = 1:9
            
            if eps < thresh_cur
                accepted = true;
            end
            
        end % end while ~accepted
            
        % Update weights
        isMk            = ModelLast == Mk;
        WeightCur(j)    = 1 / dot(WeightLast(isMk),Kdensity(ThetaLast(isMk,:),Theta,ModelNZ(Mk,:)));
        
        % Update eps, theta, model and seed
        epsTheta(j)     = eps;
        ThetaCur(j,:)   = Theta;
        ModelCur(j)     = Mk;
        SeedCur(j)      = seed - 9; % Corresponding to r = 1
    
    end % end parfor j = 1:Nsamples
    
    % Update extinct models
    extinct = UpdateExtinct(ModelCur);
    
    % Store
    Weights(:,i) = WeightCur;
    Eps(:,i)     = epsTheta;
    Models(:,i)  = ModelCur;
    Seeds(:,i)   = SeedCur;
    Samples{i}   = ThetaCur;
    
end % end for i = 1:n_thresh
    
end


%% SUBROUTINES

% DETERMINE NMAX (POPULATION THAT WOULD RESULT IN IMMEDIATE REJECTION)
function Nmax = DetNmax(Nexp,Eps,Nweight)

    Nmax = ceil(Nexp * (1 + sqrt(Eps / Nweight)));

end

% SAMPLE MODEL
function Mi = SampleModel(Extinct)
    
    % Extinct = [0,0,0,0,0] indicates that no models are yet extinct
    Mindices = 1:5;
    Malive   = Mindices(~Extinct);
    
    Mi       = Malive(randi(length(Malive)));

end

% SAMPLE PARTICLE
function Pj = SampleParticle(Mi,ModelLast,WeightLast)

    % Set weights to zero where ModelLast ~= Mi
    Weights = WeightLast .* (Mi == ModelLast);
    Wcdf    = cumsum(Weights)/ sum(Weights);
    
    Pj      = find(rand() < Wcdf, 1, 'first');
    
end

% DISCREPENCY FUNCTION
function d = d(r,NT,PCT,DT,N18,PC18,D18,NexpT,PCexpT,DexpT,Nexp18,PCexp18,Dexp18)

    Imid    = [44,51,33,32,50,42,48,44,43];
    Dlookat = (Imid(r) - 20):(Imid(r) + 20);
    d = (Nexp18(r) - N18).^2 / Nexp18(r)^2 ...
                    + sum((PC18 - PCexp18(:,r)).^2) / sum(PCexp18(:,r).^2) ...
                    + sum((D18(Dlookat) - Dexp18(Dlookat,r)).^2) / sum(Dexp18(Dlookat,r).^2) ... 
                    + (NexpT(r) - NT).^2 / NexpT(r)^2 ...
                    + sum((PCT - PCexpT(:,r)).^2) / sum(PCexpT(:,r).^2) ...
                    + sum((DT(Dlookat) - DexpT(Dlookat,r)).^2) / sum(DexpT(Dlookat,r).^2);

end

% UPDATE EXTINCT
function e = UpdateExtinct(Models)

    for i = 1:5
        e(i) = ~any(Models == i);
    end

end
