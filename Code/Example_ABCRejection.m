%% EXAMPLE ABC REJECTION

% Number samples (100,000 in supporting material document)
Nsamples    = 100;

% RNG seeds
PriorSeed   = 1;
SeedOffset  = 0;

% Simulate model multiple times, sample parameters from prior
[PriorSamples,Seeds,N,PC,D] = ABC_SimulateModel(Nsamples,PriorSeed,SeedOffset);

% Discrepency function weights
weights     = [1,1,1];

% ABC threshold (0.01 in supporting material document)
alpha       = 0.1;

% Load Experimental Results
[Nexp,PCexp,Dexp,centres,repnames] = Data_Grab();

%% Loop through summary statistics and prior samples
ErrorEach   = zeros(Nsamples,3);
Error       = zeros(Nsamples,1);

%Can vectorise this later!
for s = 1:Nsamples
    
    % N
    ErrorEach(s,1) = sum(((reshape(N(s,:,:),9,2) - Nexp)./ Nexp).^2,'all');
        
    for r = 1:9
        for t = 1:2
            
            % PC
            temp_PCsim = reshape(PC(s,r,t,:),1,20);
            temp_PCexp = reshape(PCexp(r,t,:),1,20);
            
            ErrorEach(s,2) = ErrorEach(s,2) + sum((temp_PCsim - temp_PCexp).^2) / sum(temp_PCexp.^2);
        
            % D
            Dlookat = (centres(r)-20):(centres(r)+20);
            temp_Dsim = reshape(D(s,r,t,:),1,80);
            temp_Dexp = reshape(Dexp(r,t,:),1,80);
            temp_Dsim = temp_Dsim(Dlookat);
            temp_Dexp = temp_Dexp(Dlookat);
            ErrorEach(s,3) = ErrorEach(s,3) + sum((temp_Dsim - temp_Dexp).^2) / sum(temp_Dexp.^2);
                        
        end
    end
end

clear temp_PCexp temp_PCsim temp_Dexp temp_Dsim
Error = sum(ErrorEach,2);

%% GET POSTERIOR

% Get best alpha = 0.01 samples
[~,Isort]   = sort(Error);
Ipost       = Isort(1:ceil(alpha * Nsamples));
PostSamples = PriorSamples(Ipost,:);

% Plot a quantile plot
subplot(2,3,1);
    p = 0:0.001:0.5;
    q = quantile(Error,p);
    plot(p,q,'LineWidth',3);
    title('Quantile Plot');
    xlabel('Quantile'); ylabel('Error');

% Plot Posteriors
titles = {'m','p','\gamma_m','\gamma_p','\gamma_b'};
for i = 1:5
    subplot(2,3,i+1);
        histogram(PostSamples(:,i));
        title(titles{i});
end