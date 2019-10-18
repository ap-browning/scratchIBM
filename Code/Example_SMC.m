%% EXAMPLE ABC SMC (Model selection)

% Nsamples (5000 in main document)
Nsamples = 8;

% Sigma
sigma = 12; 

% Threshold sequence ([9.6,7.3,6.3,5.5,4.9,4.6,4.4] in main document)
thresh_seq = [9.6,7.3,6.3,5.5,4.9];

% Seed
SeedOffset = 0;

% Perform ABC
[Weights,Eps,Models,Seeds,PostSamples] = ABC_SMC(Nsamples,sigma,thresh_seq,SeedOffset);

% Plot results
clf

% Plot model indices
subplot(2,3,1);
    ModelX = 1:5;
    ModelY = zeros(1,5);
    for i = 1:5
        ModelY(i) = sum(Models(:,end) == i);
    end
    bar(ModelX,ModelY);
    
% Plot posterior samples from highest density model (weighted)
Mbest  = mode(Models(:,end));
IMbest = Models(:,end) == Mbest;
PriorLimits = ABC_Prior;
for i = 1:5
    subplot(2,3,i+1);
        Plot_WeightedHistogram(PostSamples{end}(IMbest,i),Weights(IMbest,end),PriorLimits(i,:),10);
end