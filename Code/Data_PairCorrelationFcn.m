function PC = PairCorrelationFcn(Data,ydist,dr,rmax)

    % Domain size
    H           = 1440;

    % X,Y
    X = Data(:,1);
    Y = Data(:,2);
    
    % Setup
    X1          = X(Y < ydist);
    Y1          = Y(Y < ydist);
    
    X2          = X(Y > (1900 - ydist));
    Y2          = Y(Y > (1900 - ydist));

    % Storage
    PC1         = PairCorr(X1,Y1);
    PC2         = PairCorr(X2,Y2);
  
    PC          = 0.5 * (PC1 + PC2);    
    
    
%%%% CALCULATE PC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PC = PairCorr(X,Y)
        
        N = length(X);
        R = 0:dr:rmax;
        
        PC  = zeros(1,length(R)-1);
        
        % Calculate distances
        for i = 1:N

            xtg     = X(i);
            ytg     = Y(i);

            dxtg    = abs(xtg - X([1:i-1,i+1:end]));
            dytg    = abs(ytg - Y([1:i-1,i+1:end]));

            dxtg    = min(dxtg, H     - dxtg);
            dytg    = min(dytg, ydist - dytg);

            dist    = sqrt(dxtg.^2 + dytg.^2);

            PC      = PC + histcounts(dist,R);

        end

        % Rescale D
        density     = N^2 / (ydist * H);
        A           = density * pi * dr * (2 * R(1:end-1) + dr);
        PC          = PC ./ A;

    end

end
