function varargout = Data_Grab(varargin)
%% GrabData
%  Process experimental data to obtain matrices Nexp, PCexp, Dexp where
%  
%   size(Nexp)  = [  9   2  ]
%   size(PCexp) = [  9   2   80  ]
%   size(Dexp)  = [  9   2   20  ]
%                    |   |   |
%                  Reps  |   |
%                      Times |
%                           DimSS
%
%

repnames = {'8000_E2', '8000_G2', '10000_H1',  ...
            '10000_H2','8000_F2', '10000_A2', ...
            '12000_H3','12000_F3','12000_D2'};
tvec     = {'30h','48h'};

if nargin == 1
    
    IC = cell(1,9);
    for r = 1:9
        IC{r} = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_12h.csv']);
    end
    
    varargout{1} = IC;
    
else

    % Loop through reps
    for r = 1:9

        for t = 1:2

            % Load data
            X = csvread(['../Data/',repnames{r},'/PC3_',repnames{r},'_',tvec{t},'.csv']);

            % Grab summary statistics
            [Ntemp,PCtemp,Dtemp] = Data_SummaryStatistics(X);

            % Fill outputs
            N(r,t)      = Ntemp;
            PC(r,t,:)   = PCtemp;
            D(r,t,:)    = Dtemp;

        end
        
        varargout{1} = N;
        varargout{2} = PC;
        varargout{3} = D;

    end

    if nargout > 3

        centres = [44,51,33,32,50,42,48,44,43];

    end

    if nargout == 4

        varargout{4} = centres;

    elseif nargout == 5

        varargout{4} = centres;
        varargout{5} = repnames;

    end

end

end