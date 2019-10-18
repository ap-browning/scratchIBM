function out = ABC_Prior(varargin)
%SAMPLEPRIOR Sample from the prior distribution. If parameters given as an
%input, return a bool indicating whether parameter is in the prior support.
%
% Example usage:
%   PriorSample  = ABC_Prior(1);
%   PriorSamples = ABC_Prior(1e5);
%   InSupport    = ABC_Prior([1,1,1,1,1]);
%   PriorLimits  = ABC_Prior;

% Settings
p_lim = [ 0.000   5.000;    % m
          0.020   0.050;    % p
         -2.000   2.000;    % gm
          0.000   0.020;    % gp
          0.000 100.000]' ;  % gb
          %8.000  30.000]';  % sigma
      
nparams = size(p_lim,2);
      
% Output limits
if nargin == 0
    out = p_lim';
else

% Generate varargin{1} samples
if length(varargin{1}) == 1
   
    out = rand(varargin{1},nparams) .* diff(p_lim) + p_lim(1,:);
    
% Check if in support
elseif length(varargin{1}) == 5
    
    out = ~(any(varargin{1} < p_lim(1,:)) || any(varargin{1} > p_lim(2,:)));
    
end


end

