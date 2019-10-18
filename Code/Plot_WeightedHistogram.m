function f = Plot_WeightedHistogram(x,w,xlim,nbins)
%PLOT_WEIGHTEDHISTOGRAM Plot weighted histogram from data x with weights w

    % Bin edges
    edges   = linspace(xlim(1),xlim(2),nbins+1);
    bwidth  = diff(xlim) / nbins;

    % Normalise weights
    w       = w / sum(w);

    % Sum weights in each bin
    f       = zeros(size(x));
    for i = 1:nbins

        inbin   = x > edges(i) & x <= edges(i+1);
        f(i)    = sum(w(inbin));

        rectangle('Position',[edges(i),0,bwidth,f(i)],'FaceColor','b'); hold on;

    end


end