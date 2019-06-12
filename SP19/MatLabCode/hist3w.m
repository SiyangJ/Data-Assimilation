function [N,C] = hist3w(X,w,nbins)

[~,Xedges,Yedges,binX,binY] = histcounts2(X(:,1),X(:,2),nbins);

C = {(Xedges(1:end-1)+Xedges(2:end))/2,(Yedges(1:end-1)+Yedges(2:end))/2};

N = accumarray([binX,binY],w) * size(X,1);

end

