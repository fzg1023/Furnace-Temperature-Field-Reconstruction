function Phi = buildRBFDesign(cx,cy,sigma,usePR)

    if nargin<4, usePR=false; end
    n    = numel(cx);
    PhiR = zeros(n,n);
    for i = 1:n
        r2        = (cx-cx(i)).^2 + (cy-cy(i)).^2;
        PhiR(i,:) = exp(-r2/(2*sigma^2));
    end
    if ~usePR
        Phi = PhiR;
    else
        Phi = [PhiR , ones(n,1) , cx , cy];      % 16×19
    end
end
