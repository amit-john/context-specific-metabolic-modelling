function LocGini = giniTransform(mat, LT, UT)
% giniTransform
%   Compute per-gene local Gini thresholds and clamp to [LT, UT].
%
%   mat : (G x C) gene expression
%   LT  : lower bound
%   UT  : upper bound
%
%   LocGini : (G x C) matrix of thresholds

    [nGenes, nCtx] = size(mat);
    LocGini = zeros(nGenes, nCtx);

    for i = 1:nGenes
        row = sort(mat(i, :));

        if all(row == 0)
            % if completely zero, threshold = 0
            LocGini(i, :) = 0;
            continue;
        end

        n = numel(row);
        G = sum(((2*(1:n)) - n - 1) .* row) / (n * sum(row) + eps);  % Gini
        gc = G * 100;  % percentile

        LocGini(i, :) = prctile(mat(i, :), gc);
    end

    LocGini(LocGini < LT) = LT;
    LocGini(LocGini > UT) = UT;
end
