function [model, addedRxns] = addTransportsForExchanges(model, exchanges, fromSuffix, toSuffix, subsystem)
% addTransportsForExchanges
%   For each exchange reaction in "exchanges", automatically add a
%   reversible transport met[e] <=> met[c] if both mets exist and no
%   transport already exists.
%
%   model      : COBRA model
%   exchanges  : cell array of exchange rxn IDs (e.g. {'EX_glc__D_e',...})
%   fromSuffix : compartment suffix of external mets (default '_e')
%   toSuffix   : compartment suffix of internal mets (default '_c')
%   subsystem  : subsystem name for new transport reactions
%
%   Returns:
%       model      : updated model
%       addedRxns  : cell array of IDs of newly added transport reactions

if nargin < 3 || isempty(fromSuffix)
    fromSuffix = '_e';
end
if nargin < 4 || isempty(toSuffix)
    toSuffix = '_c';
end
if nargin < 5 || isempty(subsystem)
    subsystem = 'DMEM_transports';
end

addedRxns = {};

for k = 1:numel(exchanges)
    exID = exchanges{k};
    exIdx = find(strcmp(model.rxns, exID));

    if isempty(exIdx)
        fprintf('Exchange %s not found. Skipping.\n', exID);
        continue;
    end

    % Find the metabolite in the EX reaction (should be exactly 1)
    metIdx = find(model.S(:, exIdx) ~= 0);
    if numel(metIdx) ~= 1
        fprintf('Exchange %s has %d metabolites. Skipping.\n', exID, numel(metIdx));
        continue;
    end
    met_e = model.mets{metIdx};

    % Check suffix
    if ~endsWith(met_e, fromSuffix)
        fprintf('Met %s in %s does not end with %s. Skipping.\n', met_e, exID, fromSuffix);
        continue;
    end

    % Build corresponding internal metabolite ID (e.g. glc__D_c)
    base = extractBefore(met_e, strlength(met_e) - strlength(fromSuffix) + 1);
    met_c = [base toSuffix];

    if ~any(strcmp(model.mets, met_c))
        fprintf('No internal metabolite %s for %s. Skipping transport.\n', met_c, met_e);
        continue;
    end

    % Check if a transport between met_e and met_c already exists
    met_e_idx = find(strcmp(model.mets, met_e));
    met_c_idx = find(strcmp(model.mets, met_c));
    candidateRxns = find(any(model.S([met_e_idx, met_c_idx], :) ~= 0, 1));
    already = false;
    for r = candidateRxns
        metsInR = model.mets(model.S(:, r) ~= 0);
        if any(strcmp(metsInR, met_e)) && any(strcmp(metsInR, met_c))
            already = true;
            break;
        end
    end
    if already
        fprintf('Transport %s <-> %s already exists. Skipping.\n', met_e, met_c);
        continue;
    end

    % Create a new reaction ID
    rxnID = ['TR_' base];
    if any(strcmp(model.rxns, rxnID))
        rxnID = [rxnID '_1'];
    end

    % Add reversible transport: met_e <=> met_c
    model = addReaction(model, rxnID, ...
        'metaboliteList',   {met_e, met_c}, ...
        'stoichCoeffList',  [-1, 1], ...
        'reversible',       true, ...
        'subSystem',        subsystem);

    fprintf('Added transport %s: %s <=> %s\n', rxnID, met_e, met_c);
    addedRxns{end+1,1} = rxnID;
end
end
