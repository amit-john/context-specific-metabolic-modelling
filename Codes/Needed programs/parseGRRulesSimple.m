function parsed = parseGRRulesSimple(grRules)
% Parses BiGG-style grRules containing AND/OR + parentheses.
% Returns a cell array of parsed rule strings.

parsed = cell(size(grRules));

for i = 1:numel(grRules)
    rule = strtrim(grRules{i});

    if isempty(rule)
        parsed{i} = '';
        continue;
    end

    % Standardize operators
    rule = strrep(rule, 'AND', 'and');
    rule = strrep(rule, 'OR', 'or');

    rule = strrep(rule, ' and ', ' & ');
    rule = strrep(rule, ' or ', ' | ');

    parsed{i} = rule;
end
end

