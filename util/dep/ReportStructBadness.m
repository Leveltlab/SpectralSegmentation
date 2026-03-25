function ReportStructBadness(struc, required)
% ReportStructBadness(struc, required)
% Reports what a struct is missing from the requested fieldnames
% 
% input:
%   - struc (struct)
%   - required ([n x 1] cell with strings): names of fields that should be
%   in there
% 
% output:
%   - text in command window explaining the badness
% 
% Leander de Kraker
% 2026-3-25
% 

strucFieldnames = fieldnames(struc);
foundFields = ismember(required, strucFieldnames);
if any(~foundFields)
    fprintf('%d/%d fields found. Missing fields:\n', sum(foundFields), length(required))
    for i = find(~foundFields)
        fprintf('%s\n', required{i})
    end
end