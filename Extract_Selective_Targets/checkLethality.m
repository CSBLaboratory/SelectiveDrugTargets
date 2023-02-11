function IsLethal = checkLethality(model, rxns)
% checkLethality finds the lethality of a reaction set for a genome scale
% model.
%   Inputs:
%       model:      genome scale model.
%       rxns:       a reaction list which should be deleted from the model.
%   Output:
%       IsLethal:   0 if the rxns set is not lethal for the model, 1 if the 
%                   set is found to be lethal.

solWT = optimizeCbModel(model);
for i = 1 : length(rxns)
    model = changeRxnBounds(model, rxns(i), 0, 'b');
end
    sol = optimizeCbModel(model);
    if sol.f < 0.01*solWT.f
        IsLethal = true;
    else
        IsLethal = false;
    end
end