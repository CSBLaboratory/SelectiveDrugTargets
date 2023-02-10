function I = checkLethality(model, rxns)
solWT = optimizeCbModel(model);
for i = 1 : length(rxns)
    model = changeRxnBounds(model, rxns(i), 0, 'b');
end
    sol = optimizeCbModel(model);
    if sol.f < 0.01*solWT.f
        I = 1;
    else
        I = 0;
    end
end