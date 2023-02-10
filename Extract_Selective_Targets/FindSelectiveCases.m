function I = FindSelectiveCases(Models, LethalSets, CommonRxns, PathogenesID)
warning off
I = zeros(length(LethalSets), 1) - 1;
parfor i = 1 : length(LethalSets)
    I_case = zeros(length(Models), 1);
    for j = 1 : length(Models)
        I_case(j) = checkLethality(Models{j}, LethalSets(i,:));
    end
    
    if any(I_case(~PathogenesID) == 1) || any(I_case(PathogenesID) == 0)
        I(i) = 0;
    elseif all(PathogenesID) && all(I_case == 1)
        I(i) = 1;
    elseif all(~PathogenesID) && all(I_case == 0) && all(ismember(LethalSets(i,:), CommonRxns))
        I(i) = 1;
    elseif all(I_case(~PathogenesID) == 0) && all(I_case(PathogenesID) == 1) && all(ismember(LethalSets(i,:), CommonRxns))
        I(i) = 1;
    else
        I(i) = -1;
    end
end
