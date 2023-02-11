function isSelective = findSelectiveCases(Models, LethalSets, CommonRxns, PathogenesID)
% findSelectiveCases finds the selective solutions
%   Inputs:     
%       Models,         constraint based models of the understudy microorganisms
%       LethalSets,     lethal sets identified for each microorganism
%       CommonRxns,     the list of common reactions among all understudy microorganisms
%       PathogenesID,   a vector representing the status of the understudy
%                       micrtootganisms
%   Output:
%       isSelective,    0 if the case is not selective, 1 if the case is
%                       selective, -1 if the reactions of the lethal set is not in CommonRxns

warning off
isSelective = zeros(length(LethalSets), 1) - 1;
parfor i = 1 : length(LethalSets)
    isLethal = zeros(length(Models), 1);
    for j = 1 : length(Models)
        isLethal(j) = checkLethality(Models{j}, LethalSets(i,:));
    end
    
    if any(isLethal(~PathogenesID)) || any(~isLethal(PathogenesID))
        isSelective(i) = 0;
    elseif all(PathogenesID) && all(isLethal)
        isSelective(i) = 1;
    elseif all(~PathogenesID) && all(~isLethal) && all(ismember(LethalSets(i, :), CommonRxns))
        isSelective(i) = 1;
    elseif all(~isLethal(~PathogenesID)) && all(isLethal(PathogenesID)) && all(ismember(LethalSets(i,:), CommonRxns))
        isSelective(i) = 1;
    else
        isSelective(i) = -1;
    end
end
