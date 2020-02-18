function mergedModel = Cheng_mergeConditionSpecificModels(model1,model2)
% mergedModel = Cheng_mergeConditionSpecificModels(model1,model2)
% model1 will overlap model2 in all cases except for gene associations
% where the union is used
mergedModel.description = 'Merged model';
mergedModel.id = 'Merged model';
mergedModel.annotation = model1.annotation;
mergedModel.rxns = union(model1.rxns,model2.rxns);
mergedModel.genes = union(model1.genes,model2.genes);
mergedModel.mets = union(model1.mets,model2.mets);
mergedModel.S = sparse(length(mergedModel.mets),length(mergedModel.rxns));
[A1,B1] = ismember(mergedModel.mets,model2.mets);
[A2,B2] = ismember(mergedModel.rxns,model2.rxns);
mergedModel.S(A1,A2) = model2.S(B1(A1),B2(A2));
[A1,B1] = ismember(mergedModel.mets,model1.mets);
[A2,B2] = ismember(mergedModel.rxns,model1.rxns);
mergedModel.S(A1,A2) = model1.S(B1(A1),B2(A2));
rxnGeneMat1 = sparse(length(mergedModel.rxns),length(mergedModel.genes));
rxnGeneMat2 = rxnGeneMat1;
[A1,B1] = ismember(mergedModel.rxns,model2.rxns);
[A2,B2] = ismember(mergedModel.genes,model2.genes);
rxnGeneMat1(A1,A2) = model2.rxnGeneMat(B1(A1),B2(A2));
[A1,B1] = ismember(mergedModel.rxns,model1.rxns);
[A2,B2] = ismember(mergedModel.genes,model1.genes);
rxnGeneMat2(A1,A2) = model1.rxnGeneMat(B1(A1),B2(A2));
mergedModel.rxnGeneMat = logical(rxnGeneMat1+rxnGeneMat2);
[A1,B1] = ismember(mergedModel.rxns,model1.rxns);
[A2,B2] = ismember(mergedModel.rxns,model2.rxns);
mergedModel.ub = zeros(length(mergedModel.rxns),1);
mergedModel.ub(A2) = model2.ub(B2(A2));
mergedModel.ub(A1) = model1.ub(B1(A1));
mergedModel.lb = zeros(length(mergedModel.rxns),1);
mergedModel.lb(A2) = model2.lb(B2(A2));
mergedModel.lb(A1) = model1.lb(B1(A1));

mergedModel.grRules = cell(length(mergedModel.rxns),1);
for i = 1:length(mergedModel.grRules)
    mergedModel.grRules{i} = strjoin(mergedModel.genes(mergedModel.rxnGeneMat(i,:)),' or ');
end

mergedModel.rxnComps = zeros(length(mergedModel.rxns),1);
mergedModel.rxnComps(A2) = model2.rxnComps(B2(A2));
mergedModel.rxnComps(A1) = model1.rxnComps(B1(A1));
mergedModel.rxnNames = cell(length(mergedModel.rxns),1);
mergedModel.rxnNames(A2) = model2.rxnNames(B2(A2));
mergedModel.rxnNames(A1) = model1.rxnNames(B1(A1));
mergedModel.subSystems = cell(length(mergedModel.rxns),1);
mergedModel.subSystems(A2) = model2.subSystems(B2(A2));
mergedModel.subSystems(A1) = model1.subSystems(B1(A1));
mergedModel.eccodes = cell(length(mergedModel.rxns),1);
mergedModel.eccodes(A2) = model2.eccodes(B2(A2));
mergedModel.eccodes(A1) = model1.eccodes(B1(A1));
mergedModel.c = zeros(length(mergedModel.rxns),1);
mergedModel.c(A2) = model2.c(B2(A2));
mergedModel.c(A1) = model1.c(B1(A1));
mergedModel.rev = zeros(length(mergedModel.rxns),1);
mergedModel.rev(A2) = model2.rev(B2(A2));
mergedModel.rev(A1) = model1.rev(B1(A1));
[A1,B1] = ismember(mergedModel.mets,model1.mets);
[A2,B2] = ismember(mergedModel.mets,model2.mets);
if isfield(model1,'unconstrained')
    mergedModel.unconstrained = zeros(length(mergedModel.mets),1);
    mergedModel.unconstrained(A2) = model2.unconstrained(B2(A2));
    mergedModel.unconstrained(A1) = model1.unconstrained(B1(A1));
end
mergedModel.b = zeros(length(mergedModel.mets),1);
mergedModel.b(A2) = model2.b(B2(A2));
mergedModel.b(A1) = model1.b(B1(A1));
mergedModel.metNames = cell(length(mergedModel.mets),1);
mergedModel.metNames(A2) = model2.metNames(B2(A2));
mergedModel.metNames(A1) = model1.metNames(B1(A1));
mergedModel.metComps = zeros(length(mergedModel.mets),1);
mergedModel.metComps(A2) = model2.metComps(B2(A2));
mergedModel.metComps(A1) = model1.metComps(B1(A1));
mergedModel.metFormulas = cell(length(mergedModel.mets),1);
mergedModel.metFormulas(A2) = model2.metFormulas(B2(A2));
mergedModel.metFormulas(A1) = model1.metFormulas(B1(A1));
%mergedModel.metMiriams = cell(length(mergedModel.mets),1);
%mergedModel.metMiriams(A2) = model2.metMiriams(B2(A2));
%mergedModel.metMiriams(A1) = model1.metMiriams(B1(A1));
mergedModel.comps = model1.comps;
mergedModel.compNames = model1.compNames;
mergedModel.compOutside = model1.compOutside;
[A1,B1] = ismember(mergedModel.genes,model1.genes);
[A2,B2] = ismember(mergedModel.genes,model2.genes);
mergedModel.geneComps = zeros(length(mergedModel.genes),1);
mergedModel.geneComps(A2) = model2.geneComps(B2(A2));
mergedModel.geneComps(A1) = model1.geneComps(B1(A1));


end
