%Merge the 139 GBM models 
merged_model = Cheng_mergeConditionSpecificModels(model_1, model_2);   
for i = 3:139     
  a=eval(['model_' num2str(i)]);     
  merged_model = Cheng_mergeConditionSpecificModels(merged_model,a); 
end 
 
%Essentiality analysis 
EXmodel = LTM(merged_model); 
EXmodel = setParam(EXmodel,'obj','HCC_biomass',1); 
EA = FastGeneSL(EXmodel_ind,0.01,1,1); 
 
%Check removal of essential genes in normal brain model
[Lia,indeces]=ismember(SingleLethalGenes,model_normal_brain.genes);     
for i = 1:24     
  rxnsToRemove=find(model_normal_brain.rxnGeneMat(:,indeces(i)));     
  reducedModel=removeReactions(model_normal_brain,rxnsToRemove);              
  [taskReport]=checkTasks(reducedModel,'common_tasks_backup.xlsx',true ,true); 
end   


