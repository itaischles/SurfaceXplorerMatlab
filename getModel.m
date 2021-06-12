function model = getModel(modelName,experimentData,tableData)

% check which variables are fixed
fixVars = find(cell2mat(tableData(:,4))==true);

switch modelName
    case 'A->Gnd'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     -1/t(3), 1);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}];
        model.lb = [-5 0.1 0];
        model.ub = [5 5 30000];
        model.species_names = {'A'};
    case 'A->B'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)  0;
                                       1/t(3)  0;], [1,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}];
        model.lb = [-5 0.1 0];
        model.ub = [5 5 30000];
        model.species_names = {'A', 'B'};
    case 'A->B->Gnd'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)  0;
                                       1/t(3)  -1/t(4);], [1,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}];
        model.lb = [-5 0.1 0 0];
        model.ub = [5 5 30000 30000];
        model.species_names = {'A', 'B'};
    case 'A->B->C'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)    0      0;
                                       1/t(3)  -1/t(4)  0;
                                       0        1/t(4)  0;], [1,0,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}];
        model.lb = [-5 0.1 0 0];
        model.ub = [5 5 30000 30000];
        model.species_names = {'A', 'B', 'C'};
    case 'A->B->C->Gnd'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)    0      0;
                                       1/t(3)  -1/t(4)  0;
                                       0        1/t(4)  -1/t(5);], [1,0,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}, tableData{5,2}];
        model.lb = [-5 0.1 0 0 0];
        model.ub = [5 5 30000 30000 30000];
        model.species_names = {'A', 'B', 'C'};
    case 'A->B->C->D'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)    0      0           0;
                                       1/t(3)  -1/t(4)  0           0;
                                       0        1/t(4)  -1/t(5)     0;
                                       0         0      1/t(5)      0;], [1,0,0,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}, tableData{5,2}];
        model.lb = [-5 0.1 0 0 0];
        model.ub = [5 5 30000 30000 30000];
        model.species_names = {'A', 'B', 'C', 'D'};
    case 'A->B->C->D->Gnd'
        model.kinmod=@(t)buildkineticmodel(experimentData,t(1),t(2),...
                                     [-1/t(3)    0      0           0;
                                       1/t(3)  -1/t(4)  0           0;
                                       0        1/t(4)  -1/t(5)     0;
                                       0         0      1/t(5)      -1/t(6);], [1,0,0,0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}, tableData{5,2}, tableData{6,2}]; 
        model.lb = [-5 0.1 0 0 0 0];
        model.ub = [5 5 30000 30000 30000 30000];
        model.species_names = {'A', 'B', 'C', 'D'};
    case 'Biexciton decay (1D) -> Gnd'
        dt = @(k,t) 1./(1 + k(1)*sqrt(t));
        model.kinmod=@(v)buildmodel(experimentData,dt,v(1),v(2),v(3));
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}];
        model.lb = [-5 0.1 0];
        model.ub = [5 5 30000];
        model.species_names = {'A'};
    case 'Biexciton decay (1D) -> B'
        dt = @(k,t) [1./(1 + k(1)*sqrt(t)); 1-1./(1 + k(1)*sqrt(t))];
        model.kinmod=@(v)buildmodel(experimentData,dt,v(1),v(2),v(3));
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}];
        model.lb = [-5 0.1 0];
        model.ub = [5 5 30000];
        model.species_names = {'A','B'};
    case 'Biexciton decay (1D) -> Gnd with A->B'
        dadt=@(a,k,t)[  (-k(2)-0.5*k(1)/(sqrt(t+1e-2)*(1+k(1)*sqrt(t+1e-2))))*a(1);
                        k(2)*a(1);];
        model.kinmod=@(v)buildmodeldiff(experimentData,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}];
        model.lb = [-5 0.1 0 0];
        model.ub = [5 5 30000 30000];
        model.species_names = {'A','B'};
    case 'Biexciton decay (1D) -> Gnd with A->Gnd'
        dadt=@(a,k,t) (-k(2)-0.5*k(1)/(sqrt(t+1e-2)*(1+k(1)*sqrt(t+1e-2))))*a(1);
        model.kinmod=@(v)buildmodeldiff(experimentData,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1]);
        model.initialguess = [tableData{1,2}, tableData{2,2}, tableData{3,2}, tableData{4,2}];
        model.lb = [-5 0.1 0 0];
        model.ub = [5 5 30000 30000];
        model.species_names = {'A'};
end

% fix variables by enforcing upper and lower bounds for fitting equal to
% user defined value
if ~isempty(fixVars)
    model.lb(fixVars) = cell2mat(tableData(fixVars,2));
    model.ub(fixVars) = cell2mat(tableData(fixVars,2));
end

end

