%% inspect_ml_model.m
% Print what the saved housekeeping ML model (estimate_calib_from_housekeeping_v4.m) uses.
files = {'C:\Users\hervo\OneDrive\Documents\MATLAB\ALC\CL_model_bagged_tree_202502.mat', ...
         'C:\Users\hervo\OneDrive\Documents\MATLAB\CL_model_baggeg_tree_202411.mat'};
for k = 1:numel(files)
    S = load(files{k});
    fn = fieldnames(S);
    fprintf('\n== %s : variables %s\n', files{k}, strjoin(fn, ', '));
    m = S.(fn{1});
    if isstruct(m)
        disp(fieldnames(m)');
        if isfield(m, 'RequiredVariables')
            fprintf('RequiredVariables: %s\n', strjoin(m.RequiredVariables, ', '));
        end
        if isfield(m, 'RegressionEnsemble')
            e = m.RegressionEnsemble;
            fprintf('class %s, method %s, NumTrained %d\n', class(e), e.Method, e.NumTrained);
            fprintf('PredictorNames: %s\n', strjoin(e.PredictorNames, ', '));
            fprintf('ResponseName: %s, NumObservations %d\n', e.ResponseName, e.NumObservations);
        end
        if isfield(m, 'About')
            disp(m.About);
        end
    else
        disp(class(m));
    end
end
fprintf('INSPECT_DONE\n');
