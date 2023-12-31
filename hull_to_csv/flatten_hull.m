

% saves fields in structure to csv
name_of_features_wings = {'span','chord','nrml'}
name_of_features_body = {'strkPlan','X','Y','Z'}
name_of_general_features = {'frames'}
name_of_video_features = {'timeframe'}
name_of_features_wings_angle = {'phi','theta','psi'}
name_of_features_body_angle = {'pitch','yaw','roll'}

dir_of_exps = 'H:\\My Drive\\dark 2022\\'
path_to_save = [dir_of_exps,'csv_dark\\']

exps_names = dir(dir_of_exps)
for idx_exp_name = 1:1:length(exps_names)
    exp_name = exps_names(idx_exp_name).name
    path=sprintf('%s%s\\hull\\hull_Reorder\\',dir_of_exps,exp_name);
    files_in_dir = dir(path);
    if isempty(files_in_dir) == 0
        run_mov_dir(files_in_dir,exp_name,path,idx_exp_name,length(exps_names),path_to_save);
    end
end
F = findall(0,'type','figure','tag','TMWWaitbar')
delete(F)
