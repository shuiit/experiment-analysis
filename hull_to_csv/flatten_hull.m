
addpath('G:\Documents\micro-flight-lab\hull_reconstruction_git')
% saves fields in structure to csv
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
