function run_mov_dir(files_in_dir,exp_name,exp_path,idx_exp,length_exp,path_to_save)
% run on all files in movie dir, check if files exist, load hull and make
% new dir. Than, generate the csv

f = waitbar(0)
for idx = 1:length(files_in_dir)
    waitbar(idx/length(files_in_dir),f,sprintf('mov %d out of %d (%d/%d)',idx,length(files_in_dir),idx_exp,length_exp))
    if contains(files_in_dir(idx).name, 'mov')
        mov = str2double(regexp(files_in_dir(idx).name, '\d*', 'match'));
        file_to_save = sprintf('%s_mov_%d',exp_name,mov)
        
        % check if files exist, load hull and make new dir
%         if sum(cellfun(@isfile,strcat([path_to_save,exp_name,'\\',file_to_save],{'_vectors.csv','_angles_cm.csv','_ew_to_lab_rotmat.csv'}))) == 3
%             continue
%         end
        loaders = loaders_class(exp_path,mov,'','hullfile','//hull_op//');
        [hull] = loaders.load_shull(1);
        if  isstruct(hull) == 0 || isfield(hull.rightwing.vectors,'nrml') == 0 || nnz(~isnan(hull.body.angles.pitch)) < 1000
            continue
        end
        if isfolder([path_to_save,exp_name,'\\']) == 0
            mkdir([path_to_save,exp_name,'\\'])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % generate and save the csv %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hull_to_csv = HullToCsv(hull,[path_to_save,exp_name,'\\',file_to_save]);
    hull_to_csv.angles_to_csv();
    hull_to_csv.vectors_to_csv();
    
    end
end
end



