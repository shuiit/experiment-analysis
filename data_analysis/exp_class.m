classdef exp_class
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        time_vec
        mov_num
        insect_prop
        exp_mat
        t0_idx
        fps
    end
    

    methods
        function obj = exp_class(exp,fps)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.fps = fps;
            header = exp.Properties.VariableNames;
            delta_t = exp.time(2)-exp.time(1);    
            obj.time_vec = round([min(exp.time):delta_t:max(exp.time) + delta_t]*100)/100;
            movie_numbers = unique(exp.mov_num);
            obj.exp_mat = nan(length(obj.time_vec),size(exp,2),length(movie_numbers));
            for mov_idx = [1:1:length(movie_numbers)]
                movie = exp(exp.mov_num == movie_numbers(mov_idx),:);
                obj.exp_mat(ismember(obj.time_vec,round(movie.time*100)/100),:,mov_idx) = table2array(movie);
            end
            obj.insect_prop = dictionary([string(header)],[1:1:length(header)]);
            obj.mov_num = dictionary(movie_numbers',[1:1:length(movie_numbers)]);
            obj.t0_idx = find(obj.time_vec == 0);
            obj = obj.sub_t0_frame('yaw','yaw_0');

        end

        function prop = get_prop(obj,prop)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
                prop = squeeze(obj.exp_mat(:,obj.insect_prop(prop),:));
        end

        function mov = get_mov(obj,mov)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
                mov = obj.exp_mat(:,:,obj.mov_num(mov));
        end

        function prop_mov = get_prop_mov(obj,prop,mov)
            prop_mov = obj.exp_mat(:,obj.insect_prop(prop),obj.mov_num(mov));
        end

        
        function fvec = calc_ksdensity(obj,prop_name,time_vec)
            prop = obj.get_prop(prop_name);
            idx_to_violin = find(ismembertol(obj.time_vec,time_vec,1/(obj.fps/1000) -1/(obj.fps/1000)/2 ,'DataScale', 1))
            prop_time = prop(idx_to_violin,:);
            for j = 1:1:size(prop_time,1)
                [f,xi] = ksdensity(prop_time(j,:));
                fvec(j) = max(f);
            end
        end
        
        function obj = sub_t0_frame(obj,prop_to_sub,header)
            prop = obj.get_prop(prop_to_sub);
            prop_t0 = prop - prop(obj.t0_idx,:);
            size_exp_mat = size(obj.exp_mat);
            obj.exp_mat(:,size_exp_mat(2) + 1,:) = prop_t0;
            obj.insect_prop = insert(obj.insect_prop,header,size_exp_mat(2) + 1);

        end
    end


end