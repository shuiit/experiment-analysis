classdef exp_class
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        time_vec
        mov_num
        insect_prop
        exp_mat
        t0_idx
    end
    

    methods
        function obj = exp_class(exp)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
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
            obj.t0_idx = find(obj.time_vec == 0)

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
    end
end