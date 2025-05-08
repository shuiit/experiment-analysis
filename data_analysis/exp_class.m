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
        min_v
        response_time
        zero_v
        delta_angle
        name
        open_leg
        zero_v_z
        min_v_time
    end
    

    methods
        function obj = exp_class(exp,fps,insect_name)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = insect_name
            obj.fps = fps;
            header = exp.Properties.VariableNames;
            delta_t = exp.time(2)-exp.time(1);    
            obj.time_vec = round([min(exp.time):delta_t:max(exp.time) + delta_t]*100)/100;
            movie_numbers = unique(exp.mov_num);
            obj.exp_mat = nan(length(obj.time_vec),size(exp,2),length(movie_numbers));
            for mov_idx = [1:1:length(movie_numbers)]
                movie = exp(exp.mov_num == movie_numbers(mov_idx),:);
                try
                [~,exp_mat_idx] = intersect(obj.time_vec,round(movie.time*100)/100);
                [~,movie_idx] = intersect(round(movie.time*100)/100,obj.time_vec);
                obj.exp_mat(exp_mat_idx,:,mov_idx) = table2array(movie(movie_idx,:));
                catch
                    wakk = 2
                end
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
        
        function [prop_in_time,time] = get_prop_point_in_time(obj,prop_name,mov,interest_point)
            mov_name = sprintf('mov%d',mov);
            time_vec = obj.get_prop_mov('time',mov);
            idx = find(obj.(interest_point)(mov_name,:).Var1 == time_vec);
            prop = obj.get_prop_mov(prop_name,mov);
            prop_in_time = prop(idx);
            time = time_vec(idx);
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

        function percent = get_percent(obj,prop)
            prop = obj.(prop).Var1;
       
            percent = sum(prop~= 1000)*100/length(prop)
        end



        function percent = get_legs_percent(obj,prop)
                fly_open_legs = obj.(prop).Var1;
                idx_open_all_time = (fly_open_legs == 1);
                idx_not_opening = (fly_open_legs == 0);
                
                % only flies that spread legs during the video
                percent(1) = sum(idx_not_opening == 0 & idx_open_all_time == 0)*100/length(fly_open_legs(idx_open_all_time == 0));
                
                % flies that spread legs + crazy flies
                percent(2) = sum(idx_not_opening == 0 & idx_open_all_time == 0)*100/length(fly_open_legs);

                % only crazy flies
                percent(3) = sum(idx_open_all_time)*100/length(fly_open_legs);
        end

        function percent = get_percent_with_th(obj,prop,th,varargin)
            parser = inputParser;
            addParameter(parser,'th_st',true); % color of mean data
            parse(parser, varargin{:})

            prop = obj.(prop).Var1;
       
        if parser.Results.th_st == true
            percent = sum(( abs(prop) <th)*100/sum(prop ~= 1000))
        else
            percent = sum(( abs(prop) >th) & (prop ~= 1000))*100/sum(prop ~= 1000)

        end
        end

        function [mean_prop,std_prop]= get_mean_std_open_legs(obj,prop,varargin)
            idx_open_all_time = (obj.(prop).Var1 == 1);
            idx_not_opening = (obj.(prop).Var1 == 0);
            
            time_open_legs = obj.open_leg.Var1(idx_not_opening == 0 & idx_open_all_time == 0)
            mean_prop = mean(time_open_legs)
            std_prop = std(time_open_legs)
        end



        function mean_min_v = get_mean(obj,prop)
            prop = obj.(prop).Var1;
        
            mean_min_v = mean(prop(prop < 1000),1,'omitmissing')

        end
    end


end