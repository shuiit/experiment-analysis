classdef plotter
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        gray_mat
        col_mat
    end
    

    methods
        function obj = plotter()
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
           
            obj.gray_mat = colormap(gray);
            obj.col_mat = colormap(turbo);
        end

        function all_data_plot(obj,time,prop,varargin)
            parser = inputParser;
            addParameter(parser,'alpha',0.5); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})

            plot(time,prop,'LineWidth',1.5,'color',[obj.gray_mat(200,:),parser.Results.alpha],'HandleVisibility','off')
        end

        function cluster_plot_mean(obj,prop,insect,varargin)
            parser = inputParser;
            addParameter(parser,'colors',[180,85]); % number of camera pointing in Z lab axis
            addParameter(parser,'color_all_data',[180,85]); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            [data_a,data_b] = obj.cluster_plot(insect,'z_vel',0,prop,parser.Results.color_all_data,'alpha',0.03);
            
            obj.mean_plot(insect,data_a,prop,'plot_all_data',0,'color_idx',parser.Results.colors(1));
            obj.mean_plot(insect,data_b,prop,'plot_all_data',0,'color_idx',parser.Results.colors(2));

        end


        function ax = mean_plot(obj,insect,prop,label_y,varargin)
            parser = inputParser;
            addParameter(parser,'color_idx',13); % number of camera pointing in Z lab axis
            addParameter(parser,'plot_all_data',1); % number of camera pointing in Z lab axis
            addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
            addParameter(parser,'input_prop',0); % number of camera pointing in Z lab axis

            parse(parser, varargin{:})
            time = insect.time_vec;
            

            if ischar(prop) == 1
                prop = insect.get_prop(prop);
            end
            idx_to_plot = find(time > -10);

            mean_prop = mean(prop,2,'omitnan');
            ax = plot(time(idx_to_plot),mean_prop(idx_to_plot),'LineWidth',3,'color',obj.col_mat(parser.Results.color_idx,:))
            xlabel(parser.Results.xlabel)
            ylabel(label_y)
            
            
        end




        function mov_plot(obj,insect,prop_name,mov_num,label_y,varargin)
            parser = inputParser;
            addParameter(parser,'color_idx',13); % number of camera pointing in Z lab axis
            addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
            addParameter(parser,'LineWidth',3); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})

            time = insect.time_vec;
            prop = insect.get_prop(prop_name);
            idx_to_plot = find(time > -10);


            mov_prop = insect.get_prop_mov(prop_name,mov_num);
            ax1 = plot(time(idx_to_plot),mov_prop(idx_to_plot),'--','LineWidth',parser.Results.LineWidth,'color',obj.col_mat(parser.Results.color_idx,:))
            plot(time(idx_to_plot),mov_prop(idx_to_plot),'LineWidth',parser.Results.LineWidth,'color',[obj.col_mat(parser.Results.color_idx,:),0.4],'HandleVisibility','off')
            xlabel(parser.Results.xlabel)
            ylabel(label_y)

        end


        function mov_mean_plot(obj,insect,prop_name,mov_num,label_y,varargin)
            parser = inputParser;
            addParameter(parser,'color_idx',13); % number of camera pointing in Z lab axis
            addParameter(parser,'color_idx_mov',50); % number of camera pointing in Z lab axis
            addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
            addParameter(parser,'plot_cluster',0); % number of camera pointing in Z lab axis
            addParameter(parser,'cluser_mean_color',[25,85]); % number of camera pointing in Z lab axis
            addParameter(parser,'color_all_data',[25,85]); % number of camera pointing in Z lab axis
           parse(parser, varargin{:})

            time = insect.time_vec;
            prop = insect.get_prop(prop_name);
            
            if parser.Results.plot_cluster == 1
                obj.cluster_plot_mean(prop_name,insect,'colors',parser.Results.cluser_mean_color,'color_all_data',parser.Results.color_all_data);
            else
                obj.all_data_plot(time,prop);hold on
                obj.mean_plot(insect,prop,label_y,'plot_all_data',0,'color_idx',parser.Results.color_idx);
            end
            obj.mov_plot(insect,prop_name,mov_num,label_y,'color_idx',parser.Results.color_idx_mov);
            xlabel(parser.Results.xlabel);
            ylabel(label_y);
            if parser.Results.plot_cluster == 1
                legend({'mean_A','mean_B','model'},'EdgeColor','None','Orientation','horizontal','Box','off')
            else
                legend({'mean','model'},'EdgeColor','None','Orientation','horizontal','Box','off')
            end

        end


        function pert_plot(obj,pert,plot_zero_x,plot_zero_y,sp)
            if plot_zero_y == 1
            yline(0,'--k','LineWidth',1.5,'handlevisibility','off');
            end
            if plot_zero_x == 1
            xline(0,'--r','LineWidth',1.5,'handlevisibility','off','color',[0.8500 0.3250 0.0980]);
            if pert ~= 0
            xline(pert,'--','LineWidth',1.5,'handlevisibility','off','color',[0.2660 0.6740 0.1880]	);
            end
            end

            if pert == 0
                pert = sp(1).XLim(2);
            end
            
            rec =  [0,sp(1).YLim(1),pert,sp(1).YLim(2)-sp(1).YLim(1)];
            rectangle('Position',rec,...
                'EdgeColor','none','FaceColor',[0,0,0,0.1]); % Plots the rectangle
        end

        function [data_a,data_b] = cluster_by_prop(obj,prop_to_cluster,time_to_cluster,insect,prop_to_plot)

            idx_t = find(insect.time_vec == time_to_cluster);
            cluster_prop = insect.get_prop(prop_to_cluster);
            plot_prop = insect.get_prop(prop_to_plot);
            data_a = plot_prop(:,cluster_prop(idx_t,:) >= 0);
            data_b = plot_prop(:,cluster_prop(idx_t,:) < 0);

        end
        
        function [data_a,data_b] = cluster_plot(obj,insect,prop_to_cluster,time_to_cluster,prop_to_plot,color_all_data,varargin)
            parser = inputParser;
            addParameter(parser,'alpha',0.1); % number of camera pointing in Z lab axis
           parse(parser, varargin{:})
            [data_a,data_b] = cluster_by_prop(obj,prop_to_cluster,time_to_cluster,insect,prop_to_plot);

            plot(insect.time_vec,data_a,'color',[obj.col_mat(color_all_data(1),:),parser.Results.alpha],'LineWidth',2,'HandleVisibility','off');hold on
            plot(insect.time_vec,data_b,'color',[obj.col_mat(color_all_data(2),:),parser.Results.alpha],'LineWidth',2,'HandleVisibility','off');hold on

        end


        function plot_pitch_for_vel_z_vel_model_mov_mean(obj,insect,exp_name,mov,pert,color_idx,color_idx_mov,max_time_xax,cluser_mean_color,color_all_data)
            mos_flag = 0;
            figure();
            ax1 = subplot(3,1,1)
            obj.mov_mean_plot(insect.(exp_name),'pitch',mov,'pitch [deg]','color_idx',color_idx,'color_idx_mov',color_idx_mov)
            obj.pert_plot(pert,0,1,ax1)

            
            ax2 = subplot(3,1,2)
            obj.mov_mean_plot(insect.(exp_name),'forward_vel',mov,'V_f_w_d [m/s]','color_idx',color_idx,'color_idx_mov',color_idx_mov)
            obj.pert_plot(pert,0,1,ax2)
            
            ax3 = subplot(3,1,3)
            if strcmp(insect.name,'mosquito')
                mos_flag = 1;
            end
            obj.mov_mean_plot(insect.(exp_name),'z_vel',mov,'V_z [m/s]','color_idx',color_idx,'color_idx_mov',color_idx_mov,'plot_cluster',mos_flag,'cluser_mean_color',cluser_mean_color,'color_all_data',color_all_data)
            obj.pert_plot(pert,0,1,ax3)

            
            ax1.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            ax2.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            ax3.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            linkaxes([ax1,ax2,ax3],'x')
            pert = split(exp_name,'_')
            sgtitle([insect.name,' ',pert{2}])


            set([ax1,ax2,ax3], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on
            % set(gca,'fontsize',20)

        end



        function plot_mean_subplot(obj,fly,mos,exp_name,pert,color_idx,prop_name,max_time_xax)
            % 
            % figure
            ax1 = subplot(1,1,1)
            obj.mean_plot(fly.(exp_name),prop_name,'pitch [deg]','color_idx',color_idx(1))
            obj.pert_plot(pert,0,1,ax1)
            title('fly')
            
            % ax2 = subplot(2,1,2)
            % obj.mean_plot(mos.(exp_name),prop_name,'pitch [deg]','color_idx',color_idx(2))
            % obj.pert_plot(pert,0,1,ax2)
            % title('mosquito')

            ax1.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            % ax2.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            % linkaxes([ax1,ax2],'x')
            
        end

        

        function h = violin_plot(obj,insect,prop_name,time_to_violin,f_norm_vec,label_y,color,legend_val,varargin)
            parser = inputParser;
            addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
             addParameter(parser,'scatter_loc',5); % number of camera pointing in Z lab axis
              addParameter(parser,'box_xdata',5); % number of camera pointing in Z lab axis
          parse(parser, varargin{:})
            idx_to_violin = find(ismembertol(insect.time_vec,time_to_violin,1/(insect.fps/1000) -1/(insect.fps/1000)/2 ,'DataScale', 1))
            x_names = insect.time_vec(idx_to_violin)
            
            prop = insect.get_prop(prop_name);
            prop_time = prop(idx_to_violin,:);
            h = daviolinplot(prop_time',x_names,f_norm_vec,'colors',color,'box',0,'whiskers',0,'violinalpha',0.5,'bins',15, ...
                'outliers',0,'legend',legend_val,'scatter',parser.Results.scatter_loc,'scattercolors',color,'whiskers',0,'box',0,'boxwidth',20,'boxalpha',0.5,'box_xdata',parser.Results.box_xdata);
            xlabel(parser.Results.xlabel)
            ylabel(label_y)
        end




        function prop_mov = get_prop_mov(obj,prop,mov)
            prop_mov = obj.exp_mat(:,obj.insect_prop(prop),obj.mov_num(mov));
        end
    end
end