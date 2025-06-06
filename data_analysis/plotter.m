classdef plotter
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        gray_mat
        col_mat
        pertubations
        path_to_save_fig
    end


    methods
        function obj = plotter(path_to_save_fig)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here

            obj.gray_mat = colormap(gray);
            obj.col_mat = colormap(turbo);
            obj.path_to_save_fig = path_to_save_fig;
        end
        function plot_prop(obj,insect,prop,ylabel,pert,mov_num,alpha,yticks,xlim_val,ylim_val,color_struct,...
                legend_val,name_of_figure,dir_to_save,cluster, position_cm,margin,location)
           
            
            
            fig = figure()
            hold on
            propip = insect.get_prop(prop);
            time = insect.time_vec;
                            ax = subplot(1,1,1);

            if cluster == false
                obj.all_data_plot(time,propip,'all_data_alpha',alpha,'all_data_color_idx',color_struct.all_data_color_idx);hold on
                obj.mean_plot(insect,propip,prop,'plot_all_data',0,'mean_color_idx',color_struct.mean_color_idx);
                obj.mov_plot(insect,prop,mov_num,prop,'color_idx',color_struct.color_idx_mov)
                obj.mean_plot(insect,prop,ylabel,'plot_all_data',0,'mean_color_idx',color_struct.mean_color_idx);
                obj.pert_plot(pert,0,1,ax(1));
            else

                obj.cluster_plot_mean(prop,insect,'colors',color_struct.cluster_mean_color,'color_all_data',color_struct.cluster_all_data_color,'alpha',0.1);
                obj.mov_plot(insect,prop,mov_num,ylabel,'color_idx',color_struct.color_idx_mov)
                obj.pert_plot(pert,0,1,ax);
            end
            if strcmp('z_vel',prop) == false 
                [response_fv,time_response] = insect.get_prop_point_in_time(prop,mov_num,'response_time')
                [break_fv,time_break] = insect.get_prop_point_in_time(prop,mov_num,'zero_v')
               
                scatter(time_response,response_fv,'^k','HandleVisibility','off','LineWidth',2)
                scatter(time_break,break_fv,'sk','HandleVisibility','off','LineWidth',2)

     
                text(time_break - location(1),break_fv+location(2),' BT')
                text(time_response - location(1),response_fv+location(2),'RT')
                
            end
            set([ax(1)], 'LineWidth', 1,'TickLength',[0.00,0.00]);box on
            set(gca,'fontsize',10);
            xlim(xlim_val)
            ylim(ylim_val)
        
            set(gca,'YTick',yticks)
            set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
            set(gca,'fontsize',10);
            set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
            % title(titletext,'Fontsize',10,'FontWeight','Normal')
            set(gca,'units','centimeters'...
            ,'position',position_cm);

            % Get the required axes position
            axPos = get(ax, 'Position');
            
            % Add margins to account for labels, ticks, etc.
            
            % Set the figure size
            figWidth = axPos(3) + 2 * margin(1);
            figHeight = axPos(4) + 2 * margin(2);
            set(gcf, 'Units', 'centimeters', 'Position', [10, 10, figWidth, figHeight]); % [x, y, width, height]
            
            % Optional: Improve appearance
            % set(gcf, 'PaperPositionMode', 'auto');
            path = [obj.path_to_save_fig,dir_to_save,name_of_figure]
            h = findall(fig,'-property','FontName');
            set(h,'FontName','Arial');  

            h2 = findall(fig,'-property','FontSize');
            set(h2,'FontSize',10);  
            legend(legend_val,'FontSize',8)

            set(fig,'renderer','painters')
            % print(fig,'-dsvg',path)

            print(fig, path, '-dsvg')  % Adjust DPI (e.g., -r150, -r200)

            % exportgraphics(gcf, [obj.path_to_save_fig,dir_to_save,name_of_figure], 'ContentType', 'vector');

    

        
            
        
        end

        function all_data_plot(obj,time,prop,varargin)
            parser = inputParser;
            addParameter(parser,'all_data_alpha',0.5); % opacity of the plot
            addParameter(parser,'all_data_color_idx',20); % color of the plot
            parse(parser, varargin{:})

            %plot(time,prop,'LineWidth',1.5,'color',[obj.gray_mat(200,:),parser.Results.alpha],'HandleVisibility','off')
            plot(time(1:100:end),prop(1:100:end,:),'LineWidth',0.5,'color',[obj.col_mat(parser.Results.all_data_color_idx,:),parser.Results.all_data_alpha],'HandleVisibility','off')
        end

        function cluster_plot_mean(obj,prop,insect,varargin)
            parser = inputParser;
            addParameter(parser,'colors',[180,85]); % number of camera pointing in Z lab axis
            addParameter(parser,'color_all_data',[180,85]); % number of camera pointing in Z lab axis
            addParameter(parser,'alpha',0.03); % number of camera pointing in Z lab axis
           parse(parser, varargin{:})
            [data_a,data_b] = obj.cluster_plot(insect,'z_vel',0,prop,parser.Results.color_all_data,'alpha',parser.Results.alpha);

            obj.mean_plot(insect,data_a,prop,'plot_all_data',0,'mean_color_idx',parser.Results.colors(1));
            obj.mean_plot(insect,data_b,prop,'plot_all_data',0,'mean_color_idx',parser.Results.colors(2));

        end


        function ax = mean_plot(obj,insect,prop,label_y,varargin)
            parser = inputParser;
            addParameter(parser,'mean_color_idx',13); % number of camera pointing in Z lab axis
            addParameter(parser,'plot_all_data',1); % number of camera pointing in Z lab axis
            addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
            addParameter(parser,'input_prop',0); % number of camera pointing in Z lab axis
            addParameter(parser,'initial_time_to_plot',-20); % number of camera pointing in Z lab axis

            parse(parser, varargin{:})
            time = insect.time_vec;


            if ischar(prop) == 1
                prop = insect.get_prop(prop);
            end
            idx_to_plot = find(time > parser.Results.initial_time_to_plot);

            mean_prop = mean(prop,2,'omitnan');
            ax = plot(time(idx_to_plot),mean_prop(idx_to_plot),'LineWidth',2,'color',obj.col_mat(parser.Results.mean_color_idx,:))
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
            plot(time(idx_to_plot),mov_prop(idx_to_plot),'LineWidth',parser.Results.LineWidth,'color',[obj.col_mat(parser.Results.color_idx,:),0.2],'HandleVisibility','off')
            xlabel(parser.Results.xlabel)
            ylabel(label_y)

        end


        function mov_mean_plot(obj,insect,prop_name,mov_num,label_y,varargin)
            parser = inputParser;
            addParameter(parser,'mean_color_idx',13); % color of mean data
            addParameter(parser,'color_idx_mov',50); % color of "model fly"
            addParameter(parser,'all_data_alpha',0.5); % opacity pf all data plot
            addParameter(parser,'all_data_color_idx',200); % index color of all data plot
            addParameter(parser,'xlabel','time [ms]'); % x label

            addParameter(parser,'plot_cluster',0); % plot/dont plot cluster (noam Vz) (1/0),  
            addParameter(parser,'cluster_mean_color',[25,85]); % color of mean of each cluster (assuming only 2) [cluster one index, cluster 2 index]
            addParameter(parser,'cluster_all_data_color',[25,85]);  % color of the clustered data (assuming only 2) [cluster one index, cluster 2 index]
            parse(parser, varargin{:})

            time = insect.time_vec;
            prop = insect.get_prop(prop_name);

            if parser.Results.plot_cluster == 1
                obj.cluster_plot_mean(prop_name,insect,'colors',parser.Results.cluster_mean_color,'color_all_data',parser.Results.cluster_all_data_color);
            else
                obj.all_data_plot(time,prop,'all_data_alpha',parser.Results.all_data_alpha,'all_data_color_idx',parser.Results.all_data_color_idx);hold on
                obj.mean_plot(insect,prop,label_y,'plot_all_data',0,'mean_color_idx',parser.Results.mean_color_idx);
            end
            obj.mov_plot(insect,prop_name,mov_num,label_y,'color_idx',parser.Results.color_idx_mov);
            xlabel(parser.Results.xlabel);
            ylabel(label_y);
            if parser.Results.plot_cluster == 1
                legend({'mean_A','mean_B','model'},'EdgeColor','None','Orientation','horizontal','Box','off','FontSize',8)
            else
                legend({'mean','model'},'EdgeColor','None','Orientation','horizontal','Box','off','FontSize',8)
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
            %rectangle('Position',rec,...
            %    'EdgeColor','none','FaceColor',[0,0,0,0.1],'FaceAlpha',0.1); % Plots the rectangle
            xregion(0, pert, "FaceColor",[0,0,0], "FaceAlpha",0.1,'handlevisibility','off') ; % tb
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

            plot(insect.time_vec(1:100:end),data_a(1:100:end,:),'color',[obj.col_mat(color_all_data(1),:),parser.Results.alpha],'LineWidth',0.5,'HandleVisibility','off');hold on
            plot(insect.time_vec(1:100:end),data_b(1:100:end,:),'color',[obj.col_mat(color_all_data(2),:),parser.Results.alpha],'LineWidth',0.5,'HandleVisibility','off');hold on

        end


        function plot_pitch_for_vel_z_vel_model_mov_mean(obj,insect,exp_name,mov,pert,color_struct,max_time_xax)
            mos_flag = 0;
            figure();
            ax1 = subplot(3,1,1)
            obj.mov_mean_plot(insect.(exp_name),'pitch',mov,'pitch [deg]',color_struct)

            
            obj.pert_plot(pert,0,1,ax1)


            ax2 = subplot(3,1,2)
            obj.mov_mean_plot(insect.(exp_name),'forward_vel',mov,'V_f_w_d [m/s]',color_struct)
            obj.pert_plot(pert,0,1,ax2)

            ax3 = subplot(3,1,3)
            if strcmp(insect.name,'mosquito')
                mos_flag = 1;
            end
            obj.mov_mean_plot(insect.(exp_name),'z_vel',mov,'V_z [m/s]',color_struct,'plot_cluster',mos_flag)
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

        function plot_pitch_for_vel_mean(obj,insect,exp_name,pert,max_time_xax,varargin)

            parser = inputParser;
            addParameter(parser,'mean_color_idx',13); % color of mean data
            addParameter(parser,'color_idx_mov',50); % color of "model fly"
            addParameter(parser,'all_data_alpha',0.5); % opacity pf all data plot
            addParameter(parser,'all_data_color_idx',200); % index color of all data plot
            addParameter(parser,'xlabel','time [ms]'); % x label

            addParameter(parser,'plot_cluster',0); % plot/dont plot cluster (noam Vz) (1/0),  
            addParameter(parser,'cluster_mean_color',[25,85]); % color of mean of each cluster (assuming only 2) [cluster one index, cluster 2 index]
            addParameter(parser,'cluster_all_data_color',[25,85]);  % color of the clustered data (assuming only 2) [cluster one index, cluster 2 index]
            parse(parser, varargin{:})
            figure();

            time = insect.(exp_name).time_vec;
            prop = insect.(exp_name).get_prop('pitch');

            ax1 = subplot(2,1,1)
            obj.all_data_plot(time,prop,'all_data_alpha',parser.Results.all_data_alpha,'all_data_color_idx',parser.Results.all_data_color_idx);hold on
            obj.mean_plot(insect.(exp_name),prop,'pitch [deg]','plot_all_data',0,'mean_color_idx',parser.Results.mean_color_idx);
            obj.pert_plot(pert,0,1,ax1);


            prop = insect.(exp_name).get_prop('forward_vel');

            ax2 = subplot(2,1,2)
            obj.all_data_plot(time,prop,'all_data_alpha',parser.Results.all_data_alpha,'all_data_color_idx',parser.Results.all_data_color_idx);hold on
            obj.mean_plot(insect.(exp_name),prop,'V_f_w_d [m/s]','plot_all_data',0,'mean_color_idx',parser.Results.mean_color_idx);
            obj.pert_plot(pert,0,1,ax1);



            ax1.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            ax2.XLim = [ax1.XLim(1),min(ax1.XLim(2),max_time_xax)]
            linkaxes([ax1,ax2],'x')
            pert = split(exp_name,'_')
            sgtitle([insect.name,' ',pert{2}])


            set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on
            % set(gca,'fontsize',20)

        end
        
        function bar_plot(obj,data_to_plot,header,color,err,varargin)

            parser = inputParser;
            addParameter(parser,'plot_err',true); % color of mean data
            parse(parser, varargin{:})

            for k = 1:1:length(data_to_plot)               
            pertubation = split(header{k}, '_');
            pertubation_cell{k} = pertubation{end};
            bar(k,data_to_plot(k),'FaceColor',color(k,:),'LineWidth',1);hold on
            if parser.Results.plot_err == true
            e = errorbar(k,data_to_plot(k),err(k))
            e.Color = color(k,:)
            end



            end
            
            % % Bar plot formatting
            % set(gca,'units','centimeters'...
            % ,'position',position_cm);
            % 
            % % Get the required axes position
            % axPos = get(ax, 'Position');
            % 
            % % Add margins to account for labels, ticks, etc.
            % 
            % % Set the figure size
            % figWidth = axPos(3) + 2 * margin(1);
            % figHeight = axPos(4) + 2 * margin(2);
            % set(gcf, 'Units', 'centimeters', 'Position', [10, 10, figWidth, figHeight]); % [x, y, width, height]
            % 
            % pertubation_cell{end} = 'Step'
            % % ylabel(ylabel_in);
            % xticks(1:length(data_to_plot));
            % xticklabels(pertubation_cell);
            % xlabel('Dark pulse duration');
            % set(gca,'fontsize',12,'LineWidth',1);
            % if path_to_save ~= false
            %     path = [path_to_save]
            %     h = findall(fig,'-property','FontName');
            %     set(h,'FontName','San Serif');
            %     print(fig,'-dsvg',path)
            % end
        end



        function histogram_plot(obj,data_to_plot,color,path_to_save,position_cm,margin,ttl)

            fig = figure()
            ax = subplot(1,1,1);
            histogram(data_to_plot,'FaceColor',color,'FaceAlpha',1,'LineWidth',2,BinWidth=4);

            
            % Bar plot formatting
            set(gca,'units','centimeters'...
            ,'position',position_cm);

            % Get the required axes position
            axPos = get(ax, 'Position');
            
            % Add margins to account for labels, ticks, etc.
            
            % Set the figure size
            figWidth = axPos(3) + 2 * margin(1);
            figHeight = axPos(4) + 2 * margin(2);
            set(gcf, 'Units', 'centimeters', 'Position', [10, 10, figWidth, figHeight]); % [x, y, width, height]
            

            xlabel('Time [ms]');
            ylabel('Counts');
            title(ttl)

            set(gca,'fontsize',12,'LineWidth',1);
            if path_to_save ~= false
                path = [path_to_save]
                h = findall(fig,'-property','FontName');
                set(h,'FontName','Arial');
                set(fig,'renderer','painters')
                print(fig,'-dsvg',path)
            end
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