clear
close all
clc

path = 'G:\.shortcut-targets-by-id\1OA70vOJHDfV63DqG7LJCifTwW1h055ny\2024 Flight in the dark paper\data_exchange\'
fly_data = 'fly\all_data\'
mos_data = 'mosquito\'
gray_col_mat = colormap(gray)
col_mat = colormap(lines)


file_name = '60ms_all_data.csv'
pert = {'40ms','60ms','step'}
for idx_file = 1:1:length(pert)
    fly.(['pert_',pert{idx_file}]) =  exp_class(readtable([path,fly_data,pert{idx_file},'_all_data']));
    fly.fps = 16000
    fly.name = 'fly'
end

pert = {'60ms','step'}
for idx_file = 1:1:length(pert)
    mos.(['pert_',pert{idx_file}]) =  exp_class(load([path,mos_data,'mosquito_',pert{idx_file},'_shadow']).T);
    mos.fps = 20000
    mos.name = 'mosquito'
end

plotter = plotter()

%% model movie - step pertubation
exp_name = 'pert_step'
pert = 0
max_time_xax = 250
mov = 11
plotter.plot_pitch_for_vel_z_vel_model_mov(fly,exp_name,mov,pert,85,max_time_xax)

mov = 580
plotter.plot_pitch_for_vel_z_vel_model_mov(mos,exp_name,mov,pert,25,max_time_xax)
%% model movie - 60ms pertubation
exp_name = 'pert_60ms'
pert = 60
max_time_xax = 250

mov = 71
plotter.plot_pitch_for_vel_z_vel_model_mov(fly,exp_name,mov,pert,85,max_time_xax)

mov = 596
plotter.plot_pitch_for_vel_z_vel_model_mov(mos,exp_name,mov,pert,25,max_time_xax)
%% mean pitch
exp_name = 'pert_step'
prop_name = 'pitch'

plotter.plot_mean_subplot(fly,mos,exp_name,0,[85 25],prop_name,max_time_xax)

%%

exp_name = 'pert_step'
prop = 'beta'
pert = 0
figure
plotter.cluster_plot(mos.(exp_name),'z_vel',0,prop,'alpha',0.1)
plotter.mean_plot(mos.(exp_name),prop,prop,'plot_all_data',0)

%%
exp_name = 'pert_60ms'
pert = 60

figure
prop_name = 'vel_xy_ang_flat'
color = [1,0,0]
time_to_violin = linspace(-15,300,10)
ax = subplot(1,1,1)
plotter.violin_plot(mos.(exp_name),prop_name,time_to_violin,prop_name,color,mos.fps)
plotter.pert_plot(pert,1,1,ax)

exp_name = 'pert_step'
plotter.violin_plot(mos.(exp_name),prop_name,time_to_violin,prop_name,[0,1,0],mos.fps)
plotter.pert_plot(pert,1,1,ax)


%%
% function violin_plot(ylabel,color,varargin)
%     parser = inputParser;
%     addParameter(parser,'xlabel','time [ms]'); % number of camera pointing in Z lab axis
%     parse(parser, varargin{:})
% 
%     idx_to_violin = find(ismembertol(mos.(exp_name).time_vec,time_to_violin,1/fps -1/fps/2 ,'DataScale', 1))
%     x_names = mos.(exp_name).time_vec(idx_to_violin)
% 
%     prop = mos.(exp_name).get_prop(prop_name);
%     prop_time = prop(idx_to_violin,:);
%     h = daviolinplot(prop_time','xtlabels', x_names,'colors',color,'box',0,'whiskers',0,'violinalpha',0.5,'bins',15,'outliers',0);
%     xlabel(parser.Results.xlabel)
%     ylabel(ylabel)
% end

% exp_name = 'pert_step'
% exp_name = 'pert_60ms'

% fps = 20
% time_to_violin = linspace(0.79,200,10)
% idx_to_violin = find(ismembertol(fly.(exp_name).time_vec,time_to_violin,1/fps -1/fps/2 ,'DataScale', 1))
% x_names = fly.(exp_name).time_vec(idx_to_violin)
% 
% prop = fly.(exp_name).get_prop(prop_name);
% prop_time = prop(idx_to_violin,:);
% group_inx = prop_time*0+[1:1:length(idx_to_violin)]'
% h = daviolinplot(prop_time','xtlabels', x_names,'colors','r','box',0,'whiskers',0,'violinalpha',0.5,'bins',15,'outliers',0);
%%





%%


% mov(exp_obj.t0_idx,:)
% vel_xy_ang_flat = mos.(exp_name).get_prop('vel_xy_ang_flat');
% 
% vel_xy_ang_flat_500 = vel_xy_ang_flat(500,:)


% [h,L,MX,MED,bw]=violin(vel_xy_ang_flat_500)

% [experiment,insect_prop,mov_num] = generate_movie_mat(exp);
% pitch = get_prop(experiment,insect_prop,mov_num,'pitch');
% time = get_prop(experiment,insect_prop,mov_num,'time');
% t0_idx = find(time(:,:,1) == 0)
% mean_pitch = squeeze(mean(pitch,3,'omitnan'))
% 
% time = time(:,:,1);
% 
% prop2plot = 'pitch'
% 
% 
% 
% 
% 
% % hold on
% % plot(time_vec,squeeze(empty_movie(:,insect_prop(prop2plot),:)),'b');
% 
% model_movie = 73
% gray_col_mat = colormap(gray)
% col_mat = colormap(lines)
% 
% lw = 2
% sp(1)=subplot(2,1,1);hold on;box on
% yline(mean_pitch(t0_idx),'--k','LineWidth',3,'handlevisibility','off');
% plot(time,squeeze(pitch(:,1,:)),'LineWidth',3,'color',gray_col_mat(200,:))
% % plot(time_vec,empty_movie(:,insect_prop(prop2plot),mov_num(model_movie)),'LineWidth',3,'color',col_mat(13,:))
% plot(time,mean_pitch,'LineWidth',3,'color',col_mat(13,:))
% % h_lg=legend(prop2plot,'Location','southeast','edgecolor','none','Orientation','vertical','box','off');
% 
% set(gca, 'LineWidth', lw)
% set(gca,'fontsize',20)
% ylabel('Angle [°]')
% % ylim([20,75])
% 
% if pert == false
%    rec =  [0,sp(1).YLim(1),sp(1).XLim(2),sp(1).YLim(2)-sp(1).YLim(1)]
% else 
%     rec =  [0,sp(1).YLim(1),pert,sp(1).YLim(2)-sp(1).YLim(1)]
% end
% rectangle('Position',rec,...
%     'EdgeColor','none','FaceColor',[0,0,0,0.1]); % Plots the rectangle
% 



% % velocity
% sp(2)=subplot(2,1,2);
% hold on; box on
% plot(time_vec,all_body_vel{mov_ind}(:,3),'LineWidth',lw)
% plot(time_vec,all_vel_vel_x_body_xy{mov_ind},'LineWidth',lw)
% yline(0,'--k','LineWidth',3,'handlevisibility','off');
% h_lg=legend('z','forward','Location','southeast','edgecolor','none',...
%     'Orientation','vertical','box','off');
% set(gca, 'LineWidth', lw)
% set(gca,'fontsize',26)
% 
% % legend('z','forward')
% xlabel('time [ms]')
% ylabel('[m/s]')
% ylim([-0.2,0.25])
% 
% linkaxes(sp,'x')
% xlim([-23,200])
% rectangle('Position',[0,sp(2).YLim(1),sp(2).XLim(2),sp(2).YLim(2)-sp(2).YLim(1)],...
%     'EdgeColor','none','FaceColor',[0,0,0,0.1]); % Plots the rectangle
% % set(sp,'fontsize',40)





