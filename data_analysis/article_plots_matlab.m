clear
close all
clc

path = 'I:\.shortcut-targets-by-id\1OA70vOJHDfV63DqG7LJCifTwW1h055ny\2024 Flight in the dark paper\data_exchange\'
fly_data = 'fly\all_data\'
mos_data = 'mosquito\'
gray_col_mat = colormap(gray)



file_name = '60ms_all_data.csv'
pert = {'40ms','80ms','100ms','60ms','step'}
for idx_file = 1:1:length(pert)
    fly.(['pert_',pert{idx_file}]) =  exp_class(readtable([path,fly_data,pert{idx_file},'_all_data']),16000);
    fly.fps = 16000
    fly.name = 'fly'
end

pert = {'60ms','100ms','step'}
for idx_file = 1:1:length(pert)
    mos.(['pert_',pert{idx_file}]) =  exp_class(load([path,mos_data,'mosquito_',pert{idx_file},'_shadow']).T,20000);
    mos.fps = 20000
    mos.name = 'mosquito'
end

plotter_obj = plotter()

%%

figure;
for i = 1:10:size(plotter_obj.col_mat,1)
yline(i,'color',plotter_obj.col_mat(i,:),'LineWidth',5);hold on
end
%% mean mov step

exp_name = 'pert_step'
prop_name = 'pitch'
label_y = 'pitch'
pert = 0
ax1 = subplot(2,1,1)

color_idx_fly = 220
color_idx_mov_fly = 240


color_idx_mos = 40
color_idx_mov_mos = 20
cluster_mean_color = [15,60]
cluster_mean_color_all_data = [30,100]

figure
max_time_xax = 250


insect = fly
mov = 13
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_fly,color_idx_mov_fly,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])

insect = mos
mov = 580
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_mos,color_idx_mov_mos,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])
%% mean mov 100ms pertubation

exp_name = {'pert_40ms','pert_80ms','pert_100ms','pert_60ms','pert_step'}
insect = fly
mov = 101
pert = 60
color_idx = 15
max_time_xax = 250
prop_name = 'forward_vel'
for i = 1:1:length(exp_name)
plotter_obj.plot_mean_subplot(fly,mos,exp_name{i},pert,color_idx + i*40,prop_name,max_time_xax);hold on
xlim([-20,250])
end
legend(exp_name)
%%
exp_name = {'pert_40ms','pert_60ms','pert_80ms','pert_100ms','pert_step'}
prop_name = 'forward_vel'
time = fly.(exp_name{1}).time_vec;


for i = 1:1:length(exp_name)


time_data = fly.(exp_name{i}).time_vec;

time_idx = find(time_data > -10 & time_data < 250);


prop_data = fly.(exp_name{i}).get_prop(prop_name);
mean_prop = mean(prop_data(time_idx,:),2,'omitnan');
mean_cell{i} = mean_prop;
end


mean_prop_mat = cell2mat(mean_cell);
tm = repmat(time_data(time_idx),length(mean_cell),1)';

surf(mean_prop_mat(1:50:end,:),tm(1:50:end,:))
%% mean mov 60ms pertubation

exp_name = 'pert_60ms'
pert = 60
insect = fly
mov = 70
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_fly,color_idx_mov_fly,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])

insect = mos
mov = 596
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_mos,color_idx_mov_mos,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])

%% mean mov 60ms pertubation

exp_name = 'pert_40ms'
pert = 40
insect = fly
mov = 24
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_fly,color_idx_mov_fly,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])

insect = mos
mov = 596
plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_idx_mos,color_idx_mov_mos,max_time_xax,cluster_mean_color,cluster_mean_color_all_data)
xlim([-20,250])
%% Violin delta angle
fly_col_idx = 190
mos_col_idx = 40

time_to_violin = linspace(-15,230,10)
exp_name = 'pert_60ms'
prop_name = 'vel_xy_ang_flat'
fly_color = plotter_obj.col_mat(fly_col_idx,:)
mos_color = plotter_obj.col_mat(mos_col_idx,:)


exp_name_cell = {'pert_60ms'}
pert = 60
figure
fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);


f_norm_vec = max([fly_fvec;mos_fvec]);

ax1 = subplot(2,1,2)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0,'box_xdata',0)
plotter_obj.pert_plot(pert,0,1,ax1)
ylabel('delta velocity angle [deg]')
title('60ms pertubation')
ylim([0,180])


time_to_violin = linspace(-15,230,10)
exp_name = 'pert_step'
prop_name = 'vel_xy_ang_flat'
exp_name_cell = {'pert_step'}
pert = 0


fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);

f_norm_vec = max([fly_fvec;mos_fvec]);

ax2 =subplot(2,1,1)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0)
plotter_obj.pert_plot(pert,0,1,ax2)
ylabel('delta velocity angle [deg]')
title('step')
ylim([0,180])
set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on

fly_col_idx = 190
mos_col_idx = 40

time_to_violin = linspace(-15,230,10)
exp_name = 'pert_40ms'
prop_name = 'vel_xy_ang_flat'
fly_color = plotter_obj.col_mat(fly_col_idx,:)
mos_color = plotter_obj.col_mat(mos_col_idx,:)


exp_name_cell = {'pert_40ms'}
pert = 40
figure
fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);


f_norm_vec = max([fly_fvec;mos_fvec]);

ax1 = subplot(1,1,1)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
plotter_obj.pert_plot(pert,0,1,ax1)
ylabel('delta velocity angle [deg]')
title('40ms pertubation')
ylim([0,180])
xlim([-50,300])
set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on


%% Violin delta angle
fly_col_idx = 190
mos_col_idx = 40

time_to_violin = linspace(-15,230,10)
exp_name = 'pert_60ms'
prop_name = 'cone_angle'
fly_color = plotter_obj.col_mat(fly_col_idx,:)
mos_color = plotter_obj.col_mat(mos_col_idx,:)


exp_name_cell = {'pert_60ms'}
pert = 60
figure
fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);


f_norm_vec = max([fly_fvec;mos_fvec]);

ax1 = subplot(2,1,2)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0,'box_xdata',0)
plotter_obj.pert_plot(pert,0,1,ax1)
ylabel('cone angle [deg]')
title('60ms pertubation')
ylim([0,180])


time_to_violin = linspace(-15,230,10)
exp_name = 'pert_step'
prop_name = 'cone_angle'
exp_name_cell = {'pert_step'}
pert = 0


fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);

f_norm_vec = max([fly_fvec;mos_fvec]);

ax2 =subplot(2,1,1)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0)
plotter_obj.pert_plot(pert,0,1,ax2)
ylabel('cone angle [deg]')
title('step')
ylim([0,180])
set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on




%%
figure
exp_name = 'pert_step'
prop = 'z_vel'
pert = 0
insect = mos
figure
% [data_a,data_b] = plotter_obj.cluster_plot(insect.('pert_60ms'),'z_vel',0,prop,'alpha',0.1);
[data_a,data_b] = plotter_obj.cluster_plot(insect.('pert_step'),'z_vel',0,prop,'alpha',0.1);

plotter_obj.mean_plot(insect.pert_60ms,data_a,prop,'plot_all_data',0,'color_idx',1);
plotter_obj.mean_plot(insect.pert_60ms,data_b,prop,'plot_all_data',0,'color_idx',200);

%%
props = {'cone_angle','z_vel','forward_vel'};

j = 1
    ini_time_mos = find(mos.pert_step.time_vec > 50,1)
    ini_time_fly = find(fly.pert_step.time_vec > 50,1)
    skip = 100
for prop = props

    mos_prop = mos.pert_step.get_prop(prop);
    mos_cell(j) =  {reshape(mos_prop(ini_time_mos:skip:end,:,:),1,[])};

    fly_prop = fly.pert_step.get_prop(prop);
    fly_cell(j) =  {reshape(fly_prop(ini_time_fly:skip:end,:,:),1,[])};
    j = j + 1;
end


% plot3(mos_cell{1},mos_cell{2},mos_cell{3},'.');hold on
% xlabel(props{1})
% ylabel(props{2})
% zlabel(props{3})
plot3(fly_cell{1},fly_cell{2},fly_cell{3},'.')
xlabel(props{1})
ylabel(props{2})
zlabel(props{3})