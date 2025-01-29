clear
close all
clc

SAVE_FIGS = false ;

% Load data
path = 'G:\.shortcut-targets-by-id\1OA70vOJHDfV63DqG7LJCifTwW1h055ny\2024 Flight in the dark paper\data_exchange\'
path_for_figs = 'G:\.shortcut-targets-by-id\1OA70vOJHDfV63DqG7LJCifTwW1h055ny\2024 Flight in the dark paper\media\'
fly_data = 'fly\all_data\'
mos_data = 'mosquito\'

file_name = '60ms_all_data.csv'

pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step','darkan'} % Fly: pertubations to load (for all - {'40ms','80ms','100ms','60ms','step'})
for idx_file = 1:1:length(pert)
    fly.(['pert_',pert{idx_file}]) =  exp_class(readtable([path,fly_data,pert{idx_file},'_all_data']),16000,'fly');
    fly.fps = 16000;
    fly.name = 'fly';
end
pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step','darkan'} % Fly: pertubations to load (for all - {'40ms','80ms','100ms','60ms','step'})

props = [{'min_v'},{'zero_v'},{'response_time'},{'delta_angle'},{'zero_v_z'},{'min_v_time'}]
for idx_file = 1:1:length(pert)
for k = 1:1:length(props)
    path_pulses = ['\fly\pulses\',props{k},'.csv']
    data = readtable([path  ,path_pulses]);
    header = data.Properties.VariableNames(2:end);
    pertubation = ['pert_',pert{idx_file}];
    movs = cell(data.Var1);

    idx_to_keep = isnan(data.(pertubation)) == false & data.(pertubation) ~= 999
    fly.(pertubation).(props{k}) = table(data.(pertubation)(idx_to_keep));
    fly.(pertubation).(props{k}).Properties.RowNames = movs(idx_to_keep);
end
end


pert = {'60ms','100ms','step'} ;% Mosquito: pertubations to load
for idx_file = 1:1:length(pert)
    mos.(['pert_',pert{idx_file}]) =  exp_class(load([path,mos_data,'mosquito_',pert{idx_file},'_shadow']).T,20000,'mosquito');
    mos.fps = 20000;
    mos.name = 'mosquito';
end

mos.pert_60ms.response_time = table(48.6)
mos.pert_60ms.response_time.Properties.RowNames = {'mov395'}
mos.pert_60ms.zero_v = table(NaN)
mos.pert_60ms.zero_v.Properties.RowNames = {'mov395'}


mos.pert_step.response_time = table(36.6)
mos.pert_step.response_time.Properties.RowNames = {'mov169'}
mos.pert_step.zero_v = table(82)
mos.pert_step.zero_v.Properties.RowNames = {'mov169'}

%% Fly leg spreading times
time = fly.('pert_step').time_vec;


pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step','darkan'}
for k = 1:1:length(pert)
pertubation = sprintf('pert_%s',pert{k});
file_name = ['leg_open_',pert{k}];
open_leg_mat = readmatrix(['H:\My Drive\dark 2022\excel\',file_name]);
open_leg_mat = open_leg_mat(:,2)
idx_to_keep = isnan(open_leg_mat) == false & open_leg_mat~= 999% & open_leg_mat~= 1
fly.(pertubation).open_leg = table(open_leg_mat(idx_to_keep));
end

plotter_obj = plotter(path_for_figs);



%% a plot to help decide which colors to choose for the other plots (it plots the color matrix)

figure ;
for i = 1:10:size(plotter_obj.col_mat,1)
    yline(i,'color',plotter_obj.col_mat(i,:),'LineWidth',5);hold on
end
%% set colors

% maximum time (x axis) [ms]
max_time_xax = 250;
% --------------------------------------------------------------

% ---- colors (defined by the indices in color mat) --------------

% color struct - fly:----------------------------------------

color_struct_fly.mean_color_idx = 240 ;% mean of data color
color_struct_fly.color_idx_mov  = 240 ;% "model movie" color
color_struct_fly.all_data_alpha = 0.3 ;% all data opacity
color_struct_fly.all_data_color_idx = 200; % all data color

% color struct - mosquito:----------------------------------------
color_struct_mos.mean_color_idx = 15 ;% mean of data color
color_struct_mos.color_idx_mov  = 15 ;% "model movie" color
color_struct_mos.all_data_alpha = 0.15 ;% all data opacity
color_struct_mos.all_data_color_idx = 40 ;% all data color

color_struct_mos.cluster_mean_color = [15,60] ;% indices of colors for the clusters (clustered by Vz) [idx1, idx2]
color_struct_mos.cluster_all_data_color = [30,100] ;% indices of colors for the mean of clusters (clustered by Vz) [idx1, idx2]
% --------------------------------------------------------------------


%%

% 
% 
% 
% property_to_plot = 'forward_vel' ; % you can choose any property from this list: fly.pert_step.insect_prop
% property = fly.('pert_step').get_prop(property_to_plot);
% %time_mov = open_leg_mat(idx_of_movie,:);
% time_mov = open_leg_mat(:,1);
% 
% 
% 
% for k = 1:1:length(open_leg_mat(:,1))
%     try
% 
%         [v,i] =min(abs(open_leg_mat(k,1)- time));
%         mov_num = fly.pert_step.mov_num(open_leg_mat(k,2));
% 
%         plot(time,property(:,mov_num));hold on
%         scatter(time(i),property(i,mov_num),'*k');
%         property_legs_open(k) = property(i,mov_num);
% 
%     catch
%         continue
%     end
% end




%% leg spreading histogram

for k = 1:1:length(pert)
pertubation = sprintf('pert_%s',pert{k})
cell_openleg{k} = fly.(pertubation).open_leg.Var1

end
all_open_leg_time = cell2mat(cell_openleg');
all_open_leg_time= all_open_leg_time(all_open_leg_time > 0);

red = [0.8549    0.2078    0.0039];
blue = [0.2706    0.4863    0.9569];
path_to_save_fig = [path_for_figs,'\appendix\leg_spread.svg']
position_cm = [2,2 8 6]
margin = [2,2]
all_open_leg_time = all_open_leg_time(all_open_leg_time > 1)
mu = mean(all_open_leg_time(:,1))
st = std(all_open_leg_time(:,1))
ttl = ({"Fly leg spreading times \mu=" + round(mu) + "\pm" + round(st) + " ms"; " "})

plotter_obj.histogram_plot(all_open_leg_time,red,path_to_save_fig,position_cm,margin,ttl)





%% feature table

props = [{'response_time'},{'zero_v'},{'delta_angle'},{'open_leg'},{'zero_v_z'}];
pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step','darkan'};
path_to_save = [path_for_figs,'\figure5\','pulse_features.svg']


for idx_prop = 1:1:length(props)
    for k = 1:1:length(pert)
        if strcmp(props{idx_prop},'delta_angle')
            percent_feature.(props{idx_prop})(k) = fly.(['pert_',pert{k}]).get_percent_with_th(props{idx_prop},50);
        elseif strcmp(props{idx_prop},'open_leg')

            open_legs_perc = fly.(['pert_',pert{k}]).get_legs_percent(props{idx_prop});
            percent_feature.legs_nocrazy(k) = open_legs_perc(1);
            percent_feature.legs_crazy(k) = open_legs_perc(2);
            percent_feature.legs_only_crazy(k) = open_legs_perc(3);

        else
            percent_feature.(props{idx_prop})(k) = fly.(['pert_',pert{k}]).get_percent(props{idx_prop});

        end
     end
end


%% features histogram
path_to_save = [path_for_figs,'\figure5\','pulse_features.svg']
pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step'};
props = [{'response_time'},{'legs_crazy'},{'zero_v'},{'delta_angle'}];
title_cell = {'Response %','Spread legs','V_{fwd} = 0','|\Delta \gamma| < 50 [deg]'}
pert_xlable = {'5','10','20','40','60','80','100','step'};
props_color = [0.2, 0.7, 0.2;0.2, 0.1, 0.9;0.9, 0.1, 0.05;0.5, 0.1, 0.5;0.3, 0.1, 0.05]


gray =[1,0.9,0.8,0.7,0.6,0.4,0.3,0]
gray_cmap = repmat(gray,3,1)'
cmap = gray_cmap(1:length(pert),:) ;


fig = figure
t = tiledlayout("vertical","TileSpacing","compact");
ylocation  = [0.82,0.67,0.525,0.39,0.245]
xlocation = [0.2,0.2,0.2,0.63,0.63]



for k = 1:1:length(props)
tl(k) = nexttile
color = cmap .* props_color(k,:) ;
plotter_obj.bar_plot(percent_feature.(props{k})(1:end-1),pert,color,percent_feature.(props{k})(1:end-1)*0,'plot_err',false)
set(gca,'XTick',[])
ylim([0,100])
tl(k).YAxis.FontSize = 10; % Set the Y-axis tick label font size
annotation('textbox', [xlocation(k), ylocation(k), 0.3, 0.11], ...
    'String', title_cell{k}, ...
    'FontSize', 8,'LineStyle','none' );
end


for k = 1:1:length(pert)
    min_v = fly.(['pert_',pert{k}]).min_v.Var1;
    mean_min_v(k) = fly.(['pert_',pert{k}]).get_mean('min_v')
    min_v = min_v(min_v < 1000);
    std_minv(k) = std(min_v)
end


annotation('textbox', [xlocation(5), ylocation(5), 0.3, 0.11], ...
    'String', 'Minimum velocity', ...
    'FontSize', 8,'LineStyle','none' );
set(gca,'XTick',[])

lastTile  = nexttile([2,1])
color = cmap .* props_color(end,:) ;
plotter_obj.bar_plot(mean_min_v(1:end),pert,color,std_minv)
xticks(1:length(pert));
xticklabels(pert_xlable);


xlabel('Dark pulse duration [ms]','FontSize',10);
ylabel(lastTile,'V_{min} [m/s]','FontSize',10)
ylabel(tl(3),'percent [%]','FontSize',10)
tl(3).YLabel.Position = [tl(3).YLabel.Position(1)-0.5,120,-1]

set(gcf, 'Units', 'centimeters', 'Position', [1, 1, 9, 13]); % [x, y, width, height]
h = findall(fig,'-property','FontName');
set(h,'FontName','San Serif');
    print(fig,'-dsvg',path_to_save)


% plotter_obj.bar_plot(mean_min_v,pert,color,path_to_save_fig,std_minv,position_cm,margin)
%% get mean and std of properties for the timeliene plot
pert = {'5ms','10ms','20ms','40ms','60ms','80ms','100ms','step'};
props = [{'response_time'},{'open_leg'},{'zero_v'},{'delta_angle'},{'min_v_time'}];
for k = 1:1: length(pert)
    for idx_prop = 1:1:length(props)
        if strcmp('open_leg',props{idx_prop}) == true
            [mean_ol,std_ol]= fly.(['pert_',pert{k}]).get_mean_std_open_legs(props{idx_prop})
            mean_prop.(props{idx_prop})(k) = mean_ol
            std_prop.(props{idx_prop})(k) = std_ol
        else
            idx = fly.(['pert_',pert{k}]).(props{idx_prop}).Var1 ~= 1000
            mean_prop.(props{idx_prop})(k)= mean(fly.(['pert_',pert{k}]).(props{idx_prop}).Var1(idx))
            std_prop.(props{idx_prop})(k)  = std(fly.(['pert_',pert{k}]).(props{idx_prop}).Var1(idx))
        end
    end
end
%%
fig = figure();
ax = subplot(1,1,1);
path_to_save = [path_for_figs, '\figure5\', 'manouver_timeline.svg'];
xlim_val = [-25 250];
ylim_val = [-0.2 0.3];
yticks = [-0.2, 0, 0.3];
margin = [2, 1];
position_cm = [2, 1.5, 6.7, 3];
time = fly.pert_60ms.time_vec;
prop = fly.pert_60ms.get_prop('forward_vel');
label_y = {'V_{fwd} [m/s]'}
prop = 'forward_vel'

% Plot the data
plotter_obj.mean_plot(fly.pert_60ms,prop,label_y,'mean_color_idx',color_struct_fly.mean_color_idx)
plotter_obj.pert_plot(60, 0, 1, ax(1));

% Add vertical lines with specific properties
xline(mean_prop.open_leg(5) , 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.2, 0.1, 0.9]);
xline(mean_prop.response_time(5), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.2, 0.7, 0.2]);
xline(mean_prop.zero_v(5), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.9, 0.1, 0.05]);
xline(mean_prop.min_v_time(5), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.3, 0.1, 0.05]);

xline(200, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.5, 0.1, 0.5]);

% Set axis limits and labels
xlim(xlim_val);
ylim(ylim_val);
xlabel('time [ms]', 'FontSize', 10);
ylabel('V_{fwd} [m/s]', 'FontSize', 10);

text_features = {'Response','Legs spreading','Break','V_{min}','Re-orient velocity'}
text_location = [0.36, 0.29, 0.18, 0.11;0.40, 0.29, 0.3, 0.11;0.44, 0.29, 0.12, 0.11;0.53, 0.3, 0.14, 0.15;0.76, 0.29, 0.14, 0.25]
for k = 1:1:length(props)
% Add annotations with explicit font size
annotation('textbox', text_location(k,:), ...
    'String', text_features{k}, ...
    'Color', props_color(k,:), ...
    'FontSize', 8,'rotation',90,'LineStyle','none' );

end


% Ensure all text elements have the correct font size
h = findall(fig, '-property', 'FontSize');
set(h, 'FontSize', 10);

% Set axis properties and position
set(ax, 'LineWidth', 0.5, 'TickLength', [0.00, 0.00]);
box on;
set(gca, 'YTick', yticks);
set(gca, 'TickDir', 'out', 'TickLength', [0.01, 0.025]);
set(gca, 'Units', 'centimeters', 'Position', position_cm);

% Adjust figure size to include margins
axPos = get(ax, 'Position');
figWidth = axPos(3) + 2 * margin(1);
figHeight = axPos(4) + 2 * margin(2);
set(gcf, 'Units', 'centimeters', 'Position', [1, 1, figWidth, figHeight]);

% Set font for all text objects explicitly
h = findall(fig, '-property', 'FontName');
set(h, 'FontName', 'San Serif');


% Save the figure
print(fig, '-dsvg', path_to_save);


%% fly antenna
props = [{'response_time'},{'zero_v'},{'delta_angle'},{'legs_crazy'},{'zero_v_z'}];
% path_to_save = [path_for_figs,'\figure5\','pulse_features_anntena.svg']


step_and_100_props.response_time = [fly.pert_100ms.response_time.Var1;fly.pert_step.response_time.Var1];
step_and_100_props.zero_v = [fly.pert_100ms.zero_v.Var1;fly.pert_step.zero_v.Var1];
step_and_100_props.open_leg = [fly.pert_100ms.open_leg;fly.pert_step.open_leg];
step_and_100_props.zero_v_z = [fly.pert_100ms.zero_v_z.Var1;fly.pert_step.zero_v_z.Var1];

idx_open_all_time = (step_and_100_props.open_leg.Var1 == 1);
idx_not_opening = (step_and_100_props.open_leg.Var1 == 0);

step_100.response_time  = [sum(step_and_100_props.response_time ~= 1000),length(step_and_100_props.response_time) - sum(step_and_100_props.response_time ~= 1000)]
step_100.zero_v  = [sum(step_and_100_props.zero_v ~= 1000),length(step_and_100_props.zero_v) - sum(step_and_100_props.zero_v ~= 1000)]
step_100.open_leg  = [sum(idx_not_opening == 0 & idx_open_all_time == 0),length(step_and_100_props.open_leg.Var1) - sum(idx_not_opening == 0 & idx_open_all_time == 0)]
step_100.zero_v_z  = [sum(step_and_100_props.zero_v_z ~= 1000),length(step_and_100_props.zero_v_z) - sum(step_and_100_props.zero_v_z ~= 1000)]

idx_open_all_time = (fly.pert_darkan.open_leg.Var1 == 1);
idx_not_opening = (fly.pert_darkan.open_leg.Var1 == 0);


antenna.response_time = [sum(fly.pert_darkan.response_time.Var1~= 1000),length(fly.pert_darkan.response_time.Var1) - sum(fly.pert_darkan.response_time.Var1~= 1000)]
antenna.zero_v = [sum(fly.pert_darkan.zero_v.Var1~= 1000),length(fly.pert_darkan.zero_v.Var1) - sum(fly.pert_darkan.zero_v.Var1~= 1000)]
antenna.open_leg  = [sum(idx_not_opening == 0 & idx_open_all_time == 0),length(fly.pert_darkan.open_leg.Var1) - sum(idx_not_opening == 0 & idx_open_all_time == 0)]
antenna.zero_v_z = [sum(fly.pert_darkan.zero_v_z.Var1~= 1000),length(fly.pert_darkan.zero_v_z.Var1) - sum(fly.pert_darkan.zero_v_z.Var1~= 1000)]

% for k = 1:1:4
%     num_insect(k) = sum([step_100.(prop{k})',antenna.(prop{k})'],'all')
% 
% end

%%
fig = figure()
title_prop = {'Response %','V_{fwd}=0','spread leg','V_{Z}=0'}
prop = {'response_time','zero_v','open_leg','zero_v_z'}
for k = 1:1:length(prop)
subplot(1,length(prop),k)
fisher_mat = [step_100.(prop{k});antenna.(prop{k})]
fisher_table = table(fisher_mat(:,1),fisher_mat(:,2),'VariableNames',{'respond','do_not_respond'},'RowNames',{'intact','not_intact'})
    [h,p,stats] = fishertest(fisher_table)
    fisher.(prop{k}) = p
end
% for k_row = 1:1:2
%     ba = bar(k_row,[fisher_mat(k_row,1);fisher_mat(k_row,2)],'stacked');hold on
%     set(ba, 'FaceColor', 'Flat')
%     ba(1).CData = [0.4660    0.6740    0.1880];
%     ba(2).CData = [0.900    0.1250    0.0980]
%     xticks(1)
% 
% end
% 
% xticks(1:2)
% xticklabels({'intact','detached'});
% 
% title(title_prop{k})
% pval = sprintf('p_{val}=%.2e',p)
% text(0,max(num_insect) - 20,pval)
% ylim([0,max(num_insect) + 10])
% 
% end

% ylabel(lastTile,'V_{min} [m/s]')
% ylabel(tl(3),'percent [%]')
%%
for k = 1:1:4
    num_insect(k) = sum([step_100.(prop{k})',antenna.(prop{k})'],'all')

end
fig = figure()
title_prop = {'Response %','V_{fwd}=0','spread leg','V_{Z}=0'}
prop = {'response_time','zero_v','open_leg','zero_v_z'}
for k = 1:1:length(prop)
    subplot(1,length(prop),k)
    fisher_mat = [step_100.(prop{k});antenna.(prop{k})]
fisher_table = table(fisher_mat(:,1),fisher_mat(:,2),'VariableNames',{'respond','do_not_respond'},'RowNames',{'intact','not_intact'})
[h,p,stats] = fishertest(fisher_table)
ba = bar(1,100*[fisher_mat(1,1)/sum(fisher_mat(1,:)),1 - fisher_mat(1,1)/sum(fisher_mat(1,:))],'stacked');hold on
set(ba, 'FaceColor', 'Flat')

ba(1).CData = [0.4660    0.6740    0.1880];
ba(2).CData = [0.900    0.1250    0.0980]
xticks(1)


ba2 = bar(2,100*[fisher_mat(2,1)/sum(fisher_mat(2,:)),1 - fisher_mat(2,1)/sum(fisher_mat(2,:))],'stacked');hold on
set(ba2, 'FaceColor', 'Flat')

ba2(1).CData = [0.4660    0.6740    0.3880];
ba2(2).CData = [0.900    0.1250    0.2980]
xticks(1:2)
xticklabels({'intact','detached'});

title(title_prop{k})
pval = sprintf('p_{val}=%.2e',p)
text(0,130 - 20,pval)
ylim([0,100 + 30])

end

% ylabel(lastTile,'V_{min} [m/s]')
% ylabel(tl(3),'percent [%]')






%% fly step
W = 425 ;
H = 350 ;
mov = 15 %71 21
position_cm = [1.5,1.5 5 4]
margin = [1.5, 1]; % Example margins [width, height] in centimeters


location_text = [35,-0.05]
plotter_obj.plot_prop(fly.pert_step,'forward_vel','V_f_w_d [m/s]',0,mov,0.3,[-0.3, 0 0.3],[-25 250],[-0.3,0.3],...
    color_struct_fly,{'Mean','movie'},'fly_vfwd_step.svg','\figure2\',false,position_cm,margin,location_text)


location_text = [20,-5]
position_cm = [1.5,1.5 5 4]
plotter_obj.plot_prop(fly.pert_step,'pitch','pitch [deg]',0,mov,0.3,[20,50,80],[-25 250],[20,80],...
    color_struct_fly,{'Mean','Sample'},'fly_pitch_step.svg','\figure2\',false,position_cm,margin,location_text)

location_text = [0,0]
position_cm= [1.8,1.5 4.5 3.5]
plotter_obj.plot_prop(fly.pert_step,'z_vel','V_z [m/s]',0,mov,0.3,[-0.2,-0.1, 0 0.1],[-25 250],[-0.3 0.2],...
    color_struct_fly,{'Mean','Sample'},'fly_vz_step.svg','\figure3\',false,position_cm,margin,location_text)

%% mosquito step
W = 425 ;
H = 350 ;
mov = 169
location_text = [35,-0.06]
position_cm = [1.5,1.5 5 4]
plotter_obj.plot_prop(mos.pert_step,'forward_vel','V_f_w_d [m/s]',0,mov,0.3,[-0.3, 0 0.3],[-25 250],[-0.3,0.3],...
    color_struct_mos,{'Mean','Sample'},'mos_vfwd_step.svg','\figure2\',false,position_cm,margin,location_text)

location_text = [20,7]
position_cm = [1.5,1.5 5 4]
plotter_obj.plot_prop(mos.pert_step,'pitch','pitch [deg]',0,mov,0.3,[20,50,80],[-25 250],[20,80],...
    color_struct_mos,{'Mean','Sample'},'mos_pitch_step.svg','\figure2\',false,position_cm,margin,location_text)

color_struct_mos.cluster_all_data_color = [30,80] ;% indices of colors for the mean of clusters (clustered by Vz) [idx1, idx2]
color_struct_mos.cluster_mean_color = [15,50]
position_cm= [1.8,1.5 4.5 3.5]
margin = [1.3,1]
plotter_obj.plot_prop(mos.pert_step,'z_vel','V_z [m/s]',0,mov,0.3,[-0.25,-0.1, 0 0.1],[-25 250],[-0.25 0.2],...
    color_struct_mos,{'Mean a','Mean b','Sample'},'mos_vz_step.svg','\figure3\',true,position_cm,margin)


%% fly 60
W = 425 ;
H = 350 ;
mov = 80
position_cm= [1.5,1.5 5 4]
location_text = [10,0.06]

% plotter_obj.plot_prop(fly.pert_60ms,'z_vel','V_z [m/s]',60,mov,0.3,[-0.2, 0 0.3],[-25 250],[-0.2,0.3],...
%     color_struct_fly,{'mean','movie'},'fly_vz_60ms.svg','\figure4\',false,position_cm)
plotter_obj.plot_prop(fly.pert_60ms,'forward_vel','V_f_w_d [m/s]',60,mov,0.3,[-0.2, 0 0.3],[-25 250],[-0.2,0.3],...
    color_struct_fly,{'mean','Sample'},'fly_vfwd_60ms.svg','\figure4\',false,position_cm,margin,location_text)

location_text = [30,5]
position_cm = [1.5,1.5 5 4]
plotter_obj.plot_prop(fly.pert_60ms,'pitch','pitch [deg]',60,mov,0.3,[20,50,80],[-25 250],[20,80],...
    color_struct_fly,{'mean','Sample'},'fly_pitch_60ms.svg','\figure4\',false,position_cm,margin,location_text)

%% mosquito 60
W = 425 ;
H = 350 ;
mov = 395;
position_cm= [1.5,1.5 5 4]
location_text = [10,0.06]

plotter_obj.plot_prop(mos.pert_60ms,'forward_vel','V_f_w_d [m/s]',60,mov,0.3,[-0.2, 0 0.3],[-25 250],[-0.2,0.3],...
    color_struct_mos,{'mean','Sample'},'mos_vfwd_60ms.svg','\figure4\',false,position_cm,margin,location_text)

location_text = [30,5]

plotter_obj.plot_prop(mos.pert_60ms,'pitch','pitch [deg]',60,mov,0.3,[20,50,80],[-25 250],[20,80],...
    color_struct_mos,{'mean','Sample'},'mos_pitch_60ms.svg','\figure4\',false,position_cm,margin,location_text)


% plotter_obj.plot_prop(mos.pert_60ms,'z_vel','V_z [m/s]',60,mov,0.3,[W,H],[-0.2,-0.1, 0 0.1],[-25 250],[-0.3 0.2],...
%     color_struct_mos,{'mean','movie'},'60 ms','mos_vz_60ms.svg','\figure4\',true)



%% Violin delta angle
% ---- colors (defined by the indices in color mat) --------------
fly_col_idx = 190 % color of fly
mos_col_idx = 40 % color of mosquito
% time to plot violin (x axis, the initial point of the violin)
time_to_violin = linspace(-15,230,10)

%----- violin plot properties-----------
prop_name = 'vel_xy_ang_flat' % property, a list of all properies: fly.pert_60ms.insect_prop
fly_color = plotter_obj.col_mat(fly_col_idx,:) % color of the fly's violin
mos_color = plotter_obj.col_mat(mos_col_idx,:) % color of the mosquito's violin


% 60 ms pertubation -----------------------------
exp_name_cell = {'pert_60ms'}  % pertubation
exp_name = exp_name_cell{1}
pert = 60 % used to plot the pertubation as a gray box

% plot violin-------------
figure
fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);
f_norm_vec = max([fly_fvec;mos_fvec]); % normalize

%f_norm_vec(:) = max(max([fly_fvec;mos_fvec])) / 4 ;
fig = figure()
ax1 = subplot(1,1,1)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0,'box_xdata',0)
plotter_obj.pert_plot(pert,0,1,ax1)
ylabel('\Delta velocity angle ,\Delta \gamma [deg]')
ylim([0,180])
 set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
            ,'position',[6 5 20 6]);
set([ax1], 'LineWidth', 1,'TickLength',[0.00,0.00]);box on
legend('Location', 'northwest');

path = [plotter_obj.path_to_save_fig,'/figure4/','violin_60ms.svg']
h = findall(fig,'-property','FontName');
set(h,'FontName','San Serif');
print(fig,'-dsvg',path)


% step pertubation ------------
exp_name_cell = {'pert_step'}
exp_name = exp_name_cell{1}
pert = 0 % used to plot the pertubation as a gray box


fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);
f_norm_vec = max([fly_fvec;mos_fvec]);

fig = figure()

ax2 =subplot(1,1,1)
plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0)
plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0)
plotter_obj.pert_plot(pert,0,1,ax2)
ylabel('\Delta velocity angle ,\Delta \gamma [deg]')
ylim([0,180])
set([ax2], 'LineWidth', 1,'TickLength',[0.00,0.00]);box on
 set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
            ,'position',[6 5 20 6]);

legend('Location', 'northwest');
path = [plotter_obj.path_to_save_fig,'/figure4/','violin_step.svg']
h = findall(fig,'-property','FontName');
set(h,'FontName','San Serif');
print(fig,'-dsvg',path)





% %% Tsevi plot (THIS)
% exp_name = 'pert_step' % name of pertubation
% pert = 0 % % pertubation (step - 0, 100 ms - 100, 60 ms - 60...)
% subplot(2,1,1)
% prop = [{'forward_vel'},{'forward_vel'}]
% insect = fly % insect to plot
% propip = fly.(exp_name).get_prop(prop{subplot_idx});
% time = fly.(exp_name).time_vec
% 
% 
% plotter_obj.all_data_plot(time,propip,'all_data_alpha',color_struct_fly.all_data_alpha,'all_data_color_idx',color_struct_fly.all_data_color_idx,location_text);hold on
% plotter_obj.mean_plot(insect.(exp_name),propip,prop{subplot_idx},'plot_all_data',0,'mean_color_idx',color_struct_fly.mean_color_idx,location_text);
% plotter_obj.mov_plot(insect.(exp_name),prop{subplot_idx},13,prop{subplot_idx},'color_idx',color_struct_fly.color_idx_mov,location_text)
% 
% W = 425 ;
% H = 535 ;
% 
% set(gcf,'position',[200 50 W H])
% subplot(2,1,1)
% 
% ylim([0,100]) ; set(gca,'fontsize',14);
% xlim([-25 250])
% set(gca,'YTick',[30,50,70])
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% xlabel('') ;
% sgtitle('Fruit fly','Fontsize',16)
% %%
% subplot(2,1,1)
% ylim([-0.3,0.3]) ;
% set(gca,'ytick',[-0.3 0 0.3]) ;
% xlim([-25 250])
% set(gca,'fontsize',14);
% set(gca,'TickDir','out');
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% if SAVE_FIGS
%     print(gcf,'fly step response.png','-dpng','-r300')
% end
% 
% 
% insect = mos % insect to plot
% plotter_obj.plot_pitch_for_vel_mean(insect,exp_name,pert,max_time_xax,color_struct_mos)
% 
% set(gcf,'position',[200 50 W H])
% subplot(2,1,1)
% ylim([25,70]) ; set(gca,'fontsize',14);
% xlim([-25 250])
% set(gca,'YTick',[30,50,70])
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% xlabel('') ;
% sgtitle('Mosquito','Fontsize',16)
% 
% subplot(2,1,2)
% ylim([-0.3,0.3]) ;
% set(gca,'ytick',[-0.3 0 0.3]) ;
% xlim([-25 250])
% set(gca,'fontsize',14);
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% 
% if SAVE_FIGS
%     print(gcf,'mosquito step response.png','-dpng','-r300')
% end
% 
% 
% 
% 
% 
% 
% %% Tsevi plot
% exps = {'pert_5ms','pert_10ms','pert_20ms'};
% prop = {'forward_vel','pitch'};
% colors = [10,40,200];
% color_mean = [1,1,1];
% insect = fly ;% insect to plot
% 
% for subplot_idx = 1:1:2
%     for k = 1:1:3
%         exp_name = exps{k} ;% name of pertubation
%         hold on
% 
%         propip = fly.(exp_name).get_prop(prop{subplot_idx});
%         time = fly.(exp_name).time_vec;
%         ax1 = subplot(2,1,subplot_idx);
%         plotter_obj.all_data_plot(time,propip,'all_data_alpha',0.1,'all_data_color_idx',colors(k));hold on
% 
%     end
%     % subplot(2,1,1)
%     % ylim([40,70]) ; set(gca,'fontsize',14);
% 
%     for k = 1:1:3
%         exp_name = exps{k}; % name of pertubation
%         hold on
% 
%         propip = fly.(exp_name).get_prop(prop{subplot_idx});
%         time = fly.(exp_name).time_vec;
%         ax1 = subplot(2,1,subplot_idx);
% 
%         plotter_obj.mean_plot(insect.(exp_name),propip,prop{subplot_idx},'plot_all_data',0,'mean_color_idx',colors(k));
%         plotter_obj.mov_plot(insect.(exp_name),prop{subplot_idx},13,prop{subplot_idx},'color_idx',color_struct_fly.color_idx_mov)
%     end
%     legend({'5ms','10ms','20ms'});
% 
% end
% 
% 
% 
% %% mean and "model" plots
% 
% %% Tsevi plot (THIS)
% exp_name = 'pert_step' % name of pertubation
% pert = 0 % % pertubation (step - 0, 100 ms - 100, 60 ms - 60...)
% 
% insect = fly % insect to plot
% plotter_obj.plot_pitch_for_vel_mean(insect,exp_name,pert,max_time_xax,color_struct_fly)
% 
% W = 425 ;
% H = 535 ;
% 
% set(gcf,'position',[200 50 W H])
% subplot(2,1,1)
% ylim([25,70]) ; set(gca,'fontsize',14);
% xlim([-25 250])
% set(gca,'YTick',[30,50,70])
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% xlabel('') ;
% sgtitle('Fruit fly','Fontsize',16)
% 
% subplot(2,1,2)
% ylim([-0.3,0.3]) ;
% set(gca,'ytick',[-0.3 0 0.3]) ;
% xlim([-25 250])
% set(gca,'fontsize',14);
% set(gca,'TickDir','out');
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% if SAVE_FIGS
%     print(gcf,'fly step response.png','-dpng','-r300')
% end
% 
% 
% insect = mos % insect to plot
% plotter_obj.plot_pitch_for_vel_mean(insect,exp_name,pert,max_time_xax,color_struct_mos)
% 
% set(gcf,'position',[200 50 W H])
% subplot(2,1,1)
% ylim([25,70]) ; set(gca,'fontsize',14);
% xlim([-25 250])
% set(gca,'YTick',[30,50,70])
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% xlabel('') ;
% sgtitle('Mosquito','Fontsize',16)
% 
% subplot(2,1,2)
% ylim([-0.3,0.3]) ;
% set(gca,'ytick',[-0.3 0 0.3]) ;
% xlim([-25 250])
% set(gca,'fontsize',14);
% set(gca,'TickDir','out','TickLength',[0.01, 0.025]);
% 
% if SAVE_FIGS
%     print(gcf,'mosquito step response.png','-dpng','-r300')
% end
% 
% 
% %% TSEVI: add vz plot in the same format as above
% 
% % fly
% exp_name = 'pert_step' % name of pertubation
% pert = 0 % % pertubation (step - 0, 100 ms - 100, 60 ms - 60...)
% insect = fly % insect to plot
% mov = 13 % "model" movie
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_fly,max_time_xax)
% xlim([-20,250])
% xx = subplot(3,1,1) ;
% delete(xx) ;
% xx = subplot(3,1,2) ;
% delete(xx) ;
% 
% % mosquito plot
% insect = mos % insect to plot
% mov = 580 % "model" movie
% 
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_mos,max_time_xax)
% xlim([-20,250])
% xx = subplot(3,1,1) ;
% delete(xx) ;
% xx = subplot(3,1,2) ;
% delete(xx) ;
% 
% 
% 
% %%
% % --------------Step pertubation------------------------------
% exp_name = 'pert_step' % name of pertubation
% pert = 0 % % pertubation (step - 0, 100 ms - 100, 60 ms - 60...)
% insect = fly % insect to plot
% mov = 13 % "model" movie
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_fly,max_time_xax)
% xlim([-20,250])
% 
% % mosquito plot
% insect = mos % insect to plot
% mov = 580 % "model" movie
% 
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_mos,max_time_xax)
% xlim([-20,250])
% 
% % --------------60 ms ------------------------------
% 
% exp_name = 'pert_60ms'
% pert = 60 % pertubation (step - 0, 100 ms - 100, 60 ms - 60...)
% insect = fly
% mov = 70
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_mos,max_time_xax)
% xlim([-20,250])
% 
% insect = mos
% mov = 596
% plotter_obj.plot_pitch_for_vel_z_vel_model_mov_mean(insect,exp_name,mov,pert,color_struct_mos,max_time_xax)
% xlim([-20,250])
% 
% 
% %% Violin delta angle
% % ---- colors (defined by the indices in color mat) --------------
% fly_col_idx = 190 % color of fly
% mos_col_idx = 40 % color of mosquito
% % time to plot violin (x axis, the initial point of the violin)
% time_to_violin = linspace(-15,230,10)
% 
% %----- violin plot properties-----------
% prop_name = 'vel_xy_ang_flat' % property, a list of all properies: fly.pert_60ms.insect_prop
% fly_color = plotter_obj.col_mat(fly_col_idx,:) % color of the fly's violin
% mos_color = plotter_obj.col_mat(mos_col_idx,:) % color of the mosquito's violin
% 
% 
% % 60 ms pertubation -----------------------------
% exp_name_cell = {'pert_60ms'}  % pertubation
% exp_name = exp_name_cell{1}
% pert = 60 % used to plot the pertubation as a gray box
% 
% % plot violin-------------
% figure
% fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
% mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);
% f_norm_vec = max([fly_fvec;mos_fvec]); % normalize
% 
% %f_norm_vec(:) = max(max([fly_fvec;mos_fvec])) / 4 ;
% 
% ax1 = subplot(2,1,2)
% plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
% plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0,'box_xdata',0)
% plotter_obj.pert_plot(pert,0,1,ax1)
% ylabel('delta velocity angle [deg]')
% ylim([0,180])
% 
% % step pertubation ------------
% exp_name_cell = {'pert_step'}
% exp_name = exp_name_cell{1}
% pert = 0 % used to plot the pertubation as a gray box
% 
% 
% fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
% mos_fvec = get_fvec(exp_name_cell,mos,prop_name,time_to_violin);
% f_norm_vec = max([fly_fvec;mos_fvec]);
% 
% ax2 =subplot(2,1,1)
% plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0)
% plotter_obj.violin_plot(mos.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,mos_color,{'fly','mosquito'},'scatter_loc',0)
% plotter_obj.pert_plot(pert,0,1,ax2)
% ylabel('delta velocity angle [deg]')
% ylim([0,180])
% set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on
% 
% 
% 
% %% violin for 40ms pertubation
% % fly_col_idx = 190
% % mos_col_idx = 40
% %
% % time_to_violin = linspace(-15,230,10)
% % exp_name = 'pert_40ms'
% % prop_name = 'vel_xy_ang_flat'
% % fly_color = plotter_obj.col_mat(fly_col_idx,:)
% % mos_color = plotter_obj.col_mat(mos_col_idx,:)
% 
% 
% % exp_name_cell = {'pert_40ms'}
% % pert = 40
% % figure
% % fly_fvec = get_fvec(exp_name_cell,fly,prop_name,time_to_violin);
% %
% % f_norm_vec = max([fly_fvec;mos_fvec]);
% %
% % ax1 = subplot(1,1,1)
% % plotter_obj.violin_plot(fly.(exp_name),prop_name,time_to_violin,f_norm_vec,prop_name,fly_color,'fly','scatter_loc',0,'box_xdata',0)
% % plotter_obj.pert_plot(pert,0,1,ax1)
% % ylabel('delta velocity angle [deg]')
% % title('40ms pertubation')
% % ylim([0,180])
% % xlim([-50,300])
% % set([ax1,ax2], 'LineWidth', 3,'TickLength',[0.00,0.00]);box on
% 

