
path = 'I:\.shortcut-targets-by-id\1OA70vOJHDfV63DqG7LJCifTwW1h055ny\2024 Flight in the dark paper\data_exchange\fly\all_data\'
file_name = '60ms_all_data.csv'
exp = readtable([path,file_name]);
%%
frame_rate = 16000;
header = exp.Properties.VariableNames;
insect_prop = dictionary([string(header)],[1:1:length(header)]);

delta_t = time_vec(2)-time_vec(1);    
time_vec = [min(exp.time):delta_t:max(exp.time)];
movie_numbers = unique(exp.mov_num);
empty_movie = nan(length(time_vec),size(movie,2),length(movie_numbers));
for mov_num = [1:1:length(movie_numbers)]
    movie = exp(exp.mov_num == mov_num,:);
    empty_movie(ismember(time_vec,movie.time),:,mov_num) = table2array(movie);
end
prop2plot = 'vel_xy_ang'
prop = squeeze(mean(empty_movie(:,insect_prop(prop2plot),:),3,'omitnan'))
time = squeeze(mean(empty_movie(:,insect_prop('time'),:)))

hold on
for j = 1:1:length(movie_numbers)
plot(time_vec,empty_movie(:,insect_prop(prop2plot),j),'b');
end
plot(time_vec,prop,'r')

