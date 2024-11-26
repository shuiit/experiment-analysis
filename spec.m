%% get pix count
all_body_quats = 109
all_pix_counts_cam1=cell(length(all_body_quats),1);
all_pix_counts_cam2=cell(length(all_body_quats),1);
sparse_folder_path = 'J:\My Drive\dark 2022\2023_08_09_60ms\hull\hull_Reorder\'

mov_num = '109'
save_file = fullfile(sparse_folder_path,'all_pix_counts_cam2.mat')

for cam = 1:1:4

        cam_name = sprintf('_cam%d',cam)
        save_file = fullfile(sparse_folder_path,['all_pix_counts',cam_name,'.mat'])
    for mov_ind=1:all_body_quats
        dir_mov_name = sprintf('mov%d',mov_ind);
        file_name = fullfile(sparse_folder_path,[dir_mov_name,'\',dir_mov_name,cam_name,'_sparse.mat']);

        if ~exist(file_name, 'file')
            continue
        end

        disp(dir_mov_name)
        qq=load(file_name);
        pix_count=arrayfun(@(x)length(x.indIm) ,qq.frames ,'UniformOutput' ,true);

        all_pix_counts_cam2{mov_ind}=pix_count;

    end

    save(save_file,'all_pix_counts_cam2')
end
%%

sparse_folder_path = 'J:\My Drive\dark 2022\2023_08_09_60ms\hull\hull_Reorder\'
% sparse_folder_path = 'J:\My Drive\dark 2022\2024_11_12_darkan\hull\hull_Reorder\'

figure
for cam = 1:1:4
cam_name = sprintf('_cam%d',cam);
save_file = fullfile(sparse_folder_path,['all_pix_counts',cam_name,'.mat']);
load(save_file);
pix_count=all_pix_counts_cam2{mov_ind};
 [ampSpec, f] = myFFT(pix_count - mean(pix_count), Fs, false)
 subplot(2,2,cam)

        plot(f, ampSpec, '.-') ;
        title(cam_name)
        xlabel('Frequency, $f$ [Hz]') ;
        grid on ; box on ;
        axis tight
        xlim([0,500])
end



%% frquency from images
%  fourier of amount of pixels

sparse_folder_path = 'J:\My Drive\dark 2022\2023_08_09_60ms\hull\hull_Reorder\'
figure
for cam = 1:1:4
cam_name = sprintf('_cam%d',cam);
save_file = fullfile(sparse_folder_path,['all_pix_counts',cam_name,'.mat']);
load(save_file);

Fs=16000%20e3;
flap_freqs=cell(length(all_body_quats),1);
tvecs=cell(length(all_body_quats),1);
all_freqs=nan(4,length(all_body_quats));

for mov_ind=1:length(all_pix_counts_cam2)
    % for mov_ind=1:489
    % for mov_ind=8 % only 4cams
    disp(mov_ind);
    %     sparse_folder_path=folder_names(mov_ind);
    %     mov_num=num2str(mov_names(mov_ind));
    
    %     preds_path=fullfile(sparse_folder_path,['mov',mov_num,'_body.h5']);
    %     frame_inds = double(unique(squeeze(h5read(preds_path,'/frameInds'))));
    %     x_ms=round((meta_data.startFrame+double(frame_inds))*dt*1000,2);
    
    %     pix_count=all_pix_counts_cam1{mov_ind};
    pix_count=all_pix_counts_cam2{mov_ind};
    
    % Fs=20e3;
    %     figure;
    %     plot(pix_count)
    window_size=31*30;
    
    if window_size>length(pix_count)
        continue
    end
    
    [S, tvec, f] = mySpectrogram(pix_count-mean(pix_count), Fs, hann(window_size*1.2),false, true, false);
    
    flap_freq=nan(length(tvec),1);
    for time_ind=1:length(tvec)
        %         [pks,locs] = findpeaks(S(f>200,time_ind),f(f>200),'MinPeakProminence',2);
        [pks,locs] = findpeaks(S((f>1100)&(f<1300),time_ind),f((f>1100)&(f<1300)),'MinPeakProminence',2);
        [~,m_ind]=max(pks);

        if ~isempty(m_ind)
            flap_freq(time_ind,1)=locs(m_ind);
            %             flap_freq(time_ind,2)=locs(m_ind);
            
            %             if flap_freq(time_ind,1)<500
            %                 keyboard
            %             end
        end
        
        %         [pks,locs] = findpeaks(S((f>500*2)&(f<850*2),time_ind),f((f>500*2)&(f<850*2)),'MinPeakProminence',2);
        %         [~,m_ind]=max(pks);
        %         if ~isempty(m_ind)
        %             flap_freq(time_ind,1)=locs(m_ind)/2;
        %         end
        
        %         if length(locs)>1
        %             flap_freq(time_ind,1)=locs(1);
        %             flap_freq(time_ind,2)=locs(2)/2;
        %         end
    end
    tvecs{mov_ind}=tvec;%+x_ms(1)*1e-3;
    flap_freqs{mov_ind}=flap_freq;
    
    for time_ind=1:length(flap_freqs{mov_ind})
        all_freqs(time_ind,mov_ind)=flap_freqs{mov_ind}(time_ind);
    end
end


hold on;
cam1 = cell2mat(flap_freqs);
plot(tvec,cam1(:,2) - cam1(1,2))
end

% all_freqs(all_freqs==0)=nan;