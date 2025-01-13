

clear
close all
clc

experiment_path = 'H:\My Drive\dark 2022\2022_03_03\'
segmentation_path = [experiment_path,'hull\hull_Reorder\mov19\Segmentation\mov19_seg']
easywand_path = [experiment_path,'3+4_post_03_03_2022_skip5_easyWandData']


experiment_path = 'H:\My Drive\dark 2022\2024_11_12_darkan\'
segmentation_path = [experiment_path,'hull\hull_Reorder\mov11\Segmentation\mov11_seg']
easywand_path = [experiment_path,'coefs_12_11_24_easyWandData']

load(easywand_path)
load(segmentation_path)
%%
 


for j= 1:1:4
[R,K,X0,H] = decompose_dlt(easyWandData.coefs(:,j),easyWandData.rotationMatrices(:,:,j)');
camera(:,:,j) = [K,R,X0];
rotation_allcam(:,:,j) =  R; 
translation(:,:,j) = X0; 
pmdlt{j} = [K*R,-K*R*X0];
intrinsic{j} = K;
camera_center{j} = X0;

end


r = vrrotvec(rotation_allcam(3,:,1),[0,0,1])
m = vrrotvec2mat(r)

%%

% ray_ndc = [(pixels[1] - self.cx)/self.fx,(pixels[0] - self.cy)/self.fy,1]
% np.dot(self.R.T,ray_ndc) + self.X0.T
% 


for index = 1:1:4
fx_cam(index) = intrinsic{index}(1,1);
fy_cam(index) = intrinsic{index}(2,2);
cx_cam(index) = intrinsic{index}(1,3);
cy_cam(index) = intrinsic{index}(2,3);
rotation = rotation_allcam(:,:,index);
end
%%
cm_3d= []
time_vec = []
cm_idx = 1;
start_frame = seg.st_enfr(1);
for frame = 1:1:length(seg.body{index})
    
        cam_index = 1;
        ray_center_pixel = [];
        ray = [];
for index = 1:1:4

% binaryImage = ImfromSp([800, 1280], seg.all{index}(frame).indIm); % Ensure this function is optimized
% CC = bwconncomp(binaryImage); % Faster and memory-efficient labeling
% if length(CC.PixelIdxList) > 1
%     continue
% end


% labeledImage = bwlabel(binaryImage);
% stats = regionprops(labeledImage, 'Area');
% if length([stats.Area])> 1
% [~, largestBlobIndex] = max([stats.Area]);
% largestBlob = ismember(labeledImage, largestBlobIndex);
% [row,col] = find(largestBlob);
% cm(index,:) = mean([row,col]);
% else
if size(seg.all{index}(frame).indIm,1 ) > 700
    mean_cm = mean(seg.all{cam_index}(frame).indIm);
    cm(cam_index,:) = mean_cm(1:2);


    % end
    ray(cam_index,:) = [(cm(cam_index,2) - cx_cam(cam_index) )/fx_cam(cam_index),(cm(cam_index,1) - cy_cam(cam_index) )/fy_cam(cam_index),1];
    
    ray_center_pixel(cam_index,:) = rotation_allcam(:,:,cam_index)'*ray(cam_index,:)' + camera_center{cam_index};
    cam_index= cam_index + 1;


end
end
if size(ray_center_pixel,1) == 4
    


center = cell2mat(camera_center)';
ray = ray_center_pixel;
cm_3d(cm_idx,:) = lineIntersect3D(ray,center);
time_vec(cm_idx) = -(800)/16 + frame*1/16;

cm_idx = cm_idx + 1;
end
    
end




%%
% plot_camera(rotation_allcam,translation,m,'standard wand')

cm_labax = (m*cm_3d')'
figure;plot(time_vec(1:end),cm_labax)

figure;scatter3(cm_labax(:,1),cm_labax(:,2),cm_labax(:,3),5,'filled');hold on

scatter3(cm_labax(1,1),cm_labax(1,2),cm_labax(1,3),100,'filled','MarkerFaceColor','red')



