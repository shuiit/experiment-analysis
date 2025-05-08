clear
close all
clc

sparse_to_load = 'J:\My Drive\dark 2022\2023_08_09_60ms\hull\hull_Reorder\mov30\mov30_cam1_sparse.mat'
sp = load(sparse_to_load)
%%
frame = 520

sparse_coordinates = sp.frames(frame).indIm;

[im] = ImfromSp([800,1280],sparse_coordinates);

bb = [max(sparse_coordinates(:,2)) + 10,max(sparse_coordinates(:,1)) + 10;min(sparse_coordinates(:,2))-10,min(sparse_coordinates(:,1))-10]
croped_bb = im(bb(2,2):bb(1,2),bb(2,1):bb(1,1))
croped_bb(croped_bb > 0) = 1
figure;imshow(1-croped_bb);hold on
icon_final = 1-croped_bb;

imwrite(icon_final, 'fly_icon.png', 'Transparency', 1);