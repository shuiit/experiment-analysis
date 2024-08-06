
function insect_fvec = get_fvec(exp_name_cell,insect,prop_name,time_vec)
for i = 1:1:length(exp_name_cell)
    insect_fvec(i,:) = insect.(exp_name_cell{i}).calc_ksdensity(prop_name,time_vec);
end
end