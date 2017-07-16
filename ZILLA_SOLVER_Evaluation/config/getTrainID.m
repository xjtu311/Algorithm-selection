function ids = getTrainID(inst_seed_numbers, func)
if func.cheap
    ids = cell(length(inst_seed_numbers), 1);
    for i=1:length(ids)
        ids{i} = -i;
    end
else
    ids = get_inst_id(func.instance_filenames(inst_seed_numbers));
end