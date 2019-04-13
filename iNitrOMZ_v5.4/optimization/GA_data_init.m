function Data = GA_data_init(bgc,data_raw,tracers)
%GA_DATA_INIT processes data with a spline to cover the model grid
min_std_fraction = 0.1;
Data.name = tracers;
Data.val=nan(length(tracers),bgc.nz);
Data.std=nan(length(tracers),bgc.nz);

fieldname = fields(data_raw);
for i = 1 : length(tracers)
    if sum(strcmp(tracers{i},fieldname)) == 1
        Data.val(i,:) =  data_raw.(tracers{i});
        Data.std(i,:) =  data_raw.([tracers{i},'std']);
    elseif sum(strcmp(tracers{i},fieldname)) == 0 
        Data.val(i,:) = nan(1,bgc.nz);
        Data.std(i,:) = nan(1,bgc.nz);
    else
        error('strcmp(tracers,fieldname) should be 1 or 0');
    end
    maxstd = max(Data.std(i,:));minstd = min_std_fraction*maxstd;
    tmp = Data.std(i,:);tmp(tmp<minstd)=minstd;Data.std(i,:)=tmp;
end
