function [pred dat] = minmax(prein_in, dat_in)

	mmax = repmat(nanmax(nanmax(prein_in,[],2), nanmax(dat_in.val,[],2)),1,size(prein_in,2));
	mmin = repmat(nanmin(nanmin(prein_in,[],2), nanmin(dat_in.val,[],2)),1,size(prein_in,2));

	pred = (prein_in - mmin)./(mmax-mmin);
	dat.val = (dat_in.val - mmin)./(mmax-mmin);
