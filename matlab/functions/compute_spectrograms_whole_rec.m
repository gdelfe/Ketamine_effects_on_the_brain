
function [spec_rec] = compute_spectrograms_whole_rec(lfp_all, fs, start, ends, spec_par, ch_range)

pad = 2;
tot_min = size(lfp_all.B,1); % total length in min
spec_rec_B = []; spec_rec_L = []; spec_rec_M = []; spec_rec_H = [];

if nargin < 6 % if the analysis regards the whole HPC and not just a HPC subarea
    ch_range = 1:size(lfp_all.B(1,start:ends,:),3);
end

for min = 1:tot_min

    X = sq(lfp_all.B(min,start:ends,ch_range))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_B = cat(3,spec_rec_B,spec);
    
    X = sq(lfp_all.L(min,start:ends,ch_range))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_L = cat(3,spec_rec_L,spec);
    
    X = sq(lfp_all.M(min,start:ends,ch_range))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_M = cat(3,spec_rec_M,spec);
    
    X = sq(lfp_all.H(min,start:ends,ch_range))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_H = cat(3,spec_rec_H,spec);
    
end

    spec_rec.B = spec_rec_B;
    spec_rec.L = spec_rec_L;
    spec_rec.M = spec_rec_M;
    spec_rec.H = spec_rec_H;
    
    spec_rec.f = f;
    spec_rec.t = ti; 
    spec_rec.ts = linspace(0,size(X,2)/fs,size(X,2));


end 