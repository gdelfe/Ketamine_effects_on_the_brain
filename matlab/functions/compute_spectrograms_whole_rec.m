
function [spec_rec] = compute_spectrograms_whole_rec(lfp_all, start, ends, spec_par)

pad = 2;
tot_min = size(lfp_all.B,1); % total length in min
spec_rec_B = []; spec_rec_L = []; spec_rec_M = []; spec_rec_H = [];

for min = 1:tot_min

    X = sq(lfp_all.B(min,start:ends,:))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_B = cat(3,spec_rec_B,spec);
    
    X = sq(lfp_all.L(min,start:ends,:))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_L = cat(3,spec_rec_L,spec);
    
    X = sq(lfp_all.M(min,start:ends,:))';
    [spec, ~, ~] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_M = cat(3,spec_rec_M,spec);
    
    X = sq(lfp_all.H(min,start:ends,:))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_H = cat(3,spec_rec_H,spec);
    
    spec_rec.B = spec_rec_B;
    spec_rec.L = spec_rec_L;
    spec_rec.M = spec_rec_M;
    spec_rec.H = spec_rec_H;
    
    spec_rec.f = f;
    spec_rec.t = ti;
    
end