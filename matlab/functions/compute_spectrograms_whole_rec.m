
function [spec_rec_B, spec_rec_L, spec_rec_M, spec_rec_H , spec_tf] = compute_spectrograms_whole_rec(lfp_B_all,lfp_L_all,lfp_M_all,lfp_H_all,start,ends, spec_par)

pad = 2;
tot_min = size(lfp_B_all,1); % total length in min
spec_rec_B = []; spec_rec_L = []; spec_rec_M = []; spec_rec_H = [];

for min = 1:tot_min

    X = sq(lfp_B_all(min,start:ends,:))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_B = cat(3,spec_rec_B,spec);
    
    X = sq(lfp_L_all(min,start:ends,:))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_L = cat(3,spec_rec_L,spec);
    
    X = sq(lfp_M_all(min,start:ends,:))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_M = cat(3,spec_rec_M,spec);
    
    X = sq(lfp_H_all(min,start:ends,:))';
    [spec, f, ti] = tfspec(X, spec_par.tapers, spec_par.fs, spec_par.dn, spec_par.fk, pad, 0.05,1);
    spec_rec_H = cat(3,spec_rec_H,spec);
    
    spec_tf.f = f;
    spec_tf.t = ti;
    
end