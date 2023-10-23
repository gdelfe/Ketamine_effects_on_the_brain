function save_spectrograms(spec_rec, BRAIN_reg_rec_dir)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');
if ~exist(spec_dir, 'dir')
    mkdir(spec_dir)
end
    
save(strcat(spec_dir,'\spec_rec_B.mat'),'spec_rec_B')
save(strcat(spec_dir,'\spec_rec_L.mat'),'spec_rec_L')
save(strcat(spec_dir,'\spec_rec_M.mat'),'spec_rec_M')
save(strcat(spec_dir,'\spec_rec_H.mat'),'spec_rec_H')
save(strcat(spec_dir,'\spec_tf.mat'),'spec_tf')



end 