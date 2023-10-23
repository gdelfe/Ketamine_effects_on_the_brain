function load_spectrograms(BRAIN_reg_rec_dir)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');

load(strcat(spec_dir,'\spec_rec_B.mat'))
load(strcat(spec_dir,'\spec_rec_L.mat'))
load(strcat(spec_dir,'\spec_rec_M.mat'))
load(strcat(spec_dir,'\spec_rec_H.mat'))
load(strcat(spec_dir,'\spec_tf.mat'))



end 