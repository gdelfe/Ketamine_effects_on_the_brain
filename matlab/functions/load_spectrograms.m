function load_spectrograms(BRAIN_reg_rec_dir)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');

load(strcat(spec_dir,'\spec_rec.mat'))

end 