function [spec_rec] = load_spectrograms(BRAIN_reg_rec_dir)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');

load(strcat(spec_dir,'\spec_rec_all_HPC.mat'))

end 