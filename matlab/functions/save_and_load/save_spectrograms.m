function save_spectrograms(spec_rec, BRAIN_reg_rec_dir)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');
if ~exist(spec_dir, 'dir')
    mkdir(spec_dir)
end
    
save(strcat(spec_dir,'\spec_rec_all_HPC.mat'),'spec_rec')

end 