function save_spectrograms(spec_rec, BRAIN_reg_rec_dir,name_file)

spec_dir = strcat(BRAIN_reg_rec_dir,'\spectrograms');
if ~exist(spec_dir, 'dir')
    mkdir(spec_dir)
end
    
save(strcat(spec_dir, name_file),'spec_rec')

end 