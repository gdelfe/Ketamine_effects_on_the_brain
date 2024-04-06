
% save PSD into dedicated folder inside recording folder 

function save_psd(psd, dir_rec, name_file)

psd_dir = strcat(dir_rec,'\psd');
if ~exist(psd_dir, 'dir')
    mkdir(psd_dir)
end
    
save(strcat(psd_dir,sprintf('\\%s.mat',name_file)),'psd')

end 