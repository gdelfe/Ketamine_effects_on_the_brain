
function psd = load_psd(dir_rec,name_file)

psd_dir = strcat(dir_rec,'\psd');
    
load(strcat(psd_dir,sprintf('\\%s.mat',name_file)));

end 