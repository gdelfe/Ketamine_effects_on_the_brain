
function psd = load_psd(dir_rec)

psd_dir = strcat(dir_rec,'\psd');
    
load(strcat(psd_dir,'\psd.mat'));

end 