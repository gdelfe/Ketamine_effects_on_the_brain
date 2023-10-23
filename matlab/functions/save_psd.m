
% save PSD into dedicated folder inside recording folder 

function save_psd(psd, BRAIN_reg_rec_dir)

psd_dir = strcat(BRAIN_reg_rec_dir,'\psd');
if ~exist(psd_dir, 'dir')
    mkdir(psd_dir)
end
    
save(strcat(psd_dir,'\psd.mat'),'psd')

end 