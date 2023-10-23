
% save PSD into dedicated folder inside recording folder 

function save_psd(psd_B, psd_L, psd_M, psd_H, f, BRAIN_reg_rec_dir)

psd_dir = strcat(BRAIN_reg_rec_dir,'\psd');
if ~exist(psd_dir, 'dir')
    mkdir(psd_dir)
end
    
save(strcat(psd_dir,'\psd_B.mat'),'psd_B')
save(strcat(psd_dir,'\psd_L.mat'),'psd_L')
save(strcat(psd_dir,'\psd_M.mat'),'psd_M')
save(strcat(psd_dir,'\psd_H.mat'),'psd_H')
save(strcat(psd_dir,'\psd_f.mat'),'f')

end 