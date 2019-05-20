%% only works for mi_size = 13

function [ YUV_PSNR ] = calcMetricsRGB44410( REF, REC, representation_type, H, W, mi_size, output_folder )

% process REF
for j = 1:mi_size
    for i = 1:mi_size
        ref_4DLF_VIEWS(j,i,:,:,:) = bitshift(imread(strcat(REF,sprintf('%03d_%03d.ppm',i,j))),-6);
    end
end

% process REC
if representation_type == 4 % 4DLF  SAIs
    for j = 1:mi_size
        for i = 1:mi_size
            rec_4DLF_VIEWS(j,i,:,:,:) = bitshift(imread(strcat(REC,sprintf('%03d_%03d.ppm',i,j))),-6);
        end
    end
end

for j = 1:mi_size
    for i = 1:mi_size
        [Y_PSNR(j,i) YUV_PSNR(j,i) Y_SSIM(j,i)]=QM(squeeze(ref_4DLF_VIEWS(j,i,:,:,:)),squeeze(rec_4DLF_VIEWS(j,i,:,:,:)),10,10);
    end
end

fileID = fopen( strcat(output_folder,'avg_psnr_y.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'avg_psnr_yuv.txt'), 'a' );
fprintf(fileID, "%f\n",mean(YUV_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'avg_yssim.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_SSIM(:)));
fclose(fileID);
end

