%% only works for mi_size = 13

function [ ] = calcMetricsYUV4410( REF, REC, representation_type, H, W, mi_size, output_folder )

% process REF
for j = 1:mi_size
    for i = 1:mi_size
        ref_4DLF_VIEWS(j,i,:,:,:) = double(imread(strcat(REF,sprintf('%03d_%03d.ppm',i,j)))) ;
    end
end

% process REC
if representation_type == 3 % 4DLF-PVS
    cc_spiral = spiral(mi_size);
    f = fopen(REC,'r');
    for j = 1:mi_size
        for i = 1:mi_size
            [ypos, xpos] = find(cc_spiral == (j-1)*mi_size + i);
            Y = fread(f, [W H], 'uint16');
            U = fread(f, [W H], 'uint16');
            V = fread(f, [W H], 'uint16');
            rec_4DLF_VIEWS(ypos,xpos,:,:,1) = uint16(Y');
            rec_4DLF_VIEWS(ypos,xpos,:,:,2) = uint16(U');
            rec_4DLF_VIEWS(ypos,xpos,:,:,3) = uint16(V');
        end
    end
    fclose(f);
end
if representation_type == 4 % 4DLF ind SAIs
    for j = 1:mi_size
        for i = 1:mi_size
            rec_4DLF_VIEWS(j,i,:,:,:) = imread(strcat(REC,sprintf('%03d_%03d.ppm',i,j))) ;
        end
    end
end

for j = 1:mi_size
    for i = 1:mi_size
        [Y_PSNR(j,i) YUV_PSNR(j,i) Y_SSIM(j,i)]=QM_YUV44410(squeeze(ref_4DLF_VIEWS(j,i,:,:,:)),squeeze(rec_4DLF_VIEWS(j,i,:,:,:)),16,10);
    end
end

fileID = fopen( strcat(output_folder,'avg_psnr_y.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'avg_psnr_yuv.txt'), 'a' );
fprintf(fileID, "%f\n",mean(YUV_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'avg_psnr_yssim.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_SSIM(:)));
fclose(fileID);
end

