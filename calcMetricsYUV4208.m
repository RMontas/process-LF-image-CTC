%% only works for mi_size = 13
function [ ] = calcMetricsYUV4208( REF, REC, representation_type, H, W, mi_size, output_folder )

% process REF
for j = 1:mi_size
    for i = 1:mi_size
        ref_4DLF_VIEWS(j,i,:,:,:) = double(imread(strcat(REF,sprintf('%03d_%03d.ppm',i,j)))) ;
    end
end

% process REC
if representation_type == 1 % 4DLF-MI
    f = fopen(REC,'r');
    YUV{1} = (fread(f, [W H], 'uint8'))';
    YUV{2} = (fread(f, [W/2 H/2], 'uint8'))';
    YUV{3} = (fread(f, [W/2 H/2], 'uint8'))';
    YUV444=upsample(YUV);
    img_4DLF_MI(:,:,1)=uint16(YUV444(:,1:W-1,1))*4; % padding removal and bitshift
    img_4DLF_MI(:,:,2)=uint16(YUV444(:,1:W-1,2))*4;
    img_4DLF_MI(:,:,3)=uint16(YUV444(:,1:W-1,3))*4;
    rec_4DLF_VIEWS = deconstruct_lenslet_img10( img_4DLF_MI, mi_size );
end

if representation_type == 3 % 4DLF-PVS
    cc_spiral = spiral(mi_size);
    f = fopen(REC,'r');
    for j = 1:mi_size
        for i = 1:mi_size
            [ypos, xpos] = find(cc_spiral == (j-1)*mi_size + i);
            YUV{1} = (fread(f, [W H], 'uint8'))';
            YUV{2} = (fread(f, [W/2 H/2], 'uint8'))';
            YUV{3} = (fread(f, [W/2 H/2], 'uint8'))';
            YUV444=upsample(YUV);
            Y=uint16(YUV444(:,1:W-1,1))*4; % padding removal and bitshift 
            U=uint16(YUV444(:,1:W-1,2))*4; 
            V=uint16(YUV444(:,1:W-1,3))*4; 
            rec_4DLF_VIEWS(ypos,xpos,:,:,1) = Y;
            rec_4DLF_VIEWS(ypos,xpos,:,:,2) = U;
            rec_4DLF_VIEWS(ypos,xpos,:,:,3) = V;
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

fileID = fopen( strcat(output_folder,'_avg_psnr_y.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'_avg_psnr_yuv.txt'), 'a' );
fprintf(fileID, "%f\n",mean(YUV_PSNR(:)));
fclose(fileID);
fileID = fopen( strcat(output_folder,'_avg_psnr_yssim.txt'), 'a' );
fprintf(fileID, "%f\n",mean(Y_SSIM(:)));
fclose(fileID);
end

function [out] = upsample(in)

out = uint8(zeros(size(in{1}, 1), size(in{1}, 2), 3));
out(:,:,1) = in{1};
tmp = cell(3,1);

for i=2:3
    even = imfilter(double(in{i}), [-2 ; 16 ; 54 ; -4 ; 0], 'replicate', 'same');
    odd = imfilter(double(in{i}), [0 ; -4 ; 54 ; 16 ; -2], 'replicate', 'same');
    fields = zeros(2*size(even,1), size(even,2));
    fields(1:2:end,:) = even;
    fields(2:2:end,:) = odd;
    tmp{i} = fields;
end

for i=2:3
    even = (tmp{i} + 32) / 64;
    odd = imfilter(tmp{i}, [0 -4 36 36 -4], 'replicate', 'same');
    odd = (odd + 2048) / 4096; % shift = 12 i.e. (2^12); offset = 2048
    fields = zeros(size(even,1), 2*size(even,2));
    fields(:,1:2:end) = even;
    fields(:,2:2:end) = odd;
    out(:,:,i) = uint8(fields);
end
end



