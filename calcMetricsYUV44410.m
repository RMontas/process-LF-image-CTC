%% only works for mi_size = 13

function [ ] = calcMetricsYUV4410( REF, REC, representation_type, H, W, mi_size, output_folder )

% process REF
mi_size_ref = 13;
for j = 1:mi_size_ref
    for i = 1:mi_size_ref
        ref_4DLF_VIEWS(j,i,:,:,:) = double(imread(strcat(REF,sprintf('%03d_%03d.ppm',i,j)))) ;
    end
end

% process REC
if representation_type == 0 % LL
    f = fopen(REC,'r');
    Y = fread(f, [W H], 'uint16');
    U = fread(f, [W H], 'uint16');
    V = fread(f, [W H], 'uint16');
    img_LL(:,:,1)=Y';
    img_LL(:,:,2)=U';
    img_LL(:,:,3)=V';
    [img_LL_rgb] = ycbcr2rgbn(double(img_LL)./(2^10-1).*(2^16-1),16);  % to RGB444 10 bit   
    run('LFToolbox0.4/LFMatlabPathSetup')
    LFMatlabPathSetup;
    LFUtilProcessWhiteImages;
    imageName = extractAfter(REF,"RGB44410/"); imageName=extractBefore(imageName,"/");
    [~, LFCol] = YUVLenslet2LF( double(img_LL_rgb), sprintf('RAW/%s__Decoded.mat',imageName)); %16bit
     %max(max(max(max(img_LL_rgb))))
      %max(max(max(max(LFCol))))
    rec_4DLF_VIEWS = uint16(LFCol(2:14,2:14,:,:,:)); %16bit
    for j = 1:mi_size
        for i = 1:mi_size
             rec_4DLF_VIEWS(j,i,:,:,:) = rgb2ycbcrn(double(squeeze(rec_4DLF_VIEWS(j,i,:,:,:)))/(2^16-1),10); 
        end
    end
          %max(max(max(max(rec_4DLF_VIEWS))))
end

if representation_type == 1 % 4DLF-MI
    f = fopen(REC,'r');
    Y = fread(f, [W H], 'uint16');
    U = fread(f, [W H], 'uint16');
    V = fread(f, [W H], 'uint16');
    img_4DLF_MI(:,:,1)=Y';
    img_4DLF_MI(:,:,2)=U';
    img_4DLF_MI(:,:,3)=V';
    rec_4DLF_VIEWS = deconstruct_lenslet_img10( img_4DLF_MI, mi_size );
end

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
if representation_type == 33 % 4DLF-PVS-SCL
    num_Layers = 6;
    layerMask =  [ 6 6 4 6 3 6 4 6 3 6 4 6 6 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   4 6 2 6 4 6 2 6 4 6 2 6 4 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   3 6 4 6 3 6 4 6 3 6 4 6 3 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   4 6 2 6 4 6 1 6 4 6 2 6 4 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   3 6 4 6 3 6 4 6 3 6 4 6 3 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   4 6 2 6 4 6 2 6 4 6 2 6 4 ;
                   6 5 6 5 6 5 6 5 6 5 6 5 6 ;
                   6 6 4 6 3 6 4 6 3 6 4 6 6 ];
    cc_spiral = spiral(mi_size);
    f = fopen(REC,'r');
    for l = 1:num_Layers
        for j = 1:mi_size
            for i = 1:mi_size
                [ypos, xpos] = find(cc_spiral == (j-1)*mi_size + i);
                if layerMask(ypos,xpos) == l
                    Y = fread(f, [W H], 'uint16');
                    U = fread(f, [W H], 'uint16');
                    V = fread(f, [W H], 'uint16');
                    rec_4DLF_VIEWS(ypos,xpos,:,:,1) = uint16(Y');
                    rec_4DLF_VIEWS(ypos,xpos,:,:,2) = uint16(U');
                    rec_4DLF_VIEWS(ypos,xpos,:,:,3) = uint16(V');
                end
            end
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

if representation_type == 5 % 4DLF-PVS SERPENTINE
    serpentine_scaning_order = zeros(mi_size, mi_size);
    %gen serpentine scan order
    xx = 0;
    yy = 0;
    lastdir = 1; % 0 L<---R 1 L--->R
    currdir = 1;
    for j = 1:mi_size*mi_size
        yy = ceil((j)/mi_size);
        currdir = mod(yy, 2);
        if lastdir == currdir
            if currdir == 1
                xx = xx + 1; % 1 odd - add
            else
                xx = xx - 1; % 0 even - subtract
            end
        end
        lastdir = currdir;
        serpentine_scaning_order(yy,xx)=j;
    end
    f = fopen(REC,'r');
    for j = 1:mi_size
        for i = 1:mi_size
            [ypos, xpos] = find(serpentine_scaning_order == (j-1)*mi_size + i);
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

%max(max(max(max(max(ref_4DLF_VIEWS)))))
%max(max(max(max(max(rec_4DLF_VIEWS)))))
%[PSNR_Y, PSNR_U, PSNR_V, PSNR_YUV, PSNR_Y_mean, PSNR_U_mean, PSNR_V_mean, PSNR_YUV_mean] = ComputePSNR_YUV444_10bpp(rec_4DLF_VIEWS,uint16(ref_4DLF_VIEWS/65535*1023));
for j = 1:mi_size
    for i = 1:mi_size
        [Y_PSNR(j,i) YUV_PSNR(j,i) Y_SSIM(j,i)]=QM_YUV44410(squeeze(ref_4DLF_VIEWS(j,i,:,:,:)),squeeze(rec_4DLF_VIEWS(j,i,:,:,:)),16,10);
        
    end
end

mean(YUV_PSNR(:))
mean(Y_PSNR(:))

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

