% 10 bpp / 13x13 SAIs

%render SAIs from RAW LFs
clear all
close all

%% config flags
compareGeneratedWithOriginals = 0;
%% color flags
RGB44410 = 0;
YUV44410 = 1;
YUV4208 = 0;
%% representation type flag
REPSAI = 0; %% PPM format
REPPVSSPIRAL = 0; %% YUV
REPMI = 0; %% YUV
REPLL = 0; %% Lenslet
REPPVSSERPENTINE = 0; %% YUV
REPPVSSERPENTINE_N = 1;
N_MIs = 9;

mkdir('SAI')
mkdir('PVSSPIRAL')
mkdir('MI')
mkdir('LL')

mkdir('SAI/RGB44410')
mkdir('PVSSPIRAL/RGB44410')
mkdir('MI/RGB44410')
mkdir('LL/RGB44410')

mkdir('SAI/YUV44410')
mkdir('PVSSPIRAL/YUV44410')
mkdir('MI/YUV44410')
mkdir('LL/YUV44410')

mkdir('SAI/YUV4208')
mkdir('PVSSPIRAL/YUV4208')
mkdir('MI/YUV4208')
mkdir('LL/YUV4208')

mkdir('PVSSERPENTINE')
mkdir('PVSSERPENTINE/YUV44410')
mkdir('PVSSERPENTINE/RGB44410') % not implemented
mkdir('PVSSERPENTINE/YUV4208') % not implemented


run('LFToolbox0.4/LFMatlabPathSetup')
LFMatlabPathSetup;
LFUtilProcessWhiteImages;
listing = dir('RAW/*.LFR');
nfiles = size(listing);
nfiles = nfiles(1);

for i = 1:12
    if ~exist(sprintf('RAW/lensletI%02d.ppm',i),'file')
        LFUtilDecodeLytroFolder(strcat('./', listing(i).name), []); 
        movefile('tmplenslet.ppm',sprintf('RAW/lensletI%02d.ppm',i));
    end
    [~, LFCol] = PPMLenslet2LF( sprintf('RAW/lensletI%02d.ppm',i), sprintf('RAW/I%02d__Decoded.mat',i));
    LL = imread(sprintf('RAW/lensletI%02d.ppm',i));
    %%
    if REPLL == 1
        genLL(i, double(LL), RGB44410, YUV44410, YUV4208);
    end
    %%
    if REPSAI == 1
        genSAI(i, LFCol, compareGeneratedWithOriginals, RGB44410, YUV44410, YUV4208);
    end
    %%
    if REPPVSSPIRAL == 1
        genPVSSPIRAL(i, LFCol, RGB44410, YUV44410, YUV4208);
    end
    if REPPVSSERPENTINE == 1
        genPVSSERPENTINE(i, LFCol, RGB44410, YUV44410, YUV4208);
    end
    if REPPVSSERPENTINE_N == 1
        mkdir(sprintf('PVSSERPENTINE_N%d',N_MIs));
        mkdir(sprintf('PVSSERPENTINE_N%d/YUV44410',N_MIs));
        mkdir(sprintf('PVSSERPENTINE_N%d/YUV4208',N_MIs));
        mkdir(sprintf('PVSSERPENTINE_N%d/RGB44410',N_MIs));
        genPVSSERPENTINE_N(i, LFCol, RGB44410, YUV44410, YUV4208, N_MIs);
    end
    %%
    if REPMI == 1
        genMI(i, LFCol, RGB44410, YUV44410, YUV4208);
    end
end

function genLL(i, LF, RGB44410, YUV44410, YUV4208)
destinFolderRGB44410 = sprintf('LL/RGB44410');
destinFolderYUV44410 = sprintf('LL/YUV44410');
destinFolderYUV4208 = sprintf('LL/YUV4208');
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

LF = LF./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
%clipping
LF(LF<0) = 0;
LF(LF>1023) = 1023;
%rounding to integer
LF = double(uint16(LF));
%going back to uint16
if RGB44410 == 1
    I = LF./(2^10 -1) .*(2^16-1);
    I = uint16(I);
    fileID = fopen( sprintf('%s/I%02d.rgb',destinFolderRGB44410,i), 'w' );
    fwrite(fileID, permute(I, [2 1 3]), 'uint16');
    fclose(fileID);
    clear I;
end
if YUV44410 == 1
    I = rgb2ycbcrn(double(LF) / (2^10-1),10);
    fileID = fopen( sprintf('%s/I%02d.yuv',destinFolderYUV44410,i), 'w' );
    fwrite(fileID, permute(I, [2 1 3]), 'uint16');
    fclose(fileID);
    clear I;
end
if YUV4208 == 1
    I = rgb2ycbcrn(double(LF) / (2^10-1),8);
    I = downsampleChroma(I);
    fileID = fopen( sprintf('%s/I%02d.yuv',destinFolderYUV4208,i), 'w' );
    fwrite(fileID, I{1}', 'uint8'); fwrite(fileID, I{2}', 'uint8'); fwrite(fileID, I{3}', 'uint8');
    fclose(fileID);
    clear I;

end
end

function genSAI(i, LF, compareGeneratedWithOriginals, RGB44410, YUV44410, YUV4208)
destinFolderRGB44410 = sprintf('SAI/RGB44410/I%02d',i);
destinFolderYUV44410 = sprintf('SAI/YUV44410/I%02d',i);
destinFolderYUV4208 = sprintf('SAI/YUV4208/I%02d',i);
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

for xx=1:15
    for yy=1:15
        A =squeeze(double(LF(xx,yy, :, :, 1:3)));
        %A(:, end+1, :) = 0; %one-pixel padding
        A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
        %clipping
        A(A<0) = 0;
        A(A>1023) = 1023;
        %rounding to integer
        A = double(uint16(A));
        %going back to uint16
        if RGB44410 == 1
            I = A./(2^10 -1) .*(2^16-1);
            I = uint16(I);
            imwrite(I, sprintf('%s/%03d_%03d.ppm',destinFolderRGB44410, yy-1, xx-1), 'MaxValue', 1023);
            clear I;
        end
        if YUV44410 == 1
            I = rgb2ycbcrn(double(A) / (2^10-1),10);
            fileID = fopen( sprintf('%s/%03d_%03d.yuv',destinFolderYUV44410, yy-1, xx-1), 'w' );
            fwrite(fileID, permute(I, [2 1 3]), 'uint16');
            fclose(fileID);
            clear I;
        end
        if YUV4208 == 1
            I = rgb2ycbcrn(double(A) / (2^10-1),8);
            I(:, end+1, :) = 0; %one-pixel padding
            I = downsampleChroma(I);
            fileID = fopen( sprintf('%s/%03d_%03d.yuv',destinFolderYUV4208, yy-1, xx-1), 'w' );
            fwrite(fileID, I{1}', 'uint8'); fwrite(fileID, I{2}', 'uint8'); fwrite(fileID, I{3}', 'uint8');
            fclose(fileID);
            clear I;
        end
    end
end
if compareGeneratedWithOriginals == 1 && RGB44410 == 1
    if i==1 || i==2 || i==4 || i==9
        j = 1;
        for xx=1:15
            for yy=1:15
                REC=imread(sprintf('%s/%03d_%03d.ppm', destinFolder, i, yy-1, xx-1));
                REF=imread(sprintf('JPEG_I%02d/%03d_%03d.ppm', i, yy-1, xx-1));
                psnrIMGs(j,xx,yy)=psnr(REF,REC);
                j = j+1;
            end
        end
    end
    psnrIMGs
end
end

function genPVSSPIRAL(i, LF, RGB44410, YUV44410, YUV4208)
destinFolderRGB44410 = sprintf('PVSSPIRAL/RGB44410/I%02d',i);
destinFolderYUV44410 = sprintf('PVSSPIRAL/YUV44410/I%02d',i);
destinFolderYUV4208 = sprintf('PVSSPIRAL/YUV4208/I%02d',i);
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

numMIs=13;
cc_spiral = spiral(numMIs);
res=size(LF);
H=res(3);
W=res(4);

% only 13x13 in this rep
PVS_RGB44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
PVS_YUV44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
for j=1:numMIs*numMIs
	[yy, xx] = find(cc_spiral == j);
    A =squeeze(double(LF(yy+1,xx+1, :, :, 1:3))); % pad views
    %A(:, end+1, :) = 0; %one-pixel padding
    A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
    %clipping
    A(A<0) = 0;
    A(A>1023) = 1023;
    %rounding to integer
    A = double(uint16(A));
    if RGB44410 == 1
        PVS_RGB44410(:,:,:,j) = uint16(A);
    end
    if YUV44410 == 1
        PVS_YUV44410(:,:,:,j) = rgb2ycbcrn(double(A) / (2^10-1),10);
    end
    if YUV4208 == 1
        aux = rgb2ycbcrn(double(A) / (2^10-1),8);
        aux(:, end+1, :,:) = 0; %one-pixel padding
        PVS_YUV4208(:,j) = downsampleChroma(aux);
        clear aux
    end
end

if RGB44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_RGB44410.rgb',destinFolderRGB44410,i), 'w' );
    fwrite(fileID, permute(PVS_RGB44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_YUV44410.yuv',destinFolderYUV44410,i), 'w' );
    fwrite(fileID, permute(PVS_YUV44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV4208 == 1
	fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_YUV4208.yuv',destinFolderYUV4208,i), 'w' );
    for j=1:numMIs*numMIs
        fwrite(fileID, PVS_YUV4208{1,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{2,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{3,j}', 'uint8');
    end
    fclose(fileID);
end
end

function genMI(i, LF, RGB44410, YUV44410, YUV4208)
destinFolderRGB44410 = sprintf('MI/RGB44410/I%02d',i);
destinFolderYUV44410 = sprintf('MI/YUV44410/I%02d',i);
destinFolderYUV4208 = sprintf('MI/YUV4208/I%02d',i);
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

numMIs=13;
res=size(LF);
H=res(3);
W=res(4);
LF_RGB44410 = zeros(numMIs, numMIs, H, W, 3, 'uint16');
LF_YUV44410 = zeros(numMIs, numMIs, H, W, 3, 'uint16');
LF_YUV4208 = zeros(numMIs, numMIs, H, W, 3, 'uint8');
MI_RGB44410 = zeros(H*numMIs, W*numMIs, 3, 'uint16');
MI_YUV44410 = zeros(H*numMIs, W*numMIs, 3, 'uint16');

for xx=1:numMIs
    for yy=1:numMIs
        A =squeeze(double(LF(xx+1,yy+1, :, :, 1:3)));
        A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
        %clipping
        A(A<0) = 0;
        A(A>1023) = 1023;
        %rounding to integer
        A = double(uint16(A));
        %going back to uint16
        if RGB44410 == 1
            I = A./(2^10 -1) .*(2^16-1);
            I = uint16(I);
            LF_RGB44410(xx, yy, :, :, :) = I;
            clear I;
        end
        if YUV44410 == 1
            I = rgb2ycbcrn(double(A) / (2^10-1),10);
            LF_YUV44410(xx, yy, :, :, :) = I;
            clear I;
        end
        if YUV4208 == 1
            I = rgb2ycbcrn(double(A) / (2^10-1),8);
            LF_YUV4208(xx, yy, :, :, :) = I;
            clear I;
        end
    end
end

if RGB44410 == 1
    MI_RGB44410 = reconstruct_lenslet_img( LF_RGB44410, numMIs );
    fileID = fopen( sprintf('%s/I%02d_MI_RGB44410.rgb',destinFolderRGB44410,i), 'w' );
    fwrite(fileID, permute(MI_RGB44410, [2 1 3]), 'uint16');
    fclose(fileID);
end
if YUV44410 == 1
    MI_YUV44410 = reconstruct_lenslet_img( LF_YUV44410, numMIs );
    fileID = fopen( sprintf('%s/I%02d_MI_YUV44410.yuv',destinFolderYUV44410,i), 'w' );
    fwrite(fileID, permute(MI_YUV44410, [2 1 3]), 'uint16');
    fclose(fileID);
end
if YUV4208 == 1
    aux_MI_YUV4208 = reconstruct_lenslet_img( LF_YUV4208, numMIs );
    aux_MI_YUV4208(:, end+1, :) = 0; %one-pixel padding
    MI_YUV4208 = downsampleChroma(aux_MI_YUV4208);
    fileID = fopen( sprintf('%s/I%02d_MI_YUV4208.yuv',destinFolderYUV4208,i), 'w' );
    fwrite(fileID, permute(MI_YUV4208{1}, [2 1 3]), 'uint8');
    fwrite(fileID, permute(MI_YUV4208{2}, [2 1 3]), 'uint8');
    fwrite(fileID, permute(MI_YUV4208{3}, [2 1 3]), 'uint8');
    fclose(fileID);
end

end

function genPVSSERPENTINE(i, LF, RGB44410, YUV44410, YUV4208)
destinFolderRGB44410 = sprintf('PVSSERPENTINE/RGB44410/I%02d',i);
destinFolderYUV44410 = sprintf('PVSSERPENTINE/YUV44410/I%02d',i);
destinFolderYUV4208 = sprintf('PVSSERPENTINE/YUV4208/I%02d',i);
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

numMIs=13;
res=size(LF);
H=res(3);
W=res(4);

% only 13x13 in this rep
PVS_RGB44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
PVS_YUV44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
xx = 0;
yy = 0;
lastdir = 1; % 0 L<---R 1 L--->R
currdir = 1;
%usage = zeros(numMIs,numMIs); debug
for j=1:numMIs*numMIs
    yy = ceil((j)/numMIs);
    currdir = mod(yy, 2);
    if lastdir == currdir
        if currdir == 1
            xx = xx + 1; % 1 odd - add
        else
            xx = xx - 1; % 0 even - subtract
        end
    end
    %usage(yy,xx)=1;
    %imshow(usage,[])
    lastdir = currdir;
    A =squeeze(double(LF(yy+1,xx+1, :, :, 1:3))); % pad views
    %A(:, end+1, :) = 0; %one-pixel padding
    A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
    %clipping
    A(A<0) = 0;
    A(A>1023) = 1023;
    %rounding to integer
    A = double(uint16(A));
    if RGB44410 == 1
        PVS_RGB44410(:,:,:,j) = uint16(A);
    end
    if YUV44410 == 1
        PVS_YUV44410(:,:,:,j) = rgb2ycbcrn(double(A) / (2^10-1),10);
    end
    if YUV4208 == 1
        aux = rgb2ycbcrn(double(A) / (2^10-1),8);
        aux(:, end+1, :,:) = 0; %one-pixel padding
        PVS_YUV4208(:,j) = downsampleChroma(aux);
        clear aux
    end
end

if RGB44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_RGB44410.rgb',destinFolderRGB44410,i), 'w' );
    fwrite(fileID, permute(PVS_RGB44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_YUV44410.yuv',destinFolderYUV44410,i), 'w' );
    fwrite(fileID, permute(PVS_YUV44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV4208 == 1
	fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_YUV4208.yuv',destinFolderYUV4208,i), 'w' );
    for j=1:numMIs*numMIs
        fwrite(fileID, PVS_YUV4208{1,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{2,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{3,j}', 'uint8');
    end
    fclose(fileID);
end
end

function genPVSSERPENTINE_N(i, LF, RGB44410, YUV44410, YUV4208, numMIs)
destinFolderRGB44410 = sprintf('PVSSERPENTINE_N%d/RGB44410/I%02d',numMIs,i);
destinFolderYUV44410 = sprintf('PVSSERPENTINE_N%d/YUV44410/I%02d',numMIs,i);
destinFolderYUV4208 = sprintf('PVSSERPENTINE_N%d/YUV4208/I%02d',numMIs,i);
if RGB44410 == 1; mkdir(destinFolderRGB44410); end
if YUV44410 == 1; mkdir(destinFolderYUV44410); end
if YUV4208 == 1; mkdir(destinFolderYUV4208); end

res=size(LF);
H=res(3);
W=res(4);

% only 13x13 in this rep
PVS_RGB44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
PVS_YUV44410 = zeros(H, W, 3,numMIs * numMIs, 'uint16');
xx = 0;
yy = 0;
lastdir = 1; % 0 L<---R 1 L--->R
currdir = 1;
offset = (15 - numMIs)/2;
%usage = zeros(numMIs,numMIs); debug
for j=1:numMIs*numMIs
    yy = ceil((j)/numMIs);
    currdir = mod(yy, 2);
    if lastdir == currdir
        if currdir == 1
            xx = xx + 1; % 1 odd - add
        else
            xx = xx - 1; % 0 even - subtract
        end
    end
    %usage(yy,xx)=1;
    %imshow(usage,[])
    lastdir = currdir;
    A =squeeze(double(LF(yy+offset,xx+offset, :, :, 1:3))); % pad views
    %A(:, end+1, :) = 0; %one-pixel padding
    A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
    %clipping
    A(A<0) = 0;
    A(A>1023) = 1023;
    %rounding to integer
    A = double(uint16(A));
    if RGB44410 == 1
        PVS_RGB44410(:,:,:,j) = uint16(A);
    end
    if YUV44410 == 1
        PVS_YUV44410(:,:,:,j) = rgb2ycbcrn(double(A) / (2^10-1),10);
    end
    if YUV4208 == 1
        aux = rgb2ycbcrn(double(A) / (2^10-1),8);
        aux(:, end+1, :,:) = 0; %one-pixel padding
        PVS_YUV4208(:,j) = downsampleChroma(aux);
        clear aux
    end
end

if RGB44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_N%d_RGB44410.rgb',destinFolderRGB44410,i,numMIs), 'w' );
    fwrite(fileID, permute(PVS_RGB44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_N%d_YUV44410.yuv',destinFolderYUV44410,i,numMIs), 'w' );
    fwrite(fileID, permute(PVS_YUV44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV4208 == 1
	fileID = fopen( sprintf('%s/I%02d_PVSSERPENTINE_N%d_YUV4208.yuv',destinFolderYUV4208,i,numMIs), 'w' );
    for j=1:numMIs*numMIs
        fwrite(fileID, PVS_YUV4208{1,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{2,j}', 'uint8');
        fwrite(fileID, PVS_YUV4208{3,j}', 'uint8');
    end
    fclose(fileID);
end
end

function [out] = downsampleChroma(in)
    out = cell(3,1);
    for i=1:3
        out{i} = in(:,:,i);
    end
    for i=2:3
        tmp = imfilter(double(out{i}), [1 6 1], 'replicate', 'same');
        out{i} = tmp(:,1:2:end);
    end
    for i=2:3
        tmp = imfilter(out{i}, [0 ; 4 ; 4], 'replicate', 'same');
        out{i} = uint8((tmp(1:2:end,:) + 32) / 64);
    end
end













