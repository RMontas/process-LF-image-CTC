% 10 bpp / 13x13 SAIs

%render SAIs from RAW LFs
clear all
close all

%% Uses 6 layers of scalability
% L1 - 1 / 1
% L2 - 8 / 9
% L3 - 12 / 21
% L4 - 24 / 45
% L5 - 36 / 81
% L6 - 88 / 169

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

%% config flags
compareGeneratedWithOriginals = 0;
%% color flags
RGB44410 = 1;
YUV44410 = 1;
YUV4208 = 1;
%% representation type flag
REPPVSSPIRAL_SCL = 1; %% YUV
REPSAI = 0;
REPMI = 0;

mkdir('PVSSPIRAL_SCL')
mkdir('PVSSPIRAL_SCL/RGB44410')
mkdir('PVSSPIRAL_SCL/YUV44410')
mkdir('PVSSPIRAL_SCL/YUV4208')

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
    %%
    if REPSAI == 1
        %genSAI(i, LFCol, compareGeneratedWithOriginals, RGB44410, YUV44410, YUV4208);
    end
    %%
    if REPPVSSPIRAL_SCL == 1
        genPVSSPIRAL(i, LFCol, RGB44410, YUV44410, YUV4208, num_Layers, layerMask);
    end
    %%
    if REPMI == 1
        %genMI(i, LFCol, RGB44410, YUV44410, YUV4208);
    end
end

function genPVSSPIRAL(i, LF, RGB44410, YUV44410, YUV4208, num_Layers, layerMask)
destinFolderRGB44410 = sprintf('PVSSPIRAL_SCL/RGB44410/I%02d',i);
destinFolderYUV44410 = sprintf('PVSSPIRAL_SCL/YUV44410/I%02d',i);
destinFolderYUV4208 = sprintf('PVSSPIRAL_SCL/YUV4208/I%02d',i);
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
frame = 1;
for l = 1:num_Layers
    for j=1:numMIs*numMIs
        [yy, xx] = find(cc_spiral == j);
        if layerMask(yy,xx) == l % working for 13x13
            A =squeeze(double(LF(yy+1,xx+1, :, :, 1:3))); % pad views
            %A(:, end+1, :) = 0; %one-pixel padding
            A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
            %clipping
            A(A<0) = 0;
            A(A>1023) = 1023;
            %rounding to integer
            A = double(uint16(A));
            if RGB44410 == 1
                PVS_RGB44410(:,:,:,frame) = uint16(A);
            end
            if YUV44410 == 1
                PVS_YUV44410(:,:,:,frame) = rgb2ycbcrn(double(A) / (2^10-1),10);
            end
            if YUV4208 == 1
                aux = rgb2ycbcrn(double(A) / (2^10-1),8);
                aux(:, end+1, :,:) = 0; %one-pixel padding
                PVS_YUV4208(:,frame) = downsampleChroma(aux);
                clear aux
            end
            disp([frame yy xx])
            frame = frame + 1;
        end
    end
end

if RGB44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_SCL_RGB44410.rgb',destinFolderRGB44410,i), 'w' );
    fwrite(fileID, permute(PVS_RGB44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV44410 == 1
    fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_SCL_YUV44410.yuv',destinFolderYUV44410,i), 'w' );
    fwrite(fileID, permute(PVS_YUV44410, [2 1 3 4]), 'uint16');
    fclose(fileID);
end
if YUV4208 == 1
	fileID = fopen( sprintf('%s/I%02d_PVSSPIRAL_SCL_YUV4208.yuv',destinFolderYUV4208,i), 'w' );
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













