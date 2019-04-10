addpath('../4_4_4Version');

% Profile definition for easier usage
lf_path = {'/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/I01_Bikes'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/I02_Danger_de_Mort'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/I04_Stone_Pillars_Outside'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/I09_Fountain_Vincent_2'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/set2'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/poznanlab1'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/tarot'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/greek'
           '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/sideboard'};

org_bps = [10, 10, 10, 10, 10, 8, 8, 8, 8];
ppm_bps = [16, 16, 16, 16, 16, 8, 8, 8, 8];

target_bpp_small = [0.001, 0.005, 0.02, 0.1, 0.75];
target_bpp_big = [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1];

results = cell(9, 6, 6); % files, target bpp, results  

if isempty(gcp('nocreate'))
    myCluster = parcluster('local');
    myCluster.NumWorkers = 9;
    parpool(myCluster, 9)
end

parfor file = 1:9
    org_bps = [10, 10, 10, 10, 10, 8, 8, 8, 8];
    ppm_bps = [16, 16, 16, 16, 16, 8, 8, 8, 8];
    
    cur_lf = erase(lf_path{file}, '/nfs/home/lathomaz.it/JPEG_Pleno_CE6.2/WASP/DATASETS/');
    fprintf('Working on file %s: ', cur_lf);

    no_points = length(dir(strcat('./CONFIGURATIONFILES/', cur_lf, '/*.json')));
    
    if no_points == 5
        target_bpp = [target_bpp_small target_bpp_small(end)];
    else
        target_bpp = target_bpp_big;
    end

%     for p = 1:no_points
    for p = 1:6
        if p <= no_points
            work_path = strcat('./results/', cur_lf, '_', num2str(target_bpp(p)));

            aux_results = cell(1, 6);
            aux_results{1} = cur_lf;
            aux_results{2} = target_bpp(p);

            rec_LF = dir(strcat(work_path, '/output.LF'));
            aux_results{3} = rec_LF.bytes * 8;

            ref_LF = dir(strcat(lf_path{file}, '/*.ppm'));
            rec_LF = dir(strcat(work_path, '/decoded/PPM/*.ppm'));

            no_views = length(rec_LF);
            dims = size(imread(strcat(lf_path{file}, '/', ref_LF(1).name)));

            aux_results{3} = aux_results{3} / dims(1) / dims(2) / no_views;
            fprintf('%f (%4.3f) ', target_bpp(p), aux_results{3});

            aux_results{4} = 0;
            aux_results{5} = 0;
            aux_results{6} = 0;

            for v = 1:no_views
                im1 = imread(strcat(lf_path{file}, '/', rec_LF(v).name));
                im2 = imread(strcat(work_path, '/decoded/PPM/', rec_LF(v).name));
                
                ref_img = bitshift(im1,-6); 
                rec_img = bitshift(im2,-6);

                [PSNR_Y, PSNR_YUV, SSIM_Y] = QM(ref_img, rec_img, ppm_bps(p), org_bps(p));

                aux_results{4} = aux_results{4} + PSNR_Y;
                aux_results{5} = aux_results{5} + PSNR_YUV;
                aux_results{6} = aux_results{6} + SSIM_Y;
            end

            aux_results{4} = aux_results{4} / no_views;
            aux_results{5} = aux_results{5} / no_views;
            aux_results{6} = aux_results{6} / no_views;
        else
            aux_results = cell(1, 6);
            
            for i = 1:6
                aux_results{i} = nan;
            end
        end

        [results{file, p, :}] = aux_results{1, 1:6};
    end

    fprintf('\n');
end

results = reshape(permute(results, [2 1 3]), [], 6);
results(flipud(find(cell2mat(cellfun(@(x)any(isnan(x)), results(:, 1), 'UniformOutput', false)))), :) = [];

writetable(cell2table(results, 'VariableNames', {'LF', 'Target_bpp', 'bpp', 'PNSR_Y', 'PSNR_YUV', 'SSIM_Y'}), 'results_aux.csv');
