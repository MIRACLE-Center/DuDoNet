function [ma_sinogram_all, LI_sinogram_all, gt_sinogram_water_all, ma_CT_all, LI_CT_all, gt_CT] = synArtifactMulti(imgCT, CTpara, MARpara, mask_all, metal_trace_all, metal_proj_all)

% If we want Python hdf5 matrix to have size (N x H x W), 
% Matlab matrix should have size (W x H x N) 
% Therefore, we can permute (H, W, N) to (W x H x N)

n_mask = size(mask_all, 3);

%% tissue composition
MiuWater = MARpara.MiuWater;
threshWater = MARpara.threshWater;
threshBone = MARpara.threshBone;

img = imgCT / 1000 * MiuWater + MiuWater;
gt_CT = img';
imgWater = zeros(size(img));
imgBone = zeros(size(img));
bwWater = img <= threshWater;
bwBone = img >= threshBone;
bwBoth = im2bw(1 - bwWater - bwBone, 0.5);
imgWater(bwWater) = img(bwWater);
imgBone(bwBone) = img(bwBone);
imgBone(bwBoth) = (img(bwBoth) - threshWater) ./ (threshBone - threshWater) .* img(bwBoth);
imgWater(bwBoth) = img(bwBoth) - imgBone(bwBoth);


%% Metal
ma_sinogram_all = single(zeros(CTpara.sinogram_size_x, CTpara.sinogram_size_y, n_mask));
LI_sinogram_all = single(zeros(CTpara.sinogram_size_x, CTpara.sinogram_size_y, n_mask));
ma_CT_all = single(zeros(CTpara.imPixNum, CTpara.imPixNum, n_mask));
LI_CT_all = single(zeros(CTpara.imPixNum, CTpara.imPixNum, n_mask));
gt_sinogram_water_all = single(zeros(CTpara.sinogram_size_x, CTpara.sinogram_size_y, n_mask));
scatterPhoton = 20;
parfor i = 1:n_mask
    imgMetal = squeeze(mask_all(:, :, i))';
    metal_trace = squeeze(metal_trace_all(:, :, i))';
    Pmetal_kev = squeeze(metal_proj_all(:, :, i))';
    
    bwMetal = im2bw(imgMetal);
    imgWater_local = imgWater;
    imgBone_local = imgBone;
    imgWater_local(bwMetal) = 0;
    imgBone_local(bwMetal) = 0;
    % insert metal into water and bone image
    Pwater_kev = fanbeam(imgWater_local, CTpara.SOD,...
        'FanSensorGeometry', 'arc',...
        'FanSensorSpacing', CTpara.angSize, ...
        'FanRotationIncrement', 360/CTpara.angNum);
    Pwater_kev = Pwater_kev * CTpara.imPixScale;

    Pbone_kev = fanbeam(imgBone_local, CTpara.SOD,...
            'FanSensorGeometry', 'arc',...
            'FanSensorSpacing', CTpara.angSize, ...
            'FanRotationIncrement', 360/CTpara.angNum);
    Pbone_kev = Pbone_kev * CTpara.imPixScale;

    % sinogram with metal
    projkevAllLocal = zeros(CTpara.sinogram_size_y, CTpara.sinogram_size_x, 3);
    projkevAllLocal(:, :, 1) = Pwater_kev;
    projkevAllLocal(:, :, 2) = Pbone_kev;
    projkevAllLocal(:, :, 3) = Pmetal_kev;
    projkvpMetal = pkev2kvp(projkevAllLocal, MARpara.spectrum, MARpara.energies, MARpara.kev, MARpara.MiuAll);
    temp = round(exp(-projkvpMetal) .* MARpara.photonNum);
    % poisson
    temp = temp + scatterPhoton;
    ProjPhoton = poissrnd(temp);
    ProjPhoton(ProjPhoton == 0) = 1;
    projkvpMetalNoise = -log(ProjPhoton ./ MARpara.photonNum);

    % correction
    p1 = reshape(projkvpMetalNoise, CTpara.sinogram_size_y*CTpara.sinogram_size_x, 1);
    p1BHC = [p1  p1.^2  p1.^3] * MARpara.paraBHC;
    ma_sinogram = reshape(p1BHC, CTpara.sinogram_size_y, CTpara.sinogram_size_x);
    LI_sinogram = projInterp(ma_sinogram, metal_trace);

    % reconstruct   
    ma_CT = ifanbeam(ma_sinogram, CTpara.SOD,...
            'FanSensorGeometry', 'arc',...
            'FanSensorSpacing', CTpara.angSize,...
            'OutputSize', CTpara.imPixNum,...
            'FanRotationIncrement', 360 / CTpara.angNum);
    ma_CT = ma_CT / CTpara.imPixScale;
    
    LI_CT = ifanbeam(LI_sinogram, CTpara.SOD,...
            'FanSensorGeometry', 'arc',...
            'FanSensorSpacing', CTpara.angSize,...
            'OutputSize', CTpara.imPixNum,...
            'FanRotationIncrement', 360 / CTpara.angNum);
    LI_CT = LI_CT / CTpara.imPixScale;
    
    ma_sinogram_all(:, :, i) = ma_sinogram';
    LI_sinogram_all(:, :, i) = LI_sinogram';
    ma_CT_all(:, :, i) = ma_CT';
    LI_CT_all(:, :, i) = LI_CT';
    
    % sinogram groundtruth, metal pixel is set to water
    ct_gt_water = gt_CT;
    ct_gt_water(bwMetal) = MiuWater;
    gt_sinogram_water = fanbeam(ct_gt_water, CTpara.SOD,...
        'FanSensorGeometry', 'arc',...
        'FanSensorSpacing', CTpara.angSize, ...
        'FanRotationIncrement', 360/CTpara.angNum)* CTpara.imPixScale;
    gt_sinogram_water_all(:,:,i) = gt_sinogram_water';

    

end


end