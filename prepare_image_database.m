%Simulate Metal Artifact for CT
%
%     *Make sure to have (1) Deep Lesion dataset under <image_folder>
%                        (2) image_list.txt
%
% For each mode train/val, this file will generate
%    (1) mask.h5: all metal masks
%    (2) metal_trace.h5: all metal traces
%
% For each mode train/val and each CT image, this file will generate 
% an h5 file. Each h5 file includes the following attributes:
%    -gt_CT (W,H,1):
%       Monochromatic CT from Deep Lesion
%    -poly_CT (W,H,1): 
%       Polychromatic CT w/o artifact. This is the target CT during training.
%    -ma_CT (W,H,K):
%       CT images with metal artifact. K is the number of all metal masks.
%    -LI_CT (W,H,K):
%       CT reconstructed from LI sinogram
%    -poly_sinogram (W',H',1):
%       Polychromatic sinogram w/o artifact. This is the target sinogram during training.
%    -ma_sinogram (W',H',K):
%       Sinograms corrupted by metal artifact. K is the number of all metal traces.
%    -LI_sinogram (W',H',K):
%       Linear interpolated sinogram.


clear; clc

%--------------------------------------------------------------------------
% Select train/val
mode = 'train'; 
% Specify your own output folder
database_root = sprintf('database_MAR/images_%s', mode);
% Specify the location of Deep Lesion CT images
image_folder = 'images';
%--------------------------------------------------------------------------


data_list = importdata('dataset/image_list.txt');
load('dataset/CT_Samples.mat', 'CT_samples_bwMetal');
% parallel computing
poolobj = gcp('nocreate');
if isempty(poolobj)
   parpool(12) 
end

% If we want Python hdf5 matrix to have size (N x H x W), 
% Matlab matrix should have size (W x H x N) 
% Therefore, we can permute (H, W, N) to (W x H x N)

% polychromatic paramters
MARpara = getMARpara();
CTpara = getCTpara();

fprintf('Finish setting MARpara and CTpara...\n')

% image and mask indices
if strcmp(mode, 'train')
    image_indices = CTpara.train_indices;
    mask_indices = CTpara.train_mask_indices;   
else
    image_indices = CTpara.val_indices;
    mask_indices = CTpara.val_mask_indices;
end

image_size = [CTpara.imPixNum, CTpara.imPixNum, numel(mask_indices)];
sinogram_size = [CTpara.sinogram_size_x, CTpara.sinogram_size_y, numel(mask_indices)];

% metal and metal trace
selected_metal = CT_samples_bwMetal(:, :, mask_indices);
mask_all = single(zeros(image_size));
metal_trace_all = single(zeros(sinogram_size));
metal_proj_all = single(zeros(sinogram_size));

mask_name = sprintf('%s/mask.h5', database_root);
if ~exist(mask_name, 'file')
    h5create(mask_name, '/mask', image_size, 'DataType', 'single');
    fprintf('Create %s...\n', mask_name)
end

metal_trace_name = sprintf('%s/metal_trace.h5', database_root);
if ~exist(metal_trace_name, 'file')
    h5create(metal_trace_name, '/metal_trace', sinogram_size, 'DataType', 'single');
    fprintf('Create %s...\n', metal_trace_name)
end

for i = 1:numel(mask_indices)
   mask_resize = imresize(selected_metal(:, :, i), [CTpara.imPixNum, CTpara.imPixNum], 'Method', 'bilinear');
   
   mask_proj = fanbeam(mask_resize, CTpara.SOD,...
            'FanSensorGeometry', 'arc',...
            'FanSensorSpacing', CTpara.angSize, ...
            'FanRotationIncrement', 360/CTpara.angNum);
   metal_trace = single(mask_proj > 0);
   
   
   mask_all(:, :, i) = mask_resize';
   metal_trace_all(:, :, i) = metal_trace';
   
   % metal trace partial volume effect
   Pmetal_kev = mask_proj * CTpara.imPixScale;
   Pmetal_kev = MARpara.metalAtten * Pmetal_kev;
   % partial volume effect
   Pmetal_kev_bw =imerode(Pmetal_kev>0, [1 1 1]');
   Pmetal_edge = xor((Pmetal_kev>0), Pmetal_kev_bw);
   Pmetal_kev(Pmetal_edge) = Pmetal_kev(Pmetal_edge) / 4;
   metal_proj_all(:,:,i) = Pmetal_kev'; 
   
   i
end

h5write(mask_name, '/mask', mask_all);
h5write(metal_trace_name, '/metal_trace', metal_trace_all);

%% metal trace partial volume effect


%% generate data
% for i = 1:numel(image_indices)    
for i = 1 % for demo    
    fprintf('Processing %s', data_list{image_indices(i)})
    raw_image = imread([image_folder '/' data_list{image_indices(i)}]);
    
    % Deep Lesion offsets. Make sure 'image' has valid HU values after this
    % step
    image = single(raw_image) - 32768;
    
    image = imresize(image, [CTpara.imPixNum, CTpara.imPixNum], 'Method', 'bilinear');
    image(image < -1000) = -1000;
    
    % The input CT image to synArtifactMulti must have valid HU values (do not rescale to [0,1] or [0,255], etc)
    [ma_sinogram_all, LI_sinogram_all, gt_sinogram_water_all, ma_CT_all, LI_CT_all, gt_CT] = synArtifactMulti(image, CTpara, MARpara, mask_all, metal_trace_all, metal_proj_all);
    
    data_name = sprintf('%s/%04d.h5', database_root, i);
    if ~exist(data_name, 'file')
        h5create(data_name, '/ma_sinogram', sinogram_size, 'DataType', 'single');
        h5create(data_name, '/LI_sinogram', sinogram_size, 'DataType', 'single');
        h5create(data_name, '/gt_sinogram_water', sinogram_size, 'DataType', 'single');
        h5create(data_name, '/ma_CT', image_size, 'DataType', 'single');
        h5create(data_name, '/LI_CT', image_size, 'DataType', 'single');
        h5create(data_name, '/gt_CT', [CTpara.imPixNum, CTpara.imPixNum], 'DataType', 'single');
    end
    
    h5write(data_name, '/ma_sinogram', ma_sinogram_all);
    h5write(data_name, '/LI_sinogram', LI_sinogram_all);
    h5write(data_name, '/gt_sinogram_water', gt_sinogram_water_all);
    h5write(data_name, '/ma_CT', ma_CT_all);
    h5write(data_name, '/LI_CT', LI_CT_all);
    h5write(data_name, '/gt_CT', gt_CT);
    i
end
