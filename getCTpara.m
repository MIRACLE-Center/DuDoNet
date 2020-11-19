function CTpara = getCTpara()

CTpara.imPixNum = 416;
CTpara.angSize = 0.05; % degree
CTpara.angNum = 640;
CTpara.SOD = 1075;   % source to origin distance 59.5 cm in cnnmar
CTpara.imPixScale = 512 / 416 * 0.08; % cm  pixel spacing 0.08 cm in cnn-mar

CTpara.train_indices = (0:3999) * 10 + 1;
CTpara.val_indices = (0:199) * 10 + 45000;

CTpara.val_mask_indices = [1,2,20,30,36,43,63,64,98,100];
CTpara.train_mask_indices = setdiff(1:100, [1,2,20,30,36,43,63,64,98,100]);


CTpara.sinogram_size_x = 640;
CTpara.sinogram_size_y = 641;

CTpara.window = [-175, 275] / 1000 * 0.192 + 0.192;

end