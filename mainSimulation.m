close all
clear
clc

nSubpixelShifts = 10; % Denoted by k in th paper
SamplingRatio = 0.1;
ImgFusionThreshold = 0.05; % Denoted by q in the paper
Wnd = [0:1];

InputImgPath = 'girls_384x256.png';
% InputImgPath = 'USAF1951_256x256.png';
% InputImgPath = 'flowers_384x288.jpg';
% InputImgPath = 'babala_256x256.bmp';

DEBUG_MODE = 0;
RGB_MODE = 0;

InputImg = im2double(imread(InputImgPath));
[mRow, nCol, nBand] = size(InputImg);

if nBand == 3 && ~RGB_MODE
    InputImg = rgb2gray(InputImg);
end

InputImgFT = fftshift(fft2(InputImg));

if DEBUG_MODE
    figure;
    ax1 = subplot(1,2,1); imagesc(InputImg); colormap(ax1, 'gray'); axis image;
    ax2 = subplot(1,2,2); imagesc(log(1+abs(InputImgFT))); colormap(ax2, 'jet'); axis image;
end
%% Simulate under-sampled Fourier spectrum through a low-pass filter
LowPassFilter = zeros(mRow, nCol);
SpiralMat = getSpiralMat(mRow,nCol);
LowPassFilter(SpiralMat <= round(mRow * nCol * SamplingRatio)) = 1;

InputImgFTUnderSampled = InputImgFT .* LowPassFilter;
InputImgUnderSampled = real(ifft2(ifftshift(InputImgFTUnderSampled)));

if DEBUG_MODE
    figure;
    ax1 = subplot(1,2,1); imagesc(InputImgUnderSampled); colormap(ax1, 'gray'); axis image;
    ax2 = subplot(1,2,2); imagesc(log(1+abs(InputImgFTUnderSampled))); colormap(ax2, 'jet'); axis image;
end

SSIM_Original = ssim(InputImg, InputImgUnderSampled);
PSNR_Original = psnr(InputImg, InputImgUnderSampled);

tic;
%% Sub-pixel shifting amount estimation
tic;
xArr = [1:nCol];
xMax = max(xArr);

for iDimension = 1:2
    if iDimension == 1
        InputSignal2D = InputImgUnderSampled;
    else
        InputSignal2D = Signal2DShifted';
    end
    
    [ OptimizedShift, Signal2DShifted ] = estimateSubpixelShift( InputSignal2D, SamplingRatio, ...
        nSubpixelShifts, Wnd, DEBUG_MODE );
    
    if iDimension == 2
        Signal2DShifted = Signal2DShifted';
    end
    
    if DEBUG_MODE
        figure; imagesc(Signal2DShifted); axis image; colormap gray;
    end
end
ImgSubpixelShifted = Signal2DShifted;

if DEBUG_MODE
    figure; imagesc(ImgSubpixelShifted); axis image; colormap gray;
end

%% Image fusion
ImgFused = fuse( InputImgUnderSampled, ImgSubpixelShifted, ImgFusionThreshold, getRingingPeriod(SamplingRatio), 1 );
if DEBUG_MODE
    figure, imagesc(ImgFused); axis image; colormap gray;
end
toc;
SSIM_Deringinged = ssim(InputImg, ImgFused);
PSNR_Deringinged = psnr(InputImg, ImgFused);

%% Results display
figure;
subplot(2,2,1); imagesc(InputImgUnderSampled); 
TitleStr = sprintf('[Original]  SSIM = %0.2f, PSNR = %0.2f', SSIM_Original, PSNR_Original);
axis image; colormap gray; title(TitleStr);

subplot(2,2,2); imagesc(ImgFused); axis image; colormap gray;
TitleStr = sprintf('[Deringinged] SSIM = %0.2f, PSNR = %0.2f', SSIM_Deringinged, PSNR_Deringinged);
axis image; colormap gray; title(TitleStr);

%% Comparison (decovolution)
de = 1;
PSF = fspecial('gaussian',10,1.5);
ImgDecovWnr = deconvwnr(InputImgUnderSampled, PSF, de);
SSIM_DecovWnr = ssim(InputImg, ImgDecovWnr);
PSNR_DecovWnr = psnr(InputImg, ImgDecovWnr);
TitleStr = sprintf('[DecovWnr] SSIM = %0.2f, PSNR = %0.2f', SSIM_DecovWnr, PSNR_DecovWnr);
subplot(2,2,3); imagesc(ImgDecovWnr); 
axis image; colormap gray; title(TitleStr);

ImgDecovBlind = deconvblind(InputImgUnderSampled, PSF, de);
SSIM_DecovBlind = ssim(InputImg, ImgDecovBlind);
PSNR_DecovBlind = psnr(InputImg, ImgDecovBlind);
TitleStr = sprintf('[DecovBlind] SSIM = %0.2f, PSNR = %0.2f', SSIM_DecovBlind, PSNR_DecovBlind);
subplot(2,2,4); imagesc(ImgDecovBlind); 
axis image; colormap gray; title(TitleStr);
