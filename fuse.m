function [ imgFused ] = fuse( imgLowPass, imgPixelwiseShifted, ThresholdRatio, T, fuse_type)
% Replacing the pixels, whose difference is greater than the threshold
% in the pixel-wise shifted image with the corresponding pixels in the low-pass image
Diff = imgLowPass - imgPixelwiseShifted;

if fuse_type == 1
    Threshold = (max(imgLowPass(:))-min(imgLowPass(:)))*ThresholdRatio;
    imgPixelwiseShifted(abs(Diff) > Threshold) = ...
        imgLowPass(abs(Diff) > Threshold);
    imgFused = imgPixelwiseShifted;
    
end

if fuse_type == 2
    NeighbourhoodWidth = ceil(4*T);
    imgTemp = imgPixelwiseShifted;
    
    ThresholdMap = GetThresholdMap( imgLowPass, ThresholdRatio, NeighbourhoodWidth );
    imgTemp(abs(Diff)>abs(ThresholdMap)) = imgLowPass(abs(Diff)>abs(ThresholdMap));
    imgFused = imgTemp;
end

end

