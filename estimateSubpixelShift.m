function [ OptimizedShift, SignalShifted ] = estimateSubpixelShift( Signal, SamplingRatio, UpsamplingFactor, Wnd, DEBUG_MODE )
[mRow, nCol] = size(Signal);               
xArr = [1:nCol];        % Original sampling grid
xMax = max(xArr);

SignalFT = fft(Signal, [], 2);
fxArr = ifftshift(linspace(-0.5, 0.5, nCol));
fxGrid = repmat(fxArr, [mRow 1]);

RingingPeriod = getRingingPeriod(SamplingRatio);

ShiftAmountArr = linspace(-RingingPeriod/4, RingingPeriod/4, UpsamplingFactor);
nShift = length(ShiftAmountArr);
ShiftedSignalSet = zeros(mRow, nCol, nShift);

if DEBUG_MODE
    figure; 
end

for iShift = 1:nShift
    ShiftInSpatialDomain = ShiftAmountArr(iShift);
    PhaseShiftInFourierDomain  = 2*pi*fxGrid * ShiftInSpatialDomain;
    SignalShifted = abs(ifft(SignalFT .* exp(1j*PhaseShiftInFourierDomain), [], 2));
    ShiftedSignalSet(:, :, iShift) = SignalShifted;

    if DEBUG_MODE
        imagesc(SignalShifted); axis image; colormap gray; pause(1);
    end
end

xGridLocalLeftOffset = round(-fliplr(Wnd) * RingingPeriod/2);
xGridLocalRightOffset = round(Wnd * RingingPeriod/2);

LocalVariationMat = zeros(mRow, nCol, nShift, 2);

for Side = [1, 2]
    for iShift = 1:nShift
        for iCol = 1:nCol
            x = xArr(iCol);
            
            if Side == 1
                OffsetArr = xGridLocalLeftOffset;
            else
                OffsetArr = xGridLocalRightOffset;
            end
            
            xGridLocal = x + OffsetArr;
            xGridLocal(xGridLocal <= 0) = ...
                xGridLocal(xGridLocal <= 0) + xMax; % Equivalent to circular shift
            xGridLocal(xGridLocal >  xMax) = ...
                xGridLocal(xGridLocal > xMax) - xMax; % Equivalent to circular shift
            
            SignalLocal = ShiftedSignalSet(:, xGridLocal, iShift);
            LocalVariationMat(:, iCol, iShift, Side) = sum(abs(diff(SignalLocal, 1, 2)), 2);
        end
    end
end

[LocalVariationMinMat] = min(LocalVariationMat, [], 4);
[LocalVariationMinVal, OptimizedShiftIndexMat] = min(LocalVariationMinMat, [], 3);

SignalShifted = zeros(mRow,nCol);
for iRow = 1:mRow
    for jCol = 1:nCol
        SignalShifted(iRow, jCol) = ShiftedSignalSet(iRow, jCol, OptimizedShiftIndexMat(iRow,jCol));
    end
end

OptimizedShift = ShiftAmountArr(OptimizedShiftIndexMat);

return