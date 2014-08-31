function [ mOutputImage ] = NonLocalMeanITImage( mInputImage, localWinRadius, localWinStd, searchWinRadius, weightsStd )
% ----------------------------------------------------------------------------------------------- %
% [ mOutputImage ] = NonLocalMeanITImage( mInputImage, localWinRadius, searchWinRadius, localRadiusStd, weightStd )
%   Applies the Non Local Means Filter on the Input Image
% Input:
%   - mInputImage           -   Input image.
%                               Matrix, 1 Channels, Floating Point, [0, 1]
%   - localWinRadius        -   Local Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - localWinStd           -   Local Window Gaussian Kernel STD.
%                               Scalar, Floating Point [0.1, 20].
%   - searchWinRadius       -   Search Window Radius.
%                               Scalar, Floating Point, {1, 2, ..., 10}.
%   - weightsStd            -   Weights STD Factor.
%                               Scalar, Floating Point [0.1, 20].
% Output:
%   - mOutputImage          -   Input image.
%                               Matrix, 1 Channels, Floating Point, [0, 1]
% Remarks:
%   1.  Prefixes:
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  Classic implemntation by ''.
%   3.  Inve
% TODO:
%   1.  Implement `im2col` manually.
%   Release Notes:
%   -   1.0.000     22/08/2014  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

if ~libisloaded('simple')
	loadlibrary('./simple.dll', './simple.h');
end

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

imageNumRows = size(mInputImage, 1);
imageNumCols = size(mInputImage, 2);

mOutputImage = zeros(imageNumRows, imageNumCols);

gaussianKernelStd  = localWinStd;
vGaussianKernel    = exp(-( (-localWinRadius:localWinRadius) .^ 2) / (2 * gaussianKernelStd * gaussianKernelStd));
vGaussianKernel    = vGaussianKernel.' * vGaussianKernel;
vGaussianKernel    = vGaussianKernel(:) / sum(vGaussianKernel(:));

searchWinEffRadius = searchWinRadius + localWinRadius;

localWinSize               = (2 * localWinRadius) + 1;
localWinNumPixels          = localWinSize * localWinSize;

searchWindowNumRows   = (2 * searchWinRadius) + 1;
searchWindowNumCols   = searchWindowNumRows;
searchWindowNumPixels = searchWindowNumRows * searchWindowNumCols;

% The location of the Reference Pixel in 'mCurrSearchWindowPatches'
refPixelRowIdx     = (localWinNumPixels     + 1) / 2;
refPixelColIdx     = (searchWindowNumPixels + 1) / 2;

mCurrSearchWindowPatches = zeros(localWinNumPixels, searchWindowNumPixels);

mInputImage = padarray(mInputImage, [searchWinEffRadius searchWinEffRadius],'replicate');

for iRowIdx = (searchWinEffRadius + 1):(imageNumRows + searchWinEffRadius)
    for jColIdx = (searchWinEffRadius + 1):(imageNumCols + searchWinEffRadius)
        vSearchWindowRowIdx = (iRowIdx - searchWinEffRadius):(iRowIdx + searchWinEffRadius);
        vSearchWindowColIdx = (jColIdx - searchWinEffRadius):(jColIdx + searchWinEffRadius);
            
		mCurrSearchWindow        = mInputImage(vSearchWindowRowIdx, vSearchWindowColIdx);
%         mCurrSearchWindowPatchesA = im2col(mCurrSearchWindow, [localWinSize, localWinSize], 'sliding');
        mCurrSearchWindowPatches = Im2ColSliding(mCurrSearchWindow, [localWinSize, localWinSize]);
% 		mCurrSearchWindowPatches = calllib('simple', 'im2col_c', mCurrSearchWindowPatches,...
% 			mCurrSearchWindow, 2*searchWinEffRadius+1, 2*searchWinEffRadius+1, localWinSize);
% 		isequal(mCurrSearchWindowPatchesA, mCurrSearchWindowPatches)
        
        mOutputImage(iRowIdx - searchWinEffRadius, jColIdx - searchWinEffRadius) = NonLocaPatchFilter(mCurrSearchWindowPatches, refPixelRowIdx, refPixelColIdx, vGaussianKernel, weightsStd);
    end
end


end


% % % % % % function [ outputPx ] = NonLocaPatchFilter( mSearchWindowPatches, refPixelRowIdx, refPixelColIdx, vWeighingKernel, weightsStd )
% % % % % % % ----------------------------------------------------------------------------------------------- %
% % % % % % % [ outputPx ] = NonLocaPatchFilter( mSearchWindow, localWinRadius, vWeighingKernel, weightStd )
% % % % % % %   Applies the Non Local Means Filter on the Input Search Window per Local Window
% % % % % % % Input:
% % % % % % %   - mSearchWindow         -   Input Pixel.
% % % % % % %                               Matrix, 1 Channels, Floating Point, [0, 1]
% % % % % % %   - localWinRadius        -   Local Window Radius.
% % % % % % %                               Scalar, Floating Point, {1, 2, ..., 10}.
% % % % % % %   - vWeighingKernel       -   Local Window Gaussian Kernel.
% % % % % % %                               Vector, Floating Point (0, inf)
% % % % % % %   - weightStd             -   Weights STD Factor.
% % % % % % %                               Scalar, Floating Point [0.1, 20].
% % % % % % % Output:
% % % % % % %   - outputPx              -   Output Pixel.
% % % % % % %                               Scalar, Floating Point, [0, 1]
% % % % % % % Remarks:
% % % % % % %   1.  Prefixes:
% % % % % % %       -   'm' - Matrix.
% % % % % % %       -   'v' - Vector.
% % % % % % %   2.  Classic implemntation by ''.
% % % % % % %   3.  Inve
% % % % % % % TODO:
% % % % % % %   1.  Implement `im2col` manually.
% % % % % % %   Release Notes:
% % % % % % %   -   1.0.000     22/08/2014  Royi Avital
% % % % % % %       *   First release version.
% % % % % % % ----------------------------------------------------------------------------------------------- %
% % % % % % 
% % % % % % FALSE   = 0;
% % % % % % TRUE    = 1;
% % % % % % OFF     = 0;
% % % % % % ON      = 1;
% % % % % % 
% % % % % % vRefPatch = mSearchWindowPatches(:, refPixelColIdx);
% % % % % % 
% % % % % % vWeightedPacthDistnace = bsxfun(@minus, mSearchWindowPatches, vRefPatch);
% % % % % % vWeightedPacthDistnace = bsxfun(@times, vWeightedPacthDistnace, vWeighingKernel);
% % % % % % % vWeightedPacthDistnacea = sum([vWeightedPacthDistnace .^ 2]);
% % % % % % vWeightedPacthDistnace = sum([vWeightedPacthDistnace .* vWeightedPacthDistnace]);
% % % % % % vWeightedPacthDistnace = exp(-vWeightedPacthDistnace / (weightsStd * weightsStd));
% % % % % % 
% % % % % % % Cancelling the reference pixel
% % % % % % vWeightedPacthDistnace(refPixelColIdx) = 0;
% % % % % % 
% % % % % % outputPx = sum(vWeightedPacthDistnace .* mSearchWindowPatches(refPixelRowIdx, :)) ./ sum(vWeightedPacthDistnace);
% % % % % % 
% % % % % % 
% % % % % % end

