% ==============================
% INPUT SECTION
% ==============================
imgFile = 'test_img.png';

calibValue = [0.5 1e-6]; % [y, x], MeV/pixel for y and m/pixel for x 

ifRotateImg = false;

cropRect = [27 1 71 100]; % [y1 x1 y2 x2]

baseline = 0.0;
filterSpan = [3 3];

% settings of beamline
magGrad = [199.2 231.0 284.4]; % [T/m]
magLength = [0.0149 0.0352 0.0196]; % [m]
driftLength = [0.0399109, 0.0354887, 0.0141347, 0.740766]; % [m]

% beam energy at center
centerEne = 100; % [MeV]
centerIndex = 23;
specDir = 1;

dispersion = 'x';

% ==============================
% END OF INPUT SECTION
% ==============================

% read image
imgOrg = imread( imgFile );
if ifRotateImg
	imgOrg = imgOrg';
end
img = imgOrg( cropRect(1):cropRect(3), cropRect(2):cropRect(4) );

img = filterImg( img, 'median value', 'hsize', filterSpan );
img = filterImg( img, 'smooth', 'hsize', filterSpan );
img = filterImg( img, 'base substract', 'baselinetype', 'relative', 'baseline', baseline );

imagesc( img )

nSamples = size( img, 1 );
nPoints = size( img, 2 );
sigma = zeros( nSamples, 1 );
xgrid = (1:nPoints) * calibValue(2);

for ii = 1:nSamples
	lineout = double( img(ii,:) );
	[~, sigma(ii)] = getProfilePosWidth( xgrid * 1e3, lineout, 'type', 'rms gaussian fitting' );
end
sigma = sigma * 1e-3;

%
ene = centerEne + specDir * ((1:nSamples)-centerIndex) * calibValue(1);
matrixCS = zeros( nSamples, 3 );
for ii = 1:nSamples
    switch dispersion
        case 'x'
            [~, matrixTransport] = getTransportMatrix( driftLength, magLength, magGrad, ene(ii) );
        case 'y'
            [matrixTransport, ~] = getTransportMatrix( driftLength, magLength, magGrad, ene(ii) );
    end
    c = matrixTransport(1,1);
    s = matrixTransport(1,2);
    matrixCS(ii, 1) = c^2;
    matrixCS(ii, 2) = s^2;
    matrixCS(ii, 3) = -2*c*s;
end
vectorSigma2 = sigma.^2;
matrixWeight = diag( vectorSigma2.^-1 / sum( vectorSigma2.^-1 ) );

matrixLeft = matrixCS' * matrixWeight * matrixCS;
matrixRight = matrixCS' * matrixWeight * vectorSigma2;
twiss = matrixLeft \ matrixRight;

emittance = sqrt( twiss(1) * twiss(2) - twiss(3)^2 );

%%
vectorSigma2Fit = zeros( nSamples, 1 );

for ii = 1:nSamples
    switch dispersion
        case 'x'
            [~, matrixTransport] = getTransportMatrix( driftLength, magLength, magGrad, ene(ii) );
        case 'y'
            [matrixTransport, ~] = getTransportMatrix( driftLength, magLength, magGrad, ene(ii) );
    end
    c = matrixTransport(1,1);
    s = matrixTransport(1,2);
    vectorSigma2Fit(ii) = [c^2 s^2 -2*c*s] * twiss;
end

figure;
plot( ene, vectorSigma2, 'o' ); hold on
plot( ene, vectorSigma2Fit, '--' )