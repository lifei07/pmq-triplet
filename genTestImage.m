clear

twiss = [5e-6/1e-2 1e-2/5e-6 0];
eneStart = 90;
eneEnd = 110;
nSamples = 51;
nPoints = 100;

% settings of beamline
magGrad = [199.2 231.0 284.4]; % [T/m]
magLength = [0.0149 0.0352 0.0196]; % [m]
driftLength = [0.0399109, 0.0354887, 0.0141347, 0.740766]; % [m]

sigma2 = zeros( nSamples, 1 );
ene = linspace( eneStart, eneEnd, nSamples );

for ii = 1:nSamples
    matrixTransport = getTransportMatrix( driftLength, magLength, magGrad, ene(ii) );
    c = matrixTransport(1,1);
    s = matrixTransport(1,2);
    sigma2(ii) = [c^2 s^2 -2*c*s] * twiss';
end

plot(sigma2.^0.5);

img = zeros( nSamples, nPoints );