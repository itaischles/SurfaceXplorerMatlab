% species 1 parameters
mu1 = 720; % central wavelength for SAS of species 1 (nm)
sigma1 = 40; % spectral width of SAS of species 1 (nm)
k1 = 1/60; % rate constant of species 1 (ps^-1)
c0 = 1; % initial concentration (a.u.)

% species 2 parameters
mu2 = 730; % central wavelength for SAS of species 1 (nm)
sigma2 = 40; % spectral width of SAS of species 1 (nm)
k2 = 1/2500; % rate constant of species 2 (ps^-1)

% create time and wavelength axes
TA = TAsurface('real_test.csv');
TVec = TA.TVec;
WVec = TA.WVec;
TVec_linsamp = (min(TVec):TVec(2)-TVec(1):max(TVec))';

% extended TVec for convolution
dt = TVec_linsamp(2)-TVec_linsamp(1);
TVec_linsamp = -max(TVec):TVec(2)-TVec(1):max(TVec);

% instrument response function
irf_sigma = 0.3; % (ps)
irf = 1/(sqrt(2*pi)*irf_sigma)*exp(-0.5*(TVec_linsamp/irf_sigma).^2);

% create species 1
spectrum1 = 1/(sqrt(2*pi)*sigma1)*exp(-0.5*(WVec-mu1).^2/sigma1^2);
conc1 = c0*exp(-k1*TVec_linsamp);
conc1(TVec_linsamp<0) = 0;
conc1 = conv(conc1, irf, 'same');
conc1 = interp1(TVec_linsamp, conc1, TVec);

% create species 2
spectrum2 = 1/(sqrt(2*pi)*sigma2)*exp(-0.5*(WVec-mu2).^2/sigma2^2);
conc2 = 0.75*c0*k1/(k1-k2)*(exp(-k2*TVec_linsamp) - exp(-k1*TVec_linsamp));
conc2(TVec_linsamp<0) = 0;
conc2 = conv(conc2, irf, 'same');
conc2 = interp1(TVec_linsamp, conc2, TVec);

% create concentration and spectrum vectors for all species
spectVec = [spectrum1, spectrum2];
concVec = [conc1'; conc2'];

% create deltaA spectrum + noise
deltaA = spectVec*concVec + 1e-3*randn(numel(WVec), numel(TVec));

% save matrix
matrix2save = [TVec(1:end-1)'; deltaA(:, 1:end-1)];
matrix2save = [[0; WVec], matrix2save];
filter = {'*.csv'};
[FileName, PathName] = uiputfile('*.csv', 'Save test data');
csvwrite(FileName, matrix2save);

clearvars
