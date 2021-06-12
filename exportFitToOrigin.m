function exportFitToOrigin(fit, fitpar_labels)

[name,loc] = uiputfile('*.csv');

if name==0
    return;
end

wavelength = fit.wavelength;
sfit = fit.SpectraFit;
sfit = [wavelength(:), sfit.'];

fitpar = fit.FitParams;
fiterr = fit.Std';
time = fit.time;
range = fit.waverange;
kfit = fit.KineticFit;
kfit = [0,range(:)';time(:),kfit];

kin = fit.KineticData;
kdat = [0,range(:)';time(:) kin];


[~,name,~]=fileparts(name);    

csvwrite(fullfile(loc,[name,'_spectra.csv']),sfit);
csvwrite(fullfile(loc,[name,'_kinfit.csv']),kfit);
csvwrite(fullfile(loc,[name,'_kindat.csv']),kdat);

formatspec = '%s %.8f +- %.8f \n';
fid = fopen(fullfile(loc,[name,'_fitpar.txt']),'wt');
for i=1:numel(fitpar_labels)
    fprintf(fid, formatspec, fitpar_labels{i}, fitpar(i), fiterr(i));
end
fclose(fid);

return