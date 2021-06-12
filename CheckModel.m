function [prefit] = CheckModel(data,model,waverange)

spec = data.spec;
kinmodel = model.kinmod;
params = model.initialguess;
decay = kinmodel(params);

% pull out the individual wavelength components
wavelength = data.wavelength;
time = data.time;

ind=zeros(length(waverange),1);
vp=0;
for i=1:length(waverange)
    [~,v]=min(abs(wavelength-waverange(i)));
    if vp==v
        v=v+1;
    end
    ind(i)=v;
    vp=v;
end
ind=unique(ind);
%wavel=wavelength(ind);
data=spec(:,ind);

w=(decay'*decay)\(decay'*data);
prefit=decay*w;

return



