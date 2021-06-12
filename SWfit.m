function fit=SWfit(varargin)
% [params,error,kfit,kdat]=SWfit(data,waverange,model,algorithm)
% Globally fits a single or range of wavelengths using the given kinetic
% model.
%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------

if nargin<3
  error('Not enough input parameters');
end
if isstruct(varargin{1});
  data = varargin{1};
  fit = data;
  spec = data.spec;
else
  error('Provide the data structure as the first input')
end


model = varargin{3};
if ~isstruct(model)
  error('Provide the model structure as the third input')
end

% % % model = fixmodel(model);
guess = model.initialguess;

algorithm=varargin{4};

% pull out the individual wavelength components
waverange = varargin{2}(:);
wavelength = data.wavelength;

ind = zeros(length(waverange),1);
vp = 0;
for i=1:length(waverange)
  [~,v]=min(abs(wavelength-waverange(i)));
  if vp==v
    v = v+1;
  end
  ind(i) = v;
  vp = v;
end
ind = unique(ind);
wavel = wavelength(ind);
data = spec(:,ind);

if strncmp(algorithm,'l',1)
  res = @(params)fitkin(params,model,data);
  
  [bestfit, ~, iter] = LMFnlsq(res,guess,'MaxIter', 5000,'FunTol',1e-10,'XTol',1e-10);
  if iter==1 || iter ==-100
    disp('Bad initial guess or poor model');return
  end
elseif strncmp(algorithm,'s',1)
  res=@(params)fitkinSSE(params,model,data);
  
  options=optimset('MaxFunEvals',1000, 'MaxIter', 5000,'TolFun',1e-10,'TolX',1e-10);

  [bestfit,~,flag] = fminsearchbnd(res,guess,model.lb,model.ub,options);
  if flag<1
    disp('Bad initial guess or poor model');return
  end
end

% % % bestfit = model.normfun(bestfit);
% sort of a work around at the moment for dealing with the constraints
% % % if isfield(model,'lb') && isfield(model,'ub')
% % %   model = rmfield(model,['lb'; 'ub']);
% % % end
% % % res=@(params)fitkin(params,model,data);

% Calculate the fit based on minimum obtained through fitting.
[~,w,decay] = res(bestfit);
FUN = @(params)fitkin(params,model,data);
[mse,std] = lmstats(FUN,bestfit);

sf = (decay'*decay)\(decay'*spec);

kfit = decay*w;
residual = data-kfit;

fit.FitParams = bestfit;
fit.KineticFit = kfit;
fit.KineticData = data;
fit.SpectraFit = sf;
fit.waverange = wavel;
fit.Residual = residual;
fit.MeanSquareError = mse;
fit.Std = std;
fit.Population = decay;
fit.ResidualMap = sqrt(reshape(fitkin(bestfit,model,spec),size(spec)).^2);

return

function model = fixmodel(model)
if isfield(model,'lb') && isfield(model,'ub')
  lb = model.lb;
  ub = model.ub;
  guess = model.initialguess;
  % rescale the guess to take into account the bounds
  x = [lb;ub];
  range = max(x) - min(x);
  mx = mean(x);
  model.initialguess = 2*(guess - mx )./range;
  model.normfun = @(y) range(:).*y(:)/2+mx(:);
else
  model.normfun = @(y) y(:);
end
return

function [res,w,decay]=fitkin(params,model,spec)

% % % if isfield(model,'lb') && isfield(model,'ub')
% % %   ind = abs(params)>1;
% % %   penalty = sum(params(ind).^2);
% % %   params = model.normfun(params);
% % % else
% % %   penalty = 0;
% % % end

kinmodel = model.kinmod;
decay = kinmodel(params);

% Linear combination of the decays traces used to fit the decay
w=(decay'*decay)\(decay'*spec);
fit=decay*w;
% calculate the residual
res=spec-fit;

% For Levenberg-Maquardt optimization we need a column of the residuals to
% to calculate the jacobian
res=res(:);

return

function [res,w,decay]=fitkinSSE(params,model,spec)

% % % if isfield(model,'lb') && isfield(model,'ub')
% % %   ind = abs(params)>1;
% % %   penalty = sum(params(ind).^2);
% % %   params = model.normfun(params);
% % % else
% % %   penalty = 0;
% % % end

kinmodel = model.kinmod;
decay = kinmodel(params);

% Linear combination of the decays traces used to fit the decay
w=(decay'*decay)\(decay'*spec);
fit=decay*w;
% calculate the residual
res=spec-fit;

% Sum of Squares of residuals
res=res(:)'*res(:);

return