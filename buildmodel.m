function S = buildmodel(varargin)
% S = buildexp(time,fcn,t0,irf,KineticParams)
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2017
% ----------------------------
if nargin~=5
  error('not enough inputs, check the function input')
end

if isstruct(varargin{1});
  data = varargin{1};
  time = data.time;
else
  time = varargin{1};
end
if length(time)<2;
  error('Provide a time vector of more than a single value')
end

model = varargin{2};
t0 = varargin{3};
irf = varargin{4};
K = varargin{5};

if ~isa(model, 'function_handle')
  error('Provide a model in function form')
end
if length(t0)~=1
  error('t0 should be a single value')
end
if length(irf)~=1
  error('IRF should be a single value')
end
te = time(:).';


S = model(K,te);

% step function at time zero
indx = (te)<0;
S(:,indx) = 0;
S = S.';

% shift t0
for i = 1:size(S,2)
  S(:,i) = interp1(te+t0,S(:,i).',te,'pchip');
end

% perform the convolution
if irf ~=0
  ind = te<(t0+100*irf);
  Sp = convIRF(S(ind,:),te(ind),irf);
  S(ind,:) = Sp;
end


return