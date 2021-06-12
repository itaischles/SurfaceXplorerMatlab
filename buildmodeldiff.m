function S = buildmodeldiff(varargin)
% S = builddiffmodel(time,diffeqn,t0,irf,KineticParams,population)
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2017
% ----------------------------
if nargin~=6
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
dt = varargin{2};
t0 = varargin{3};
irf = varargin{4};
K = varargin{5};
conc = varargin{6}(:);

if ~isa(dt, 'function_handle')
  error('Provide a model in function form')
end
if length(t0)~=1
  error('t0 should be a single value')
end
if length(irf)~=1
  error('IRF should be a single value')
end
time = time(:);
tdim = length(time);
sdim = length(conc);
conc = conc/max(conc); 


te = [0;time;];
% need to add a zero and resort the data since we technically start from
% zero and it may not be there.
te = sort(te);
te = unique(te);
    
% setup the differential equations
ind = (te>=0);
odefun = @(t,a)dt(a,K,t);
options = odeset('NonNegative', 1:sdim);
S = zeros(tdim,sdim);
% solver, depending on the difficult may need to swap it to a different one
% [~,S(ind,:)] = ode45(odefun,te(ind),conc,options);
[~,S(ind,:)] = ode15s(odefun,te(ind),conc,options);
if length(te)~=tdim
    in = te~=0;
    S = S(in,:);
    te = te(in);
end

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
