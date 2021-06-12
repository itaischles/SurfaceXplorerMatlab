function S = buildkineticmodel(varargin)
% S = buildkineticmodel(time,t0,irf,Kmatrix,population)
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2017
% ----------------------------
if nargin<4
  error('not enough inputs, check the function input')
end

if isstruct(varargin{1})
    data = varargin{1};
    time = data.time;
else
    time = varargin{1};
end
if length(time)<2;
  error('Provide a time vector of more than a single value')  
end
time = time(:);

t0 = varargin{2};
irf = varargin{3};
K = varargin{4};

if length(t0)~=1
  error('t0 should be a single value')
end
if length(irf)~=1
  error('IRF should be a single value')
end
if size(K,1)~=size(K,2)
  error('The kinetic matrix needs to be a square matrix')
end

sdim = length(K);
tdim = length(time);

% Concentration vector if its not provided function will use 1 for the 
% first component and zeros for the rest. 
if nargin>4
    conc = varargin{5}(:);
    if length(conc)~=length(K)
      error('Starting populations need to be the same length as the kinetic matrix')
    end
    ind = conc<0;
    conc(ind) = 0;
    conc = conc/max(conc);
else
    conc = [1;zeros(sdim-1,1)];
end

te = time;

% Solve the differential equation based on the rate matrix K
[P,V] = eig(K);
Pr = P\conc;

% due to the form of the kinetic function we can't vectorize the time
% calculation and have to perform the calculation stepwise good candidate
% to parallelize 
S = zeros(sdim,tdim);
indx = find(te>=0,1);
for i=indx:tdim
    S(:,i) = P*expm(V*te(i))*Pr;
end
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

if nargin>5
  ncom = varargin{6};
  S = S(:,ncom);
end
return
