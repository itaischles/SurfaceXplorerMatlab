function conS = convIRF(S,te,irf)
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2017
% ----------------------------

tdim = length(te);
sdim = size(S);
te = te(:);
% ensure we're working with columns
if sdim(1)==tdim
  sdim = sdim(2);
else
  sdim = sdim(1);
  S = S.';
end

% generate new time axis for interpolation
dt = min(te(2:end)-te(1:end-1));
range = max(te)-min(te);

% work in powers of two to help avoid issues that could arise in the FFT
nt = ceil(range/dt);
nt = pow2(nextpow2(nt));

t = (0:dt:nt*dt-dt).'+min(te);

% setup the guassian 
% work in an integer base to avoid dealing with the t0 and negative times
% and such

n = (0:2*nt-1).';
mid = round((length(n))/2)+1;

irf = irf/dt;
b = 4*log(2)/irf^2;
rf = exp(-b.*(n-n(mid)).^2);
% remove values below machine precision not sure how much this helps us in the end 
ind = abs(rf)<eps;
rf(ind) = 0;
% normalize the gaussian
rf = rf/sum(rf);

rfi = ifft(fftshift(rf));

conS = zeros(tdim,sdim);
Si =  zeros(length(n),sdim);
con = zeros(length(n),sdim);

for i=1:sdim
  % interpolate to the new time base
  Si(1:length(t),i) = interp1(te,S(:,i).',t,'pchip');
  % convolute
  con(:,i) = length(rfi)*real(fft(ifft(Si(:,i),length(n)).*rfi));
  % interpolate back to the original time base
  conS(:,i) = interp1(t,con(1:length(t),i),te,'pchip');
end

return