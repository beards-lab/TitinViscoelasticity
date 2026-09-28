function u = deconvSensor(y, dt, f0, zeta, fc)
% deconvSensor  Inverse of sensorFilter.m: recover the force acting on the
% transducer from the recorded force y (uniform step dt), i.e. apply
%   H^-1(s) = (s^2 + 2 zeta w0 s + w0^2)/w0^2,  w0 = 2 pi f0,
% times a zero-phase low-pass |1/(1 + (f/fc)^8)| (default fc = 3 kHz) that
% bounds the noise gain of the inverse. Utility for DataProcessing/AnalyzeRinging.m
% (display/diagnostics only - fits should keep filtering the MODEL with
% sensorFilter.m, which does not amplify data noise).
% The transducer is a mass-spring-damper m x'' + c x' + k x = F_fibre with
% reading k x: f0 = sqrt(k/m)/(2 pi), zeta = c/(2 sqrt(k m)).
if nargin < 5 || isempty(fc), fc = 3000; end
y = y(:);
N = numel(y);
Np = 4*N;                                   % pad with the last value (hold)
fr = (0:floor(Np/2))'/(Np*dt);
w = 2*pi*fr; w0 = 2*pi*f0;
Hinv = (w0^2 - w.^2 + 2i*zeta*w0*w)/w0^2;
lp = 1./(1 + (fr/fc).^8);
X = fft([y; repmat(y(end), Np - N, 1)] - y(1));
X = X(1:numel(fr)).*Hinv.*lp;
u = real(ifft([X; conj(X(end - 1 + mod(Np, 2):-1:2))]));
u = u(1:N) + y(1);
end
