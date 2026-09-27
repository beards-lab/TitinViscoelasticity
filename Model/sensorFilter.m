function y = sensorFilter(u, dt, f0, zeta)
% sensorFilter  Second-order (underdamped) force-transducer response
%   H(s) = w0^2/(s^2 + 2 zeta w0 s + w0^2), w0 = 2 pi f0,
% applied to the model force u sampled uniformly with step dt (zero-order hold,
% exact discretization). Utility for FitFirstStretch.m: the ~1.6 kHz ringing
% after a ~2.5 ms stretch is phase-locked across repeats and absent in the
% length trace, i.e. measurement dynamics that also inflate the recorded peak.
w0 = 2*pi*f0;
A = [0 1; -w0^2 -2*zeta*w0];
B = [0; w0^2];
Ad = expm(A*dt);
Bd = A\((Ad - eye(2))*B);
x = [u(1); 0];                 % start at rest in equilibrium with u(1)
y = zeros(size(u));
for k = 1:numel(u)
    y(k) = x(1);
    x = Ad*x + Bd*u(k);
end
end
