% Function calculates velocity & acceleration
% from position by fitting a cubic spline

% Inputs:
%   kin = position (m, deg, rad); posFs = sampling frequency (Hz)
% Outputs:
%   splkin(:,1) = position
%   splkin(:,2) = velocity
%   splkin(:,3) = acceleration

function [splkin] = splineKin_140916(kin, posFs)

[nrow,ncol] = size(kin);

for icol = 1:ncol
   postime = [1/posFs:1/posFs:nrow/posFs]';                    % Time
   pp1 = spline(postime,kin(:,icol));                          % Position function
   coefs = [pp1.coefs(:,1)*3,pp1.coefs(:,2)*2,pp1.coefs(:,3)]; % Coefficients of dy/dx (velocity)
   pp2 = mkpp(pp1.breaks,coefs);                               % Velocity function
   coefs = [pp2.coefs(:,1)*2,pp2.coefs(:,2)];                  % Coefficients of d2y/dx2 (acceleration)
   pp3 = mkpp(pp2.breaks,coefs);                               % Acceleration function
   splkin(:,1,icol) = ppval(pp1,postime);    % Setting position velocity and acceleration to output array
   splkin(:,2,icol) = ppval(pp2,postime);
   splkin(:,3,icol) = ppval(pp3,postime);
end

end

