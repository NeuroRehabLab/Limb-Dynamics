% Function calculates velocity
% & acceleration from position by ftting a cubic spline
% Inputs: position (m, deg, rad), sampling frequency (Hz)

function [splkin] = splineKin(kin, posFs)

[nrow,ncol] = size(kin);
for icol = 1:ncol
   postime = [1/posFs:1/posFs:nrow/posFs]';
   pp1 = spline(postime,kin(:,icol));
   coefs = [pp1.coefs(:,1)*3,pp1.coefs(:,2)*2,pp1.coefs(:,3)];
   pp2 = mkpp(pp1.breaks,coefs);
   coefs = [pp2.coefs(:,1)*2,pp2.coefs(:,2)];
   pp3 = mkpp(pp2.breaks,coefs);
   splkin(:,1,icol) = ppval(pp1,postime);
   splkin(:,2,icol) = ppval(pp2,postime);
   splkin(:,3,icol) = ppval(pp3,postime);
end

