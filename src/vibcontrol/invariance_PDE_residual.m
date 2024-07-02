% function res = invariance_PDE_residual(DS,t,z,u)
function [res, res_A, res_R] = invariance_PDE_residual(DS,t,z,u)
% Absolute error and relative error
%% left side
ts = t(2:end-1);
zs = z(:,2:end-1);
us = u(:,2:end-1);
dt = t(2)-t(1);
% dot(z)
zdot = (z(:,3:end)-z(:,1:end-2))/(2*dt); % center difference
lhs  = DS.B*zdot;

%% right side
nt  = numel(ts);
fnl = zeros(size(DS.A,1),nt);
for k=1:nt
    fnl(:,k) = DS.evaluate_Fnl(zs(:,k));
end
Dus = DS.D*us;
if ~isempty(DS.E); Dus = Dus+DS.E(ts); end
rhs = DS.A*zs+fnl+DS.epsilon*[Dus; zeros(size(Dus))];

%% residual
res1 = sqrt(sum((lhs-rhs).^2)); 
res_A = res1;
% ref = sqrt(sum(Dus.^2)); ref = max(ref);
% res = res/(DS.epsilon*ref);
% res = res./max(sqrt(sum(zs.^2)));
res2 = res1./sqrt(sum((lhs).^2));
res_R = res2;

if norm(res1,'inf')<0.05
    fprintf('Plot absolute invariance error\n');
    res = res1;
else
    fprintf('Plot relative invariance error\n');
    res = res2;
end
figure;
% plot(ts,res);
plot(ts,res2);
xlabel('t'); ylabel('normalized invariance error');

end
