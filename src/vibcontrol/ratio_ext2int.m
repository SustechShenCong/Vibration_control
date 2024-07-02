function ratio = ratio_ext2int(DS,t,z,u)


%% right side
nt  = numel(t);
fnl = zeros(size(DS.A,1),nt);
for k=1:nt
    fnl(:,k) = DS.evaluate_Fnl(z(:,k));
end
Ext = DS.D*u;
if ~isempty(DS.E); Ext = Ext+DS.E(t); end
Int = DS.A*z+fnl;
Ext = sum(sqrt(Ext.^2));
Int = sum(sqrt(Int.^2)); 

ratio = DS.epsilon*Ext./Int;


%% residual
figure;
plot(t,ratio);
xlabel('t'); ylabel('normalized external force');

end
