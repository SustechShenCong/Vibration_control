function p0 = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap)
% GET_INITIAL_PO This function returns initial p0 for given z0.
%
% p0 = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap)
%
% DS: Dynamical system
% masterModes: master modes
% proj: linear or nonlinear
% autData: mapping data
% Wmap: autonomous part of SSM


% WW = W(masterModes,:);
% pl = WW*z0; % linear projection
U  = DS.spectrum.W(:,masterModes);
pl = U'*DS.B*z0; % linear projection
disp("pl")
disp(pl)
switch proj
    case 'linear'
        p0 = pl; 
    case 'nonlinear'
        % mapping complex pl to real form
        % setup mapping data
        dim     = numel(masterModes);
        lambda  = autData.lamd;
        mapData = map_spectrum(lambda);
        compx   = mapData.compx;
        realx   = mapData.realx;
        idxReal = mapData.idxReal;
        idxComp = mapData.idxComp;
        pl(compx(2:2:end),:) = [];
        % map initial guess to real form
        pv                   = zeros(dim, 1);
        pv(realx,:)          = real(pl(idxReal));
        pv(compx(1:2:end-1)) = real(pl(idxComp));
        pv(compx(2:2:end))   = imag(pl(idxComp));
        % call optimization routine to return p0 (complex form)
        [p0, fval]                   = proj2SSM(z0,'nonlinear',Wmap,mapData,pv);
end
disp("p0")
disp(p0)
end