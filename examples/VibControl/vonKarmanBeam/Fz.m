function u = Fz(t)
    load("uz.mat","uz")
    load("time.mat","time")
    u = interp1(time,uz,t,'linear','extrap') * 0.001;
end
