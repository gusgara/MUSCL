function [sl] = VanAlbada(r)
phi = (r.^2 +r)./(r.^2+1);
 sl = max(0,phi);
end