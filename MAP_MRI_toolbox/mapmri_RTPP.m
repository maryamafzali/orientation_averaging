function RTPP = mapmri_RTPP(K_par,ux,a)

RTPP = K_par'*a/sqrt(2*pi*ux^2);

end