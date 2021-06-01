function NG_per = mapmri_NG_per(p0_par,p_par,a)

NG_per = sqrt(1-norm(p0_par.*a)^2/norm(p_par.*a)^2);

end