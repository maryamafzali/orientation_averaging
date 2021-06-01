function NG_par = mapmri_NG_par(a0_ver,a_ver,a)

NG_par = sqrt(1-norm(a0_ver.*a)^2/norm(a_ver.*a)^2);

end       