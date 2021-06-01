function [PA_DTI , theta_DTI] = mapmri_PA_DTI(ux,uy,uz,u0)

cos_DTI = (8*u0^3*ux*uy*uz)/((ux^2 + u0^2)*(uy^2 + u0^2)*(uz^2 + u0^2));

sin_DTI = sqrt(1-cos_DTI);

theta_DTI = acos(sqrt(cos_DTI));

PA_DTI = sin_DTI^(3*0.4)/(1-3*sin_DTI^0.4 +3*sin_DTI^0.8);

end


