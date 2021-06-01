function [PA , theta_DTI] = mapmri_PA(O,T,A,u0,a)

cos_PQ = O'*T*a/(norm(a)*norm(O))*((A(1)*A(2)*A(3))/u0^6)^0.25;

sin_PQ = sqrt(1-cos_PQ^2);

theta_DTI = acos(cos_PQ);

PA = sin_PQ^(3*0.4)/(1-3*sin_PQ^0.4 +3*sin_PQ^0.8);

end







