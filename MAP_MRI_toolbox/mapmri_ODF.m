function I = mapmri_ODF(ind,A,beta,mag_rho,a)

C = zeros(size(ind,1),1);
for m = 1:size(ind,1)
%     if ind(m,3)/2 == round(ind(m,3)/2) && ind(m,2)/2 == round(ind(m,2)/2) && ind(m,1)/2 == round(ind(m,1)/2)
        coeff = 0;
        for i = 0:2:ind(m,1)
            for j = 0:2:ind(m,2)
                for k = 0:2:ind(m,3)
                    coeff = coeff + (-1)^((i+j+k)/2)*gamma((3+ind(m,1)+ind(m,2)+ind(m,3)-i-j-k)/2)*...
                        beta(1)^(ind(m,1)-i)*beta(2)^(ind(m,2)-j)*beta(3)^(ind(m,3)-k)/...
                        (factorial(ind(m,1)-i)*factorial(ind(m,2)-j)*factorial(ind(m,3)-k)*factd(i)*factd(j)*factd(k));
                end
            end
        end
        C(m) = sqrt(factorial(ind(m,1))*factorial(ind(m,2))*factorial(ind(m,3)))*coeff;
%     end
end
I = mag_rho^(-3)*C'*a/sqrt(4*pi^3*A(1)*A(2)*A(3));

end