close all
clear all
clc

SNR = [5, 10, 25, 100, 500];
sigma_g = 1./SNR/sqrt(2);

noise_gen = zeros(9,123,2,5,3,100);
[s1, s2, s3, s4, s5, s6] = size(noise_gen);

% m = 0;
% for i = 1:s1 % b-values
%     for j = 1:s2 % gradient directions (61+43+19)
%         for k = 1:s3 % real and imaginary part of the noise
%             for l = 1:s4 % sigma_g
%                 for n = 1:s5 % kappa
%                     for p = 1:s6 % noise realization
%                         rng(m) % set the seed
%                         m = m+1;
%                         noise_gen(i,j,k,l,n,p) = normrnd(0,sigma_g(l));
%                     end
%                 end
%             end
%         end
%     end
% end

% save('noise_new.mat','noise_gen')
load noise_new.mat

% load X_061_Qdistr_03632__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat % 61 directions for one shell
% bvec = allvec(:,1:2:122);

load X_344_Qdistr_03753__1.00_2.00_2.00_1.00_0.00_J0_S1_4.00_1.50_0.23.mat % non_shelled 
% load X_488_Qdistr_03760__1.00_2.00_2.00_1.00_0.00_J0_S1_4.00_1.50_0.21.mat % shelled

bvec = rivec1(:,1:2:end);

ndir = 43;

for i = 1:length(bvec)
    norm_bvec(i) = (norm(bvec(:,i)).^2);
    bvec(:,i) = bvec(:,i)/sqrt(norm_bvec(i));
end

b = norm_bvec'/max(norm_bvec)*12;

% [leb_tmp] = getLebedevSphere(38);
%
% x = leb_tmp.x;
% y = leb_tmp.y;
% z = leb_tmp.z;
%
% vec = [x, y, z];
%
% w_leb = leb_tmp.w;
%
% bvec = vec';
% bvec = bvec';
% bvec1 = bvec;
% for i = 1:length(bvec)
%     for j = i:length(bvec)
%         if bvec1(i,:) == -bvec(j,:)
%             bvec1(i,:) = 0;
%         end
%     end
% end
%
% w_leb(all(~bvec1,2)) = [];
%
% bvec1( all(~bvec1,2), : ) = [];
%
% bvec = bvec1';

% gradDWI = bvec';
% bvec = bvec';
%
% ndir = length(bvec);
%
% anglesIn  = cart2sph_bis(gradDWI);
%
% theta = anglesIn(:,1);
% phi = anglesIn(:,2);
%
fdir = [0.4,0.6,-sqrt(1-0.16-0.36)];
D_par = 1;
D_per = 0.14;

L_max = 6;
kappa = [1,9]; % 1, 9
Nmax = [6,8];
% M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);

method = 'FREE';
m = 10;

% R = (L_max+1)*(L_max+2)/2;
%
% Y = zeros(ndir,R,L_max/2);
% for L = 2:2:L_max
%     for k=0:2:L
%         Pkm=legendre(k, cos(theta))';
%         for m=-k:k
%             j=k*(k+1)/2+m+1;
%             if m<0
%                 Y(:,j,L/2)=sqrt(((2*k+1)/(2*pi))*factorial(k+m)/factorial(k-m))*Pkm(:,-m+1).*cos(m*phi);
%             elseif m==0
%                 Y(:,j,L/2)=sqrt((2*k+1)/(4*pi))*Pkm(:,1).*ones(length(phi),1);
%             else
%                 Y(:,j,L/2)=(-1)^m*sqrt(((2*k+1)/(2*pi))*factorial(k-m)/factorial(k+m))*Pkm(:,m+1).*sin(m*phi);
%             end
%         end
%     end
% end
bvec = [zeros(3,ndir) bvec];
b = [zeros(ndir,1); b];
S = zeros(ndir,s1,4,3);
for i = 1:s1 % b-values
%     b = (i-1)*1.5*ones(ndir,1);
    for l = 1:s4 % sigma_g
        for n = 1:3 % kappa = \infty, 1, 9
            if n == 1
                S(:,i,l,n) = exp(-b(ndir*(i-1)+1:ndir*i)*D_per-b(ndir*(i-1)+1:ndir*i)...
                    *(D_par - D_per).*(bvec(:,ndir*(i-1)+1:ndir*i)'*fdir').^2);
            elseif n == 2 || n == 3
                S(:,i,l,n) = slow_exchange_b_tensor_SynthMeasWatsonHinderedDiffusion_PGSE([D_par D_per kappa(n-1)], ...
                    bvec(:,ndir*(i-1)+1:ndir*i)', b(ndir*(i-1)+1:ndir*i), fdir',1);
            end
        end
    end
end

Sn = zeros(ndir,s1,4,3,100);
for i = 1:s1 % b-values
    for k = 1:2 % real and imaginary part of the noise
        for l = 1:s4 % sigma_g
            for n = 1:3 % kappa
                for p = 1:100 % noise realization
                    Sn(:,i,l,n,p) = sqrt((squeeze(S(:,i,l,n)) + squeeze(noise_gen(i,61+1:61+43,1,l,n,p))').^2 + ...
                        squeeze(noise_gen(i,61+1:61+43,2,l,n,p))'.^2);
%                     Sn(:,i,l,n,p) = squeeze(S(:,i,l,n)) + squeeze(noise_gen(i,61+1:61+43,1,l,n,p))';
                end
            end
        end
    end
end

S = Sn;

% [Vwgts,~,~,~] = Optimize_SPH_filters(bvec');
% % these are the weights
% w_hans = Vwgts{1,1};
%
% g = bvec;
% G = [g(:,1).^2, g(:,2).^2, g(:,3).^2, 2*g(:,1).*g(:,2), 2*g(:,2).*g(:,3), 2*g(:,1).*g(:,3)];
%
% dir_ave_spharm = zeros(s1,L_max/2,s4,3,100);
% arithmatic_sum = zeros(s1,s4,3,100);
% S_mean_SH = zeros(s1,s4,3,100);
% mean_Hans = zeros(s1,s4,3,100);
% mean_Leb = zeros(s1,s4,3,100);
% for L = 2:2:L_max % L_max
%     R = (L+1)*(L+2)/2;
%     for n = 1:3 % kappa
%         for i = 1:s1 % b-value
%             for l = 1:s4 % sigma_g
%                 for p = 1:100 % noise realization
%                     a = inv(Y(:,1:R,L/2)'*Y(:,1:R,L/2))*Y(:,1:R,L/2)'*squeeze(S(:,i,l,n,p));
%                     dir_ave_spharm(i,L/2,l,n,p) = a(1)/sqrt(4*pi);
%                     arithmatic_sum(i,l,n,p) = mean(S(:,i,l,n,p));
%
%                     ungama = lscov(G,squeeze(S(:,i,l,n,p)));
%                     Dlls = [ungama(1),ungama(4),ungama(6);ungama(4),ungama(2),ungama(5);ungama(6),ungama(5),ungama(3)];
%                     S_mean_SH(i,l,n,p) = trace(Dlls)/3;
%
%                     mean_Hans(i,l,n,p) = sum(w_hans.*squeeze(S(:,i,l,n,p)))/sum(w_hans);
% %                     mean_Leb(i,l,n,p) = sum(w_leb.*squeeze(S(:,i,l,n,p)))/sum(w_leb);
%                 end
%             end
%         end
%     end
% end

% b = [0*ones(length(bvec),1); 1.5*ones(length(bvec),1); 3*ones(length(bvec),1);...
%     4.5*ones(length(bvec),1);  6*ones(length(bvec),1);...
%     7.5*ones(length(bvec),1);  9*ones(length(bvec),1); ...
%     10.5*ones(length(bvec),1);  12*ones(length(bvec),1)];

% bvecs = [bvec; bvec; bvec; bvec; bvec; bvec; bvec; bvec; bvec];

b = b'*1000;
bvecs = bvec';
% b = b*1000;

big_delta = 43.1*10^-3;
little_delta = 10.6*10^-3;

% [b_sort, order] = sort(b);
% [C, ia, ic] = unique(b_sort);
% pl = length(C);

GT(b~=0) = sqrt(pi)/2*exp(-D_per*b(b~=0)/1000).*erf(sqrt((D_par - D_per)*b(b~=0)/1000))...
    ./sqrt((D_par - D_per)*b(b~=0)/1000);

GT(b==0) = 1;

S_iso = zeros(length(b),3,2,3, s4, 100);

kappa_iso = zeros(95,3,2,3, s4, 100);

for p = 1:100 % noise realization
    for l = 1:s4 % sigma_g
        for j = 1:3 % kappa
            for i = 1:2 % Nmax
                for r = 1:3 % method
                    if r == 1
                        method = 'MAP';
                    elseif r == 2
                        method = 'MAPL';
                    elseif r == 3
                        method = 'FREE';
                    end
                    S1 = squeeze(Sn(:,:,l,j,p));
                    S1 = reshape(S1,[ndir*s1,1]);
                    M_sym = 1/6*(Nmax(i)/2+1)*(Nmax(i)/2+2)*(2*Nmax(i)+3);
                    M_sym = round(M_sym);
                    [kappa_iso(1:M_sym,j,i,r,l,p),E,u0(j,i,r,l,p)] = ...
                        mapmri_iso(S1,b,bvecs', big_delta,little_delta,Nmax(i),method);
                    
                    S_iso(:,j,i,r,l,p) = E*squeeze(kappa_iso(1:M_sym,j,i,r,l,p));
                    

                end
            end
        end
    end
    p
end



td = big_delta-little_delta/3;
% bvals1 = 0:1500:12000;
% bvecs1 = repmat([0,0,1],length(bvals1),1);
% Siso_shell = zeros(length(bvals1),3,2,3, s4, 100);
% q = bvecs1'.*repmat(sqrt(bvals1/td)/(2*pi),3,1);
% 
p = 1;
for N = 0:2:Nmax(1)
    for j = 1:Nmax(1)
        for l = 0:2:Nmax(1)
            for m = -l:l
                if l + 2*j - 2 == N
                    ind_sph6(p,:) = [j l m];
                    p = p + 1;
                end
            end
        end
    end
end

p = 1;
for N = 0:2:Nmax(2)
    for j = 1:Nmax(2)
        for l = 0:2:Nmax(2)
            for m = -l:l
                if l + 2*j - 2 == N
                    ind_sph8(p,:) = [j l m];
                    p = p + 1;
                end
            end
        end
    end
end

bvals2 = 0:1000:20000;
bvecs2 = repmat([0,0,1],length(bvals2),1);
Siso_new = zeros(length(bvals2),3,2,3, s4, 100);
q = bvecs2'.*repmat(sqrt(bvals2/td)/(2*pi),3,1);
 
clear E
for p = 1:100 % noise realization
    for l = 1:s4 % sigma_g
        for j = 1:3 % kappa
            for i = 1:2 % Nmax
                if i == 1;
                    ind_sph = ind_sph6;
                else
                    ind_sph = ind_sph8;
                end
                for r = 1:3 % method
                    if r == 1
                        method = 'MAP';
                    elseif r == 2
                        method = 'MAPL';
                    elseif r == 3
                        method = 'FREE';
                    end
                    M_sym = 1/6*(Nmax(i)/2+1)*(Nmax(i)/2+2)*(2*Nmax(i)+3);
                    M_sym = round(M_sym);
                    
                    for k = 1:size(ind_sph,1)
                        E(:,k) = mapmri_Xi(ind_sph(k,1),ind_sph(k,2),ind_sph(k,3),u0(j,i,r,l,p),q,bvecs2');
                    end
                    
                    Siso_new(:,j,i,r,l,p) = E(:,1:M_sym)*squeeze(kappa_iso(1:M_sym,j,i,r,l,p));
                    
                end
            end
        end
    end
    p
end

Siso_MAP6(:,:,:,:) = Siso_new(:,:,1,1,:,:);

bvals2 = 0:1500:12000;
bvecs2 = repmat([0,0,1],length(bvals2),1);
Siso_non_shell = zeros(length(bvals2),3,2,3, s4, 100);
q = bvecs2'.*repmat(sqrt(bvals2/td)/(2*pi),3,1);
 
clear E
for p = 1:100 % noise realization
    for l = 1:s4 % sigma_g
        for j = 1:3 % kappa
            for i = 1:2 % Nmax
                if i == 1;
                    ind_sph = ind_sph6;
                else
                    ind_sph = ind_sph8;
                end
                for r = 1:3 % method
                    if r == 1
                        method = 'MAP';
                    elseif r == 2
                        method = 'MAPL';
                    elseif r == 3
                        method = 'FREE';
                    end
                    M_sym = 1/6*(Nmax(i)/2+1)*(Nmax(i)/2+2)*(2*Nmax(i)+3);
                    M_sym = round(M_sym);
                    
                    for k = 1:size(ind_sph,1)
                        E(:,k) = mapmri_Xi(ind_sph(k,1),ind_sph(k,2),ind_sph(k,3),u0(j,i,r,l,p),q,bvecs2');
                    end
                    
                    Siso_non_shell(:,j,i,r,l,p) = E(:,1:M_sym)*squeeze(kappa_iso(1:M_sym,j,i,r,l,p));
                    
                end
            end
        end
    end
    p
end

% td = big_delta-little_delta/3;
% bvals1 = 0:1500:12000;
% bvecs1 = repmat([0,0,1],length(bvals1),1);
% Siso_shell = zeros(length(bvals1),3,2,3, s4, 100);
% q = bvecs1'.*repmat(sqrt(bvals1/td)/(2*pi),3,1);
% 
% % ind_sph = zeros(round(M_sym),3);
% p = 1;
% for N = 0:2:Nmax(1)
%     for j = 1:Nmax(1)
%         for l = 0:2:Nmax(1)
%             for m = -l:l
%                 if l + 2*j - 2 == N
%                     ind_sph6(p,:) = [j l m];
%                     p = p + 1;
%                 end
%             end
%         end
%     end
% end
% 
% p = 1;
% for N = 0:2:Nmax(2)
%     for j = 1:Nmax(2)
%         for l = 0:2:Nmax(2)
%             for m = -l:l
%                 if l + 2*j - 2 == N
%                     ind_sph8(p,:) = [j l m];
%                     p = p + 1;
%                 end
%             end
%         end
%     end
% end
% 
% clear E
% for p = 1:100 % noise realization
%     for l = 1:s4 % sigma_g
%         for j = 1:3 % kappa
%             for i = 1:2 % Nmax
%                 if i == 1;
%                     ind_sph = ind_sph6;
%                 else
%                     ind_sph = ind_sph8;
%                 end
%                 for r = 1:3 % method
%                     if r == 1
%                         method = 'MAP';
%                     elseif r == 2
%                         method = 'MAPL';
%                     elseif r == 3
%                         method = 'FREE';
%                     end
%                     S1 = squeeze(Sn(:,:,l,j,p));
%                     S1 = reshape(S1,[ndir*s1,1]);
%                     M_sym = 1/6*(Nmax(i)/2+1)*(Nmax(i)/2+2)*(2*Nmax(i)+3);
%                     M_sym = round(M_sym);
%                     u0(j,i,r,l,p) = mapmri_u0(S1,b,bvecs',td);
%                     
%                     for k = 1:size(ind_sph,1)
%                         E(:,k) = mapmri_Xi(ind_sph(k,1),ind_sph(k,2),ind_sph(k,3),u0(j,i,r,l,p),q,bvecs1');
%                     end
%                     
%                     Siso_shell(:,j,i,r,l,p) = E(:,1:M_sym)*squeeze(kappa_iso(1:M_sym,j,i,r,l,p));
%                     
%                 end
%             end
%         end
%     end
%     p
% end
% 
% Siso_non_shell = Siso_shell;
% 
% bvals2 = 0:1000:20000;
% bvecs2 = repmat([0,0,1],length(bvals2),1);
% Siso_new = zeros(length(bvals2),3,2,3, s4, 100);
% q = bvecs2'.*repmat(sqrt(bvals2/td)/(2*pi),3,1);
% 
% clear E
% for p = 1:100 % noise realization
%     for l = 1:s4 % sigma_g
%         for j = 1:3 % kappa
%             for i = 1:2 % Nmax
%                 if i == 1;
%                     ind_sph = ind_sph6;
%                 else
%                     ind_sph = ind_sph8;
%                 end
%                 for r = 1:3 % method
%                     if r == 1
%                         method = 'MAP';
%                     elseif r == 2
%                         method = 'MAPL';
%                     elseif r == 3
%                         method = 'FREE';
%                     end
% %                     S1 = squeeze(Sn(:,:,l,j,p));
% %                     S1 = reshape(S1,[ndir*s1,1]);
%                     M_sym = 1/6*(Nmax(i)/2+1)*(Nmax(i)/2+2)*(2*Nmax(i)+3);
%                     M_sym = round(M_sym);
% %                     u0(j,i,r,l,p) = mapmri_u0(S1,b,bvecs',td);
%                     
%                     for k = 1:size(ind_sph,1)
%                         E(:,k) = mapmri_Xi(ind_sph(k,1),ind_sph(k,2),ind_sph(k,3),u0(j,i,r,l,p),q,bvecs2');
%                     end
%                     
%                     Siso_new(:,j,i,r,l,p) = E(:,1:M_sym)*squeeze(kappa_iso(1:M_sym,j,i,r,l,p));
%                     
%                 end
%             end
%         end
%     end
%     p
% end
% 
% Siso_MAP6(:,:,:,:) = Siso_new(:,:,1,1,:,:);
% 
% d0_dir_ave_spharm = zeros(s1,L_max/2,s4,s6);
% d0_arithmatic_sum = zeros(s1,s4,s6);
% d0_S_mean_SH = zeros(s1,s4,s6);
% d0_mean_Hans = zeros(s1,s4,s6);
% d0_mean_Leb = zeros(s1,s4,s6);
% d0_Siso_shell = zeros(s1,2,3,s4,s6);
% 
% for n = 1:3 % kappa
%     for i = 1:s1 % b-value
%         for l = 1:s4 % sigma_g
%             for p = 1:100 % noise realization
%                 for L = 2:2:L_max % L_max
%                     R = (L+1)*(L+2)/2;
%                     d0_dir_ave_spharm(i,L/2,l,p) = 1/3*sum((dir_ave_spharm(i,L/2,l,:,p)-GT(i)));
%                 end
%                 d0_arithmatic_sum(i,l,p) = 1/3*sum((arithmatic_sum(i,l,:,p)-GT(i)));
%                 d0_S_mean_SH(i,l,p) = 1/3*sum((S_mean_SH(i,l,:,p) - GT(i)));
%                 d0_mean_Hans(i,l,p) = 1/3*sum((mean_Hans(i,l,:,p) - GT(i)));
%                 d0_mean_Leb(i,l,p) = 1/3*sum((mean_Leb(i,l,:,p) - GT(i)));
%                 for nm = 1:2 % Nmax
%                     for r = 1:3 % method
%                         d0_Siso_shell(i,nm,r,l,p) = 1/3*sum((Siso_shell(i,:,nm,r,l,p) - GT(i)));
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% d0 = zeros(s1,s4,s6,13);
% d0(:,:,:,4) = d0_dir_ave_spharm(:,1,:,:);
% d0(:,:,:,5) = d0_dir_ave_spharm(:,2,:,:);
% d0(:,:,:,6) = d0_dir_ave_spharm(:,3,:,:);
% d0(:,:,:,1) = d0_arithmatic_sum;
% d0(:,:,:,7) = d0_S_mean_SH;
% d0(:,:,:,3) = d0_mean_Hans;
% d0(:,:,:,2) = d0_mean_Leb;
% d0(:,:,:,8) = d0_Siso_shell(:,1,1,:,:);
% d0(:,:,:,9) = d0_Siso_shell(:,1,2,:,:);
% d0(:,:,:,10) = d0_Siso_shell(:,1,3,:,:);
% d0(:,:,:,11) = d0_Siso_shell(:,2,1,:,:);
% d0(:,:,:,12) = d0_Siso_shell(:,2,2,:,:);
% d0(:,:,:,13) = d0_Siso_shell(:,2,3,:,:);
% 
% d1_dir_ave_spharm = zeros(L_max/2,s4,s6);
% d1_arithmatic_sum = zeros(s4,s6);
% d1_S_mean_SH = zeros(s4,s6);
% d1_mean_Hans = zeros(s4,s6);
% d1_mean_Leb = zeros(s4,s6);
% d1_Siso_shell = zeros(2,3,s4,s6);
% 
% for l = 1:s4 % sigma_g
%     for p = 1:100 % noise realization
%         for L = 2:2:L_max % L_max
%             R = (L+1)*(L+2)/2;
%             d1_dir_ave_spharm(L/2,l,p) = 1/8*sum(abs(d0_dir_ave_spharm(2:end,L/2,l,p)));
%         end
%         d1_arithmatic_sum(l,p) = 1/8*sum(abs(d0_arithmatic_sum(2:end,l,p)));
%         d1_S_mean_SH(l,p) = 1/8*sum(abs(d0_S_mean_SH(2:end,l,p)));
%         d1_mean_Hans(l,p) = 1/8*sum(abs(d0_mean_Hans(2:end,l,p)));
%         d1_mean_Leb(l,p) = 1/8*sum(abs(d0_mean_Leb(2:end,l,p)));
%         for nm = 1:2 % Nmax
%             for r = 1:3 % method
%                 d1_Siso_shell(nm,r,l,p) = 1/8*sum(abs(d0_Siso_shell(2:end,nm,r,l,p)));
%             end
%         end
%     end
% end
% 
% d1 = zeros(s4,s6,13);
% d1(:,:,4) = d1_dir_ave_spharm(1,:,:);
% d1(:,:,5) = d1_dir_ave_spharm(2,:,:);
% d1(:,:,6) = d1_dir_ave_spharm(3,:,:);
% d1(:,:,1) = d1_arithmatic_sum;
% d1(:,:,7) = d1_S_mean_SH;
% d1(:,:,3) = d1_mean_Hans;
% d1(:,:,2) = d1_mean_Leb;
% d1(:,:,8) = d1_Siso_shell(1,1,:,:);
% d1(:,:,9) = d1_Siso_shell(1,2,:,:);
% d1(:,:,10) = d1_Siso_shell(1,3,:,:);
% d1(:,:,11) = d1_Siso_shell(2,1,:,:);
% d1(:,:,12) = d1_Siso_shell(2,2,:,:);
% d1(:,:,13) = d1_Siso_shell(2,3,:,:);
% 
% m_d1 = mean(d1,2);
% m_d1 = squeeze(m_d1);
% 
% std_d1 = zeros(s4,13);
% for i = 1:13 % b-value
%     for l = 1:s4 % sigma_g
%         std_d1(l,i) = std(squeeze(d1(l,:,i)));
%     end
% end
% 
% Rcc = zeros(2,2,s4,100,13);
% Pcc = zeros(2,2,s4,100,13);
% RL = zeros(2,2,s4,100,13);
% RU = zeros(2,2,s4,100,13);
% 
% m_Rcc = zeros(s4,13);
% std_Rcc = zeros(s4,13);
% for i = 1:13 % b-value
%     for l = 1:s4 % sigma_g
%         for p = 1:100
%             [Rcc(:,:,l,p,i),Pcc(:,:,l,p,i),RL(:,:,l,p,i),RU(:,:,l,p,i)] = ...
%                 corrcoef(C(2:end)',squeeze(d0(2:end,l,p,i)));
%         end
%         m_Rcc(l,i) = mean(Rcc(1,2,l,:,i),4);
%         std_Rcc(l,i) = std(squeeze(Rcc(1,2,l,:,i)));
%     end
% end
% 
% x = [1, 2, 3, 4, 4.25, 4.5, 5, 6, 6.25, 6.5, 7, 7.25, 7.5];
% 
% figure
% ha = tight_subplot(1,1,[.02 .06],[.26 .01],[.08 .01]);
% errorbar(x, m_Rcc(1,:), std_Rcc(1,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% hold on
% errorbar(x+0.03, m_Rcc(2,:), std_Rcc(2,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.06, m_Rcc(3,:), std_Rcc(3,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.09, m_Rcc(4,:), std_Rcc(4,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.09, m_Rcc(5,:), std_Rcc(5,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% plot(x,x*0,'k--','LineWidth',3)
% set(gca, 'FontSize', 36)
% set(gca,'xtick',[])
% ylabel('d_2')
% xlim([0.9,7.7])
% xtickangle(45)
% ylim([-1,1])
% xticks(x)
% xticklabels({'Arithmetic sum','Lebedev','Hans','SH (L = 2)','SH (L = 4)','SH (L = 6)','trace(M)/3',...
%     'MAP (N = 6)','MAPL (N = 6)','no constraint (N = 6)','MAP (N = 8)','MAPL (N = 8)',...
%     'no constraint (N = 8)'})
% legend('\sigma_g = 1/(5\surd{2})','\sigma_g = 1/(10\surd{2})','\sigma_g = 1/(25\surd{2})','\sigma_g = 1/(100\surd{2})',...
%     '\sigma_g = 1/(500\surd{2})')
% 
% figure
% ha = tight_subplot(1,1,[.02 .06],[.26 .01],[.08 .01]);
% errorbar(x, m_d1(1,:), std_d1(1,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% hold on
% errorbar(x+0.03, m_d1(2,:), std_d1(2,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.06, m_d1(3,:), std_d1(3,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.09, m_d1(4,:), std_d1(4,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% errorbar(x+0.09, m_d1(5,:), std_d1(5,:),'o','LineWidth',4,'MarkerSize',9,'CapSize',10)
% plot(x,x*0,'k--','LineWidth',3)
% set(gca, 'FontSize', 36)
% set(gca,'xtick',[])
% ylabel('d_1')
% xlim([0.9,7.7])
% xtickangle(45)
% ylim([0,0.125])
% xticks(x)
% xticklabels({'Arithmetic sum','Lebedev','Hans','SH (L = 2)','SH (L = 4)','SH (L = 6)','trace(M)/3',...
%     'MAP (N = 6)','MAPL (N = 6)','no constraint (N = 6)','MAP (N = 8)','MAPL (N = 8)',...
%     'no constraint (N = 8)'})
% legend('\sigma_g = 1/(5\surd{2})','\sigma_g = 1/(10\surd{2})','\sigma_g = 1/(25\surd{2})','\sigma_g = 1/(100\surd{2})',...
%     '\sigma_g = 1/(500\surd{2})')
