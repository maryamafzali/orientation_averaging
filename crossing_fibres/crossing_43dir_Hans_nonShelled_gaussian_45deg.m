close all
clear all
clc

SNR = [5, 10, 25, 100, 500];
sigma_g = 1./SNR/sqrt(2);

noise_gen = zeros(9,123,2,5,3,100);
[s1, s2, s3, s4, s5, s6] = size(noise_gen);

% m = 0;
% for i = 1:s1 % b-values
%     for j = 1:s2 % gradient directions (61+23+19)
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
% 
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
fdir1 = [0.4,0.6,-sqrt(1-0.16-0.36)];
te = acos(-sqrt(1-0.16-0.36));
ph = acos(0.4/sin(te));
fdir2 = [sin(te+pi/4)*cos(ph) sin(te+pi/4)*sin(ph) cos(te+pi/4)];

D_par = 1;
D_per = 0.14;

L_max = 6;
kappa = [1,9]; % 1, 9
Nmax = [6,8];
% M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);

m = 10;

bvec = [zeros(3,ndir) bvec];
b = [zeros(ndir,1); b];
S = zeros(ndir,s1,4,3);
for i = 1:s1 % b-values
%     b = (i-1)*1.5*ones(length(bvec),1);
    for l = 1:s4 % sigma_g
        for n = 1:s5 % kappa = \infty, 1, 9
            if n == 1
                S(:,i,l,n) = 0.5*exp(-b(ndir*(i-1)+1:ndir*i)*D_per-b(ndir*(i-1)+1:ndir*i)...
                    *(D_par - D_per).*(bvec(:,ndir*(i-1)+1:ndir*i)'*fdir1').^2) + ...
                    0.5*exp(-b(ndir*(i-1)+1:ndir*i)*D_per-b(ndir*(i-1)+1:ndir*i)...
                    *(D_par - D_per).*(bvec(:,ndir*(i-1)+1:ndir*i)'*fdir2').^2); 
            elseif n == 2 || n == 3
                S(:,i,l,n) = 0.5*slow_exchange_b_tensor_SynthMeasWatsonHinderedDiffusion_PGSE([D_par D_per kappa(n-1)], ...
                    bvec(:,ndir*(i-1)+1:ndir*i)', b(ndir*(i-1)+1:ndir*i), fdir1',1) + ...
                    0.5*slow_exchange_b_tensor_SynthMeasWatsonHinderedDiffusion_PGSE([D_par D_per kappa(n-1)], ...
                    bvec(:,ndir*(i-1)+1:ndir*i)', b(ndir*(i-1)+1:ndir*i), fdir2',1);
            end
        end
    end
end

Sn = zeros(ndir,s1,s4,s5,s6);

for i = 1:s1 % b-values
    for l = 1:s4 % sigma_g
        for n = 1:3 % kappa
            for p = 1:100 % noise realization
%                 Sn(:,i,l,n,p) = sqrt((squeeze(S(:,i,l,n)) + squeeze(noise_gen(i,61+43+1:end,1,l,n,p))').^2 + ...
%                     squeeze(noise_gen(i,61+43+1:end,2,l,n,p))'.^2);
                Sn(:,i,l,n,p) = squeeze(S(:,i,l,n)) + squeeze(noise_gen(i,61+1:61+43,1,l,n,p))';
%                 Sn(:,i,l,n,p) = squeeze(S(:,i,l,n)) + squeeze(noise_gen(i,61+43+1:end,1,l,n,p))';
            end
        end
    end
end

S = Sn;

b = b'*1000;
bvecs = bvec';

big_delta = 43.1*10^-3;
little_delta = 10.6*10^-3;

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
