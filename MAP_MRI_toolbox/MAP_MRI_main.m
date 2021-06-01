close all
clear all
clc

big_delta = 43.1*10^-3;
little_delta = 10.6*10^-3;

load bvals2 % bvals
load bvecs2 % bvecs

load mask72.mat % mask
mask = A; 

load yd72.mat % image file
Simg = yd1;

S = zeros(size(Simg,4),size(Simg,1),size(Simg,2),size(Simg,3));

for i = 1:size(Simg,4)
    S(i,:,:,:) = Simg(:,:,:,i);
end

Nmax = 6;
M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);

S0 = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
x = zeros(50,size(Simg,1),size(Simg,2),size(Simg,3));
RTOP = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
RTAP = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
RTPP = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
NG = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
NG_per = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
NG_par = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
PA_DTI = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
PA = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
delta_theta = zeros(size(Simg,1),size(Simg,2),size(Simg,3));
n = 10;
I = zeros(size(Simg,1),size(Simg,2),size(Simg,3),(n+1)*(n+1),3);

% a_coeff = zeros(size(Simg,1),size(Simg,2),size(Simg,3),50);
% V = zeros(3,3,size(Simg,1),size(Simg,2));
% D = zeros(3,size(Simg,1),size(Simg,2));
% FA = zeros(size(Simg,1),size(Simg,2));

method = 'MAPL'; % It can be 'MAP' or 'MAPL'

td = big_delta - little_delta/3;
  
S(S == 0) = 1e-7;

for i = 1:size(Simg,1)
    for j = 1:size(Simg,2)
        for k = 1:size(Simg,3)
            if mask(i,j,k)~=0
                [RTOP(i,j,k),RTAP(i,j,k),RTPP(i,j,k),NG(i,j,k),...
                    NG_par(i,j,k),NG_per(i,j,k),PA_DTI(i,j,k), PA(i,j,k), ...
                    delta_theta(i,j,k),I(i,j,k,:,:),a_coeff,O_coeff,phi,phi_iso,S0,D,V,FA] = mapmri(S(:,i,j,k),bvals2,bvecs2, ...
                    big_delta,little_delta,Nmax,method,n);
                i
                j
                k
            end
        end
    end
end

% for i = 1:size(Simg,1)
%     for j = 1:size(Simg,2)
%         for k = 1:size(Simg,3)
%             if mask(i,j,k)~=0
%                 [RTOP(i,j,k),RTAP(i,j,k),RTPP(i,j,k),NG(i,j,k),...
%                     NG_par(i,j,k),NG_per(i,j,k),PA_DTI(i,j,k), PA(i,j,k), ...
%                     delta_theta(i,j,k),I(i,j,k,:,:)] = mapmri_stiffness(S(:,i,j,k),bvals2,bvecs2, ...
%                     big_delta,little_delta,Nmax,method,n);
%                 i
%                 j
%                 k
%             end
%         end
%     end
% end


figure,colormap(gray),imagesc(rot90(flipud(PA(:,:,72)))), axis off, axis equal, colorbar, title('PA') 
figure,colormap(gray),imagesc(rot90(flipud(PA_DTI(:,:,72)))), axis off, axis equal, colorbar, title('PA_{DTI}') 
figure,colormap(gray),imagesc(rot90(flipud(delta_theta(:,:,72)))), axis off
axis equal, colorbar, title('\Delta\theta_{PO}'),caxis([-0.05,0.15]) 

[X,Y,Z] = sphere(n);
% xyz = [X(:),Y(:),Z(:)];
sx = size(X);

figure
for i = 1:(size(Simg,1))
    for j = 1:(size(Simg,2))
        for k = 72%1:size(Simg,3)
            if mask(i,j,k)~=0
            xyze(:,:) = I(i,j,k,:,:);
            m_xyze = sqrt(xyze(:,1).^2 + xyze(:,2).^2 + xyze(:,3).^2);
%             max_xyze = max(m_xyze);
%             xyze = [xyze(:,1)./max_xyze xyze(:,2)./max_xyze xyze(:,3)./max_xyze];
            Xe = reshape(xyze(:,1),sx);
            Ye = reshape(xyze(:,2),sx);
            Ze = reshape(xyze(:,3),sx);
            surf(Xe,Ye,Ze,reshape(m_xyze/max_xyze,sx(1),sx(1)));
%             surf(Xe+2*i,Ye+2*j,Ze,reshape(m_xyze/max_xyze,sx(1),sx(1)));
            hold on
            end
        end
    end
end

axis equal
axis vis3d
camlight right
lighting phong
shading interp
xlabel X
ylabel Y
zlabel Z
% view(0,90)
title('diffusion-tensor')

