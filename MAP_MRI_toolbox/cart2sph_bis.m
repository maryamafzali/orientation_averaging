function angles = cart2sph_bis(angles0, rev)

if ~exist('rev', 'var') || ~rev
    [TH, PHI, R] = cart2sph(angles0(:,1), angles0(:,2), angles0(:,3));

    angles(:,1) = pi/2-PHI;
    angles(:,2) = TH;
    angles(:,3) = R;
else
    angles = (angles0(:,3)*[1 1 1]).*[sin(angles0(:,1)).*cos(angles0(:,2)), sin(angles0(:,1)).*sin(angles0(:,2)), cos(angles0(:,1))];
end
