function [ MX, MY ] = getTransportMatrix( driftL, magL, magG, ene )

gamma = ene / 0.511;

k = sqrt( 586.0806 * magG / sqrt( gamma^2 - 1 ) ); % e/(mc)=586.0806 (SI unit)
magPhi = k .* magL;

md = cell(4, 1);
mx = cell(3, 1);
my = cell(3, 1);

for ii = 1:4
    md{ii} = [1 driftL(ii); 0 1];
end

mx{1} = [cos(magPhi(1)) sin(magPhi(1))/k(1); -k(1)*sin(magPhi(1)) cos(magPhi(1))];
mx{2} = [cosh(magPhi(2)) sinh(magPhi(2))/k(2); k(2)*sinh(magPhi(2)) cosh(magPhi(2))];
mx{3} = [cos(magPhi(3)) sin(magPhi(3))/k(3); -k(3)*sin(magPhi(3)) cos(magPhi(3))];

my{1} = [cosh(magPhi(1)) sinh(magPhi(1))/k(1); k(1)*sinh(magPhi(1)) cosh(magPhi(1))];
my{2} = [cos(magPhi(2)) sin(magPhi(2))/k(2); -k(2)*sin(magPhi(2)) cos(magPhi(2))];
my{3} = [cosh(magPhi(3)) sinh(magPhi(3))/k(3); k(3)*sinh(magPhi(3)) cosh(magPhi(3))];

MX = md{1}; MY = md{1};
for ii = 1:3
    MX = md{ii+1} * mx{ii} * MX;
    MY = md{ii+1} * my{ii} * MY;
end

end

