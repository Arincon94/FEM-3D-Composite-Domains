clc, clearvars, close all;

% Global dimensions 
Lx = 52e-3; Ly = 12e-3; Lz = 52e-3;             % Cell dimensions
nelx = 22; nely = 22; nelz = 20;                % Cell elements
ndof = 3*(nelx + 1)*(nely + 1)*(nelz + 1);      % Degrees of freedom
nele = nelx*nely*nelz;                          % Number of elements

% Material properties
load('DM_Mid.mat');
omega = 2*pi*(0:5:2000);
E = 1e9*fE_2(omega/2/pi); eta = fL_2(omega/2/pi);   % Resonator properties in frequency domain
wt = 1000*2*pi; Fi = find(wt == omega);             % Frequency of interest
nu = 0.3;                                           % Poisson's ratio
par.E(1) = 2.35e9; par.E(2) = E(Fi);                % Stiffness: (1) walls, (2) resonator
par.eta(1) = 0.01; par.eta(2) = eta(Fi);            % Loss factor
par.rho(1) = 1160; par.rho(2) = 1130;               % Density
clear E eta fE_Int fL_Int

% Void material properties
xmin = 1e-11;                                       
par.E(3) = xmin*par.E(2);
par.eta(3) = xmin*par.eta(2);
par.rho(3) = xmin*par.rho(2);

% Top material properties
E1 = 144.8e9; E2 = 9.65e9; E3 = E2;
G12 = 4.14e9; G13 = G12; G23 = 3.45e9;
nu12 = 0.3; nu13 = 0.3; nu23 = 0.4;
par.eta(4) = 0.01;
par.rho(4) = 1389;
D = orthocl(E1,E2,E3,G23,G13,G12,nu23,nu13,nu12);

% Mesh definition
node_1 = [];

for k = 1:nely
    for j = 1:nelx
        for i = 1:nelz
            node_1 = [node_1; (j-1)*(nelz+1) + k*(nelz+1)*(nelx+1) + (i + 1)];
        end
    end
end

node_2 = node_1 + (nelz + 1);
node_3 = node_2 - (nelx + 1)*(nelz + 1);
node_4 = node_1 - (nelx + 1)*(nelz + 1);
node_5 = node_1 - 1;
node_6 = node_2 - 1;
node_7 = node_3 - 1;
node_8 = node_4 - 1;

edofMat = [3*node_1-2, 3*node_1-1, 3*node_1,...
    3*node_2-2, 3*node_2-1, 3*node_2,...
    3*node_3-2, 3*node_3-1, 3*node_3,...
    3*node_4-2, 3*node_4-1, 3*node_4,...
    3*node_5-2, 3*node_5-1, 3*node_5,...
    3*node_6-2, 3*node_6-1, 3*node_6,...
    3*node_7-2, 3*node_7-1, 3*node_7,...
    3*node_8-2, 3*node_8-1, 3*node_8];

% Element dimensions 
l1 = 1e-3*[1, 1, 1];            
l2 = 1e-3*[2.5, 1, 1];          
l3 = 1e-3*[1, 2.5, 1];
l4 = 1e-3*[2.5, 2.5, 1];

% Element regions 
[r1_i, r1_j, r1_k] = meshgrid(0, 0:nelz-1, 0);
r11_id = r1_k*(nelx)*(nelz) + r1_i*(nelz) + r1_j + 1;
[r1_i, r1_j, r1_k] = meshgrid(nelx-1, 0:nelz-1, 0);
r12_id = r1_k*(nelx)*(nelz) + r1_i*(nelz) + r1_j + 1;
[r1_i, r1_j, r1_k] = meshgrid(0, 0:nelz-1, nely-1);
r13_id = r1_k*(nelx)*(nelz) + r1_i*(nelz) + r1_j + 1;
[r1_i, r1_j, r1_k] = meshgrid(nelx-1, 0:nelz-1, nely-1);
r14_id = r1_k*(nelx)*(nelz) + r1_i*(nelz) + r1_j + 1;
first_region = [r11_id(:); r12_id(:); r13_id(:); r14_id(:)];
first_region = sort(first_region);

[r2_i, r2_j, r2_k] = meshgrid(1:nelx-2, 0:nelz-1, 0);
r21_id = r2_k*(nelx)*(nelz) + r2_i*(nelz) + r2_j + 1;
[r2_i, r2_j, r2_k] = meshgrid(1:nelx-2, 0:nelz-1, nely-1);
r22_id = r2_k*(nelx)*(nelz) + r2_i*(nelz) + r2_j + 1;
second_region = [r21_id(:); r22_id(:)];
second_region = sort(second_region);

[r3_i, r3_j, r3_k] = meshgrid(0, 0:nelz-1, 1:nely-2);
r31_id = r3_k*(nelx)*(nelz) + r3_i*(nelz) + r3_j + 1;
[r3_i, r3_j, r3_k] = meshgrid(nelx-1, 0:nelz-1, 1:nely-2);
r33_id = r3_k*(nelx)*(nelz) + r3_i*(nelz) + r3_j + 1;
third_region = [r31_id(:); r33_id(:)];
third_region = sort(third_region);

[r4_i, r4_j, r4_k] = meshgrid(1:nelx-2, 0:nelz-1, 1:nely-2);
r4_id = r4_k*(nelx)*(nelz) + r4_i*(nelz) + r4_j + 1;
fourth_region = r4_id(:);
fourth_region = sort(fourth_region);

% Region indices
[st_i, st_j, st_k] = meshgrid(0:nelx-1, 0, 0:nely-1);
st_id = st_k*(nelx)*(nelz) + st_i*(nelz) + st_j + 1;
[sb_i, sb_j, sb_k] = meshgrid(0:nelx-1, nelz-1, 0:nely-1);
sb_id = sb_k*(nelx)*(nelz) + sb_i*(nelz) + sb_j + 1;
stiffness_indices = [st_id(:); sb_id(:)];

% Resonator's topology definition
load('Resonators_topologies.mat');      
% Load the matrix X, which contains 36 different resonator topologies.
% Each column (with 4000 elements) represents a distinct topology,
% using 0 and 1 to indicate void and solid elements, respectively.

Xresonator = reshape(X(:,23),[10,20,20]);
clear X
xresonator = zeros(nelz,nelx,nely);

for i = 2:nely - 1
    B = xresonator(:,:,i);
    B(6:end-5,2:end-1) = Xresonator(:,:,i-1);
    xresonator(:,:,i) = B;
end

resonator_indices = find(xresonator == 1);
void_indices = setdiff(fourth_region, resonator_indices);

% Global matrices assembly
K = sparse(ndof, ndof); M = sparse(ndof, ndof);

KEC1 = ElementStiffness_comp(D, l1(1), l1(2), l2(3));
KEC2 = ElementStiffness_comp(D, l2(1), l2(2), l2(3));
KEC3 = ElementStiffness_comp(D, l3(1), l3(2), l3(3));
KEC4 = ElementStiffness_comp(D, l4(1), l4(2), l4(3));

tic;
for el = 1:nele
    edofs = edofMat(el,:);

    if ismember(el, first_region)
        KE = ElementStiffness(1, nu, l1(1), l1(2), l1(3));
        ME = ElementMass(1, l1(1), l1(2), l1(3));

        if ismember(el, stiffness_indices)
            K(edofs, edofs) = K(edofs, edofs) + (1 + 1j*par.eta(4))*KEC1;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(4)*ME;
        else
            K(edofs, edofs) = K(edofs, edofs) + par.E(1)*(1 + 1j*par.eta(1))*KE;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(1)*ME;
        end

    elseif ismember(el, second_region)
        KE = ElementStiffness(1, nu, l2(1), l2(2), l2(3));
        ME = ElementMass(1, l2(1), l2(2), l2(3));

        if ismember(el, stiffness_indices)
            K(edofs, edofs) = K(edofs, edofs) + (1 + 1j*par.eta(4))*KEC2;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(4)*ME;
        else
            K(edofs, edofs) = K(edofs, edofs) + par.E(1)*(1 + 1j*par.eta(1))*KE;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(1)*ME;
        end

    elseif ismember(el, third_region)
        KE = ElementStiffness(1, nu, l3(1), l3(2), l3(3));
        ME = ElementMass(1, l3(1), l3(2), l3(3));

        if ismember(el, stiffness_indices)
            K(edofs, edofs) = K(edofs, edofs) + (1 + 1j*par.eta(4))*KEC3;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(4)*ME;
        else
            K(edofs, edofs) = K(edofs, edofs) + par.E(1)*(1 + 1j*par.eta(1))*KE;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(1)*ME;
        end

    else 
        KE = ElementStiffness(1, nu, l4(1), l4(2), l4(3));
        ME = ElementMass(1, l4(1), l4(2), l4(3));

        if ismember(el, stiffness_indices)
            K(edofs, edofs) = K(edofs, edofs) + (1 + 1j*par.eta(4))*KEC4;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(4)*ME;
        elseif ismember(el, resonator_indices)
            K(edofs, edofs) = K(edofs, edofs) + par.E(2)*(1 + 1j*par.eta(2))*KE;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(2)*ME;
        else 
            K(edofs, edofs) = K(edofs, edofs) + par.E(3)*(1 + 1j*par.eta(3))*KE;
            M(edofs, edofs) = M(edofs, edofs) + par.rho(3)*ME;
        end

    end

    percent = (el/nele)*100;
    disp(['Percent: ' sprintf('%4.4f%%',percent)]);

end
toc;

%% Ploting the unit cell

axe_x = 1e-3*[0 1 3.5 6 8.5 11 13.5 16 18.5 21 23.5 26 28.5 31 33.5 36 38.5 41 43.5 46 48.5 51 52];
axe_y = sort(axe_x, 'descend');
axe_z = 1e-3*[20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0];

ui = [];

for j = 1:length(axe_y)
    for i = 1:length(axe_x)
        for k = 1:length(axe_z)
            ui = [ui; axe_x(i); axe_z(k); axe_y(j)];
        end
    end
end

figure;
corners = setdiff(first_region, stiffness_indices);
walls_l = setdiff(second_region, stiffness_indices);
walls_2 = setdiff(third_region, stiffness_indices); 

face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];

ui = ui(edofMat);
UF = ui;

for t = 1:nelx*nely*nelz
    Element = UF(t,:);
    vert = [Element(1:3); Element(4:6); Element(7:9); Element(10:12); 
        Element(13:15); Element(16:18); Element(19:21); Element(22:24)];
    vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
    if ismember(t,corners)
        patch('Faces',face,'Vertices',vert,'FaceColor','[0.9290 0.6940 0.1250]','FaceAlpha',0.1,'EdgeColor','none');
        hold on
    elseif ismember(t,walls_l)
        patch('Faces',face,'Vertices',vert,'FaceColor','[0.9290 0.6940 0.1250]','FaceAlpha',0.1,'EdgeColor','none');
        hold on
    elseif ismember(t,walls_2)
        patch('Faces',face,'Vertices',vert,'FaceColor','[0.9290 0.6940 0.1250]','FaceAlpha',0.1,'EdgeColor','none');
        hold on
    elseif ismember(t,stiffness_indices)
        patch('Faces',face,'Vertices',vert,'FaceColor','magenta','FaceAlpha',0.2,'EdgeColor','none');
        hold on
    elseif ismember(t,resonator_indices)
        patch('Faces',face,'Vertices',vert,'FaceColor','green','EdgeColor','#4DBEEE','EdgeAlpha','0.7');
        hold on
    else
    end
    clear vert
end
axis equal; axis tight; box off; 

xticks([0 0.025 0.05])
yticks([-0.05 -0.025 0])
zticks([0 0.01 0.02])
xlabel('x')
ylabel('y')
zlabel('z')

at = gca;
set(at,'FontName','Times New Roman','FontSize',12)

view(30,10)

%% FRF - Strain energy

% For the frequency response, it is assumed imposed displacement
% in the external nodes of all walls.

% Imposed displacement - Nodes definition
[ii, ij, ik] = meshgrid(0,0:nelz,0:nely); % Left 
left_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
left_dofs = [3*left_id(:) - 2; 3*left_id(:) - 1; 3*left_id(:)];
clear ii ij ik

[ii, ij, ik] = meshgrid(nelx,0:nelz,0:nely); % Right
right_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
right_dofs = [3*right_id(:) - 2; 3*right_id(:) - 1; 3*right_id(:)];
clear ii ij ik

[ii, ij, ik] = meshgrid(0:nelx,0:nelz,0); % Back
back_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
back_dofs = [3*back_id(:) - 2; 3*back_id(:) - 1; 3*back_id(:)];
clear ii ij ik

[ii, ij, ik] = meshgrid(0:nelx,0:nelz,nely); % Front
front_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
front_dofs = [3*front_id(:) - 2; 3*front_id(:) - 1; 3*front_id(:)];
clear ii ij ik

bound = union(union(left_dofs,right_dofs),union(back_dofs,front_dofs));
alldofs = 1:ndof;
internal = setdiff(alldofs, bound);

% Imposed displacement - Nodal displacements vector definition
U = zeros(ndof,1);
Ub = 5e-3;

% Face I
for i = 0:nely
    clear ii ij ik nodes nodes_id
    [ii, ij, ik] = meshgrid(0,0:nelz,i);
    nodes_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
    nodes = 3*nodes_id(:);
    U(nodes) = (-0.25*Ub/nely)*i + Ub;
end

% Face II
for i = 0:nelx
    clear ii ij ik nodes nodes_id
    [ii, ij, ik] = meshgrid(i,0:nelz,0);
    nodes_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
    nodes = 3*nodes_id(:);
    U(nodes) = (-0.25*Ub/nelx)*i + Ub;
end

% Face III
for i = 0:nely
    clear ii ij ik nodes nodes_id
    [ii, ij, ik] = meshgrid(nelx,0:nelz,i);
    nodes_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
    nodes = 3*nodes_id(:);
    U(nodes) = (-0.25*Ub/nely)*i + 0.75*Ub;
end

% Face IV
for i = 0:nelx
    clear ii ij ik nodes nodes_id
    [ii, ij, ik] = meshgrid(i,0:nelz,nely);
    nodes_id = ik*(nelx+1)*(nelz+1) + ii*(nelz+1)+(nelz+1-ij);
    nodes = 3*nodes_id(:);
    U(nodes) = (-0.25*Ub/nelx)*i + 0.75*Ub;
end

K = (K + K.')/2;
M = (M + M')/2;

Mii = M; Mii(:,bound) = []; Mii(bound,:) = [];
Mib = M; Mib(:,internal) = []; Mib(bound,:) = [];
Kii = K; Kii(:,bound) = []; Kii(bound,:) = [];
Kib = K; Kib(:,internal) = []; Kib(bound,:) = [];

Energy = zeros(length(omega),1);

for i = 1:length(omega)
    Dib = (-omega(i)^2*Mib + Kib);
    Dii = (-omega(i)^2*Mii + Kii);
    U(internal) = -Dii\(Dib*U(bound));
    Energy(i) = 0.5*abs(U'*K*U);
    percent = (i/length(omega))*100;
    disp(['Percent: ' sprintf('%4.4f%%',percent)]);
end

figure;
semilogy(omega/2/pi, Energy,'k')
legend('show','Location','best')
xlabel('Frequency (Hz)')
ylabel(' Strain energy (N/m)')
grid on
