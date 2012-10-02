%--------------------------------------------------------------------------
% DRIVER
%--------------------------------------------------------------------------
function cylinders()

% Parameters
params.N = 2; % number of cylinders
params.a = 1; % radius of cylinders (all the same)
params.rhoc = 1; % density of cylinders (all the same)
params.rhof = 1; % density of fluid

% Set ODE solver options
ti = 0.0; % start time
tf = 10.0; % end time
odeopts = odeset(...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-4, ...
    'Mass', @(t, z) massMatrix(t, z, params), ...
    'MassSingular', 'no', ...
    'Stats', 'on');

% Set interpolation options
dt = 0.05;
tout = ti:dt:tf;

% Initialize cylinders and velocities
x = [0; 0]; % position
y = [-4; 4]; % position
u = [1; 1]; % velocity
v = [0; 0]; % velocity

z0 = zeros(4*params.N,1);
for i = 1:params.N
    z0(4*(i-1)+1) = x(i);
    z0(4*(i-1)+2) = y(i);
    z0(4*(i-1)+3) = u(i);
    z0(4*(i-1)+4) = v(i);
end

% Solve using variable time step RK45, interpolate to fixed points in time
sol = ode45(@(t, z) computeOde(t, z, params), [ti, tf], z0, odeopts); 
z = deval(sol, tout);

% Reshape back into vectors and plot
figure(1); clf; hold on;

x = [];
y = [];
dxdt = [];
dydt = [];
for i = 1:params.N
    start = 4*(i-1);
    x = [x; z(start+1,:)];
    y = [y; z(start+2,:)];
    
    % Do some plotting
    plot(x(i,:), y(i,:));
end
axis equal;

display('Done...')

%--------------------------------------------------------------------------
% ODE - RHS
%--------------------------------------------------------------------------
function dz = computeOde(t, z, params)

a = params.a;
N = params.N;

dz = zeros(4*N,1);

x = z(1:4:4*N);
y = z(2:4:4*N);
dxdt = z(3:4:4*N);
dydt = z(4:4:4*N);

for i = 1:N
    xi = [x(i); y(i)];
    dxdti = [dxdt(i); dydt(i)];
    UA = sqrt(dot(dxdti,dxdti));
    thetaA = atan2(dxdti(2), dxdti(1));
    
    term1 = zeros(2,1);
    UIvec = zeros(2,1);
    dudx = 0;
    dudy = 0;
    
    for j = 1:N
        if j ~= i
            xj = [x(j); y(j)];
            dxdtj = [dxdt(j); dydt(j)];
            deltaji = xj - xi;
            deltaji2 = dot(deltaji, deltaji);
            deltaji4 = deltaji2^2;
            deltaji6 = deltaji2^3;
            
            % Calculate term 1
            term1a = 2*dot(dxdtj-dxdti,deltaji)/deltaji4 * dxdtj;
            term1b = 2*dot(dxdtj,dxdtj-dxdti)/deltaji4 * deltaji;
            term1c = 2*dot(dxdtj,deltaji)/deltaji4 * (dxdtj - dxdti);
            term1d = -8*dot(dxdtj,deltaji)*dot(dxdtj-dxdti,deltaji)/deltaji6 * deltaji;
            term1 = term1 + term1a + term1b + term1c + term1d; % 2 components
            
            % Calculate term 2
            
            % Accumulate UI
            UIvec = UIvec + -dxdtj/deltaji2 + (2*dot(dxdtj,deltaji)/deltaji4 * deltaji);
            
            % Accumulate dudx, dudy to calculate thetaE and E
            term2a = (- 6*deltaji(1)*dxdtj(1) - 2*deltaji(2)*dxdtj(2))/deltaji4;
            term2b = (8*dot(dxdtj, deltaji)*deltaji(1)^2)/deltaji6;
            dudx = dudx + term2a + term2b;
            
            term3a = (- 2*deltaji(1)*dxdtj(2) - 2*(deltaji(2))*dxdtj(1))/deltaji4;
            term3b = (8*dot(dxdtj,deltaji)*deltaji(1)*deltaji(2))/deltaji6;
            dudy = dudy + term3a + term3b;
        end
    end
    
    UI = sqrt(dot(UIvec,UIvec));
    thetaI = atan2(UIvec(2),UIvec(1));
    
    E = sqrt(dudx*dudx + dudy*dudy);
    thetaE = 0.5*atan2(dudy, dudx);
    
    term2x = -2*pi*a^2*E*(UI*cos(2*thetaE-thetaI) - UA*cos(2*thetaE-thetaA));
    term2y = 2*pi*a^2*E*(UI*sin(2*thetaE-thetaI) - UA*sin(2*thetaE-thetaA));
    
    ddxdt = term1(1) + term2x;
    ddydt = term1(2) + term2y;
    
    start = 4*(i-1);
    
    dz(start+1) = dxdti(1);
    dz(start+2) = dxdti(2);
    dz(start+3) = ddxdt;
    dz(start+4) = ddydt;
end

%--------------------------------------------------------------------------
% MASS MATRIX - LHS
%--------------------------------------------------------------------------

function m = massMatrix(t, z, params)

a = params.a;
rhoc = params.rhoc;
rhof = params.rhof;
N = params.N;

x = z(1:4:4*N);
y = z(2:4:4*N);

m = zeros(4*N);

factor = 2*pi*rhof*a^2;

for i = 1:N
    xi = [x(i); y(i)];
    
    for j = 1:N    
        xj = [x(j); y(j)];
        deltaji = xj - xi;        
        deltaji2 = dot(deltaji, deltaji);
        deltaji4 = deltaji2^2;
        
        if i == j
            m(4*(i-1)+1,4*(j-1)+1) = 1; % dx = dx 
            m(4*(i-1)+2,4*(j-1)+2) = 1; % dy = dy
            m(4*(i-1)+3,4*(j-1)+3) = pi*a^2*(rhoc + rhof); % mass and added mass
            m(4*(i-1)+4,4*(j-1)+4) = pi*a^2*(rhoc + rhof);
        else
            m(4*(i-1)+3,4*(j-1)+3) = factor*(2.0/deltaji4)*deltaji(1)*deltaji(1) - 1.0/deltaji2; % ijTerm for ddXiXj
            m(4*(i-1)+3,4*(j-1)+4) = factor*(2.0/deltaji4)*deltaji(1)*deltaji(2); % ijTerm for ddXiYj
            m(4*(i-1)+4,4*(j-1)+3) = factor*(2.0/deltaji4)*deltaji(2)*deltaji(1); % ijTerm for ddYiXj
            m(4*(i-1)+4,4*(j-1)+4) = factor*(2.0/deltaji4)*deltaji(2)*deltaji(2) - 1.0/deltaji2; % ijTerm for ddYiYj
        end
    end
    
end
