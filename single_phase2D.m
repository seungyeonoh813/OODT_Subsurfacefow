%% solve for Quarter five-spot with compressible flow
%  \farc{\partial(\phi \rho)}{\partial t} + \nabla \cdot (\rho v) = q,
%  v = -K/\mu (\nabla p -g\rho \nabla z)

clear all;
close all;

addpath(genpath('../mrst'));
run startup
% Add modules
mrstModule add ad-core ad-blackoil


logk = 1+ (2-1)*rand(1000,1);
q_inj = 100 + (200-100) *rand (1000,1);
q_prod = zeros(1000,100);

for n = 1:1000
    perm = 10^(logk(n)) * milli * darcy;
    inj_rate = q_inj(n) * meter^3/day;  
    [G, schedule, ws,states] = simulate_single_phase_2D(perm,inj_rate);
    filename = sprintf('data/states_%d.mat', n);
    % save(filename, "states","ws");
    for i = 1:numel(ws)
        for j = 1:numel(ws{i})
            if strcmp(ws{i}(j).name, 'Producer')
                q_prod(n,i) = ws{i}(j).qWs;
            end
        end
    end
end
% save('data/QoI.mat','logk','q_inj','q_prod');

% Plot final pressure
figure;
plotCellData(G, states{end}.pressure);
title('Final Pressure with Wells');
axis equal tight; colorbar;

% Extract injection/production rate over time
prod_rates = zeros(numel(ws), 1);
inj_rates = zeros(numel(ws), 1);
for i = 1:numel(ws)
    for j = 1:numel(ws{i})
        if strcmp(ws{i}(j).name, 'Injector')
            inj_rates(i) = ws{i}(j).qWs;
        end
        if strcmp(ws{i}(j).name, 'Producer')
            prod_rates(i) = ws{i}(j).qWs;
        end
    end
end
% Plot injection rate
figure;
plot(cumsum(schedule.step.val)/day, inj_rates * day);
xlabel('Time [days]');
ylabel('Injection Rate [m^3/day]');
title('SW Producer Water Rate');
grid on;


% Plot production rate
figure;
plot(cumsum(schedule.step.val)/day, prod_rates * day);
xlabel('Time [days]');
ylabel('Production Rate [m^3/day]');
title('NE Producer Water Rate');
grid on;




function [G, schedule, ws,states] = simulate_single_phase_2D(perm,inj_rate)
    % Create grid and rock
    G = cartGrid([32, 32], [500, 500]); G = computeGeometry(G);
    rock = makeRock(G, perm, 0.2);
    
    % Fluid: Compressible water
    fluid = initSimpleADIFluid('phases', 'W', ...
                               'mu', [1 * centi * poise], ...
                               'rho', [1014 * kilogram/meter^3], ...
                               'c', [1e-5 / barsa]);  % compressibility
    
    % Set up model
    model = WaterModel(G, rock, fluid);
    
    % Set initial state
    state0 = initResSol(G, 100 * barsa, [1]);  % 100 bar and 100% water
    
    % Define wells
    W = [];
    % Injector at SW corner (cell 1) with specified rate (input)
    W = addWell(W, G, rock, 1, 'Type', 'rate', ...
                     'Val', inj_rate, ...
                     'Name', 'Injector', ...
                     'Sign', 1, ...
                     'comp_i', [1]);
    
    % % Injector at SW corner (cell 1) with BHP control
    % inj_bhp = 110 * barsa;   % bottom-hole pressure
    % W = addWell(W, G, rock, 1, 'Type', 'bhp', ...
    %                  'Val', inj_bhp, ...
    %                  'Name', 'Injector', ...
    %                  'Sign', 1, ...
    %                  'comp_i', [1]);
    
    
    % Producer at NE corner (last cell) with BHP control
    prod_bhp = 90 * barsa;  % bottom-hole pressure
    W = addWell(W, G, rock, G.cells.num, 'Type', 'bhp', ...
                     'Val', prod_bhp, ...
                     'Name', 'Producer', ...
                     'Sign', -1, ...
                     'comp_i', [1]);
    
    % Set schedule
    T = 30 * day;
    dt = repmat(T / 100, 100, 1);  % 500 uniform timesteps
    schedule = simpleSchedule(dt, 'W', W);
    
    % Simulate
    [ws, states, reports] = simulateScheduleAD(state0, model, schedule);
end




