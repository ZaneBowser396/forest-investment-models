clear, clc;

IFL_DATA = readtable('tabulated_ifl_data.xlsx');
% IFL_DATA = IFL_DATA(1:10, :); % To select only the best static sorted IFLs

I = 2e11; % Resources allocated per timestep

z = 0.3; % Species preserved scalar
T = 50; % Time period in years
dt = 0.1; % Timestep Size
tol = 1e-6; % iteration tolerance
maxiters = 10; % number of iterations

SpeciesProtected = optimal_control(IFL_DATA, I, z, T, dt, tol, maxiters);


function SpeciesProtected = optimal_control(IFL_DATA, I, z, T, dt, tol, maxiters)

A_2016 = IFL_DATA.A_2016; % Area in 2016
A_2020 = IFL_DATA.A_2020; % Area in 2020
T_i = A_2016; % Total Area of each IFL

B = IFL_DATA.Species_Richness; % Number of species per IFL
c = IFL_DATA.Travel_Time; % travel time (cost) for each IFL
d = IFL_DATA.LossRate; % Loss rate of available land per IFL
p = IFL_DATA.LocalSupport; % Chance of local support based on loss rate and carbon
q = IFL_DATA.ChanceOfSuccess; % Chance of less than two incidents
N = length(A_2016); % Number of IFLs in consideration

Tvec = 0:dt:T; % Timestep vector for time considered
Nt = length(Tvec); % Number of timesteps

AL = zeros(N, Nt); % Available Land over time per IFL
PL = zeros(N, Nt); % Protected Land over time per IFL
LL = zeros(N, Nt); % Lost Land over time per IFL

AL(:, 1) = A_2020; % Initially 2020 land is available
LL(:, 1) = (A_2016-A_2020); % Lost the difference since 2016
PL(:, 1) = 0; % No land protected

u = zeros(N, Nt); % Percentage of I to allocate to each IFL per timestep
u(:, :) = 1/N; % begin with equal allocation to all over all time
alpha = (p .* q * I) ./ (c .* T_i); % Static ranking coefficient

% costate variables for each IFL for all timesteps
% No benefit to conserve lost land, so always zero (ignored)
lambda_A = zeros(N, Nt); % scaling the value of available land in optimising Hamiltonian
lambda_P = zeros(N, Nt); % scaling the value of protected land in optimising Hamiltonian

% Iterate until a solution is found
% Seems to flip between solutions with this, so we run one iteration
for iter = 1:maxiters
    u_prev = u;

    % Calculate Land values at each time based on resource allocation in
    % differential equations
    for k = 1:Nt-1
        A = AL(:, k);
        P = PL(:, k);
        L = LL(:, k);

        Lost = min(A, dt * d .* A);
        A = A - Lost;
        L = L + Lost;

        u_k = u(:, k);
        Prot = min(A, dt * alpha .* u_k);
        A = A - Prot;
        P = P + Prot;

        AL(:, k+1) = A;
        PL(:, k+1) = P;
        LL(:, k+1) = L;
    end

    % Calculate costate variables based on land values
    % Determines how much benefit each gives
    A_T = AL(:, end);
    P_T = PL(:, end);
    common_term = z * B .* ((A_T + P_T)./T_i).^(z - 1)./T_i;
    lambda_A(:, Nt) = common_term; % .* alpha .* u(:, end);
    lambda_P(:, Nt) = common_term; % .* d .* A_T;

    % Backwards iteration to calculate new costate variables across time
    for t = Nt:-1:2
        lambda_A(:, t-1) = lambda_A(:, Nt) .* exp(d .* (t - Nt));
        lambda_P(:, t-1) = lambda_P(:, t);
    end

    % Using new costate variables determine Delta, score of how much
    % benefit allocating resources gives, then give all resources at that
    % timestep to biggest Delta IFL
    for t = 1:Nt
        Delta = (lambda_P(:, t) - lambda_A(:, t));
        Delta(AL(:, t) <= 0) = -Inf;
        u(:, t) = 0;
        [~, idx] = max(Delta);
        u(idx, t) = 1;
    end

    % Check if u has converged to optimal solution or if reached maxiters,
    % and if so break the loop
    du = max(max(abs(u - u_prev)));
    if du < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
    if iter == maxiters
        fprintf('Reached maximum iterations without full convergence.\n');
    end
end

% Calculate final land variables based on latest resource allocation u
for k = 1:Nt-1
    A = AL(:, k);
    P = PL(:, k);
    L = LL(:, k);

    Lost = min(A, dt * d .* A);
    A = A - Lost;
    L = L + Lost;

    u_k = u(:, k);
    Prot = min(A, dt * alpha .* u_k);
    A = A - Prot;
    P = P + Prot;

    AL(:, k+1) = A;
    PL(:, k+1) = P;
    LL(:, k+1) = L;
end

SpeciesProtected = sum(B.*((PL(:, Nt) + AL(:, Nt))./T_i).^z);

% Plot results: changes in land over time, with respect to changing u
figure(1), clf;
time = Tvec;

subplot(2,1,1), hold on;
for i=1:length(A_2020)
    plot(time, AL(i, :)/T_i(i), 'Color', [0 0 1 0.2], 'LineWidth', 1);
    plot(time, PL(i, :)/T_i(i), 'Color', [0 1 0 1], 'LineWidth', 1);
    plot(time, LL(i, :)/T_i(i), 'Color', [1 0 0 0.2], 'LineWidth', 1);
end
legend({'Available', 'Protected', 'Lost'}, 'FontSize', 16);
ylabel('Mean Land Fraction', 'FontSize', 16);

subplot(2, 1, 2), hold on;
for i = 1:length(A_2020)
    plot(time, u(i, :), 'Color', [0 0 0 0.2], 'LineWidth', 1);
end
ylim([0 1]);
xlabel('Time', 'FontSize', 16);
ylabel('Investment', 'FontSize', 16);

title(sprintf('Total Species Protected: %.2f', SpeciesProtected));

% Plot cumulative u over time
figure
legend_entries = {};
plot_handles = [];
hold on;
for i=1:length(A_2020)
    cum_u = zeros(Nt,1);
    for t = 1:Nt
        cum_u(t) = sum(u(i,1:t));
    end
    if sum(cum_u) > 1
        h = plot(Tvec, cum_u, 'LineWidth', 1.5);
        plot_handles = [plot_handles, h];
        label = replace(string(IFL_DATA.IFL_ID{i}), "_", "\_");
        legend_entries{end+1} = sprintf("%s", label);
    end
end

if ~isempty(legend_entries)
    legend(plot_handles, legend_entries, 'Location', 'best');
end
xlabel('Time');
ylabel('Total Allocated Resources');
title('Cumulative Allocated Resources per IFL over time');
hold off;


end