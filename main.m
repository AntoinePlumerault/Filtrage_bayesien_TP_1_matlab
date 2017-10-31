%% DATA 

% map 
map = load('mnt.data');
[N1, N2] = size(map);

% Real trajectory
load('traj.mat', 'rtrue', 'vtrue');
nmax = size(rtrue, 2); 

% Inertial measurments
load('ins.mat', 'a_INS');

% Altitude measurments
load('alt.mat', 'h_ALT');

% Constants initialization
X1MIN = -10000;
X1MAX = 10000;
X2MIN = -10000;
X2MAX = 10000;

r_0 = [-6000;-2000];
v_0 = [120;0];

delta = 1; 

sigma_r_0 = 100;
sigma_v_0 = 10;

sigma_INS = 7;
sigma_ALT = 10;
sigma_BAR = 20;

%% PREDICTION

% Initialization of the display
set(gcf, 'DoubleBuffer', 'on')
figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1);
% map
imagesc([X1MIN,X1MAX], [X2MIN,X2MAX], transpose(map));
hold on;
colormap('default');
axis square;
axis off;

% Initialization of prediction using inertial data only
r_INS(:,1) = r_0;
v_INS(:,1) = v_0;

% Initization of the partcles for the adaptive SIR filter 
N = 1000; % Number of particles
c = 0.01;
N_thr = c*N;
r_part = randn(2,N) * sigma_r_0;
v_part = randn(2,N) * sigma_v_0;
w_part = ones(1,N) / N;
r_SIR = zeros(2,nmax);
r_SIR(:,1) = r_0;

alt = zeros(1,nmax);
alt(1) = map(floor((rtrue(1,1) - X1MIN) / 20000.0 * N1), ...
             floor((rtrue(2,1) - X2MIN) / 20000.0 * N2));
for n = 2:nmax
    % Prediction using inertial data only
    r_INS(:,n) = r_INS(:,n-1) + delta * v_INS(:,n-1);
    v_INS(:,n) = v_INS(:,n-1) + delta * a_INS(:,n-1);
    
    % Predicton using the adaptive SIR filter
    N_eff = 1/sum(w_part.^2);
    if N_eff < N_thr
        % Resampling
        new_r_part = zeros(2,N);
        new_v_part = zeros(2,N);
        new_w_part = ones(1,N)/N;
        
        [cumulative_p, cumulative_i] = sort(w_part);
        cumulative_p = cumsum(cumulative_p);
        for k = 1:N
            r = rand;
            j = 2;
            while cumulative_p(j) <= r
                j = j + 1;
            end
            index = cumulative_i(j-1);
            new_r_part(:,k) = r_part(:, index);  
            new_v_part(:,k) = v_part(:, index);
        end
        
        r_part = new_r_part;
        v_part = new_v_part;
        w_part = new_w_part;
    end
    
    % Actualization
    r_part = r_part + delta * v_part;
    v_part = v_part - delta * randn(2,N) * sigma_INS;
    
    % Correction
    for j = 1:N
        x = floor((r_INS(1,n) + r_part(1,j) - X1MIN) / 20000.0 * N1);
        y = floor((r_INS(2,n) + r_part(2,j) - X2MIN) / 20000.0 * N2);
        if x < 1 || x > N1 || y < 1 || y > N2
            w_part(:,j) = 0; % if a particle is out of the map we assign it a weight of zero
        else
            w_part(:,j) = w_part(:,j) * exp(-1.0/(2.0 * (sigma_BAR^2 + sigma_ALT^2)) * (h_ALT(n) - map(x, y))^2);
        end
    end
    w_part = w_part / sum(w_part);
    
    % Prediction
    r_SIR(:,n) = r_INS(:,n) + sum(r_part .* repmat(w_part,2,1),2);
    
    % display
    %real_alt
    alt(n) = map(floor((rtrue(1,n) - X1MIN) / 20000.0 * N1), ...
                 floor((rtrue(2,n) - X2MIN) / 20000.0 * N2));
    clf;
    
    subplot(1,2,1);
    imagesc([X1MIN, X1MAX], [X2MIN, X2MAX], transpose(map));
    hold on;
    axis square;
    axis off;
    plot(r_SIR(1,1:n), r_SIR(2,1:n), 'm-');
    plot(r_INS(1,1:n), r_INS(2,1:n), 'g-');
    plot(rtrue(1,1:n), rtrue(2,1:n), 'r-');
    scatter(r_INS(1,n)+r_part(1,:), r_INS(2,n)+r_part(2,:), 'k+');
    title('Trajectories of the plane on the relief map');
    legend(['SIR (c = ',num2str(c),')'],'Inertial estimation', 'Real trajectory', ['Particles (N = ', num2str(N), ')'])
    
    subplot(1,2,2);
    plot(0:n-1,alt(1:n));
    title('Real relief and measured relief over time');
    xlabel('time');
    ylabel('relief');
    axis square;
    hold on;
    scatter(0:n-1,h_ALT(1:n), 'x');
    legend('Real relief', 'Measured relief');
    
    drawnow;
    
end

%% OLD VERSION OF THE CODE
% for n = 2:nmax
%     r_INS(:,n) = r_INS(:,n-1) + delta * v_INS(:,n-1);
%     v_INS(:,n) = v_INS(:,n-1) + delta * a_INS(:,n-1);
% end
% 
% %% ALTITUDE DURING FLIGHT
% alt = zeros(1,nmax);
% 
% for n = 1:nmax
%     alt(n) = map(floor((rtrue(1,n) - X1MIN) / 20000.0 * N1), ...
%                  floor((rtrue(2,n) - X2MIN) / 20000.0 * N2));
% end
% 
% %% SIS FILTER
% N = 1000;
% r_part = randn(2,N) * sigma_r_0;
% v_part = randn(2,N) * sigma_v_0;
% 
% w_part = ones(1,N);
% 
% r_SIS = zeros(2,nmax);
% 
% for i = 1:nmax
%     % mutation
%     r_part = r_part + delta * v_part;
%     v_part = v_part - delta * randn(2,N) * sigma_INS;
%     
%     %correction
%     for j = 1:N
%         w_part(:,j) = w_part(:,j)* exp(-1.0/(2.0 * (sigma_BAR^2 + sigma_ALT^2)) * (h_ALT(i) - map(min(N1,max(1,floor((r_INS(1,i) + r_part(1,j) - X1MIN) / 20000.0 * N1))), ...
%                                                                                      min(N2,max(1,floor((r_INS(2,i) + r_part(2,j) - X2MIN) / 20000.0 * N2)))))^2);
%     end
%     w_part = w_part / sum(w_part);
%    
%     r_SIS(:,i) = r_INS(:,i) + sum(r_part .* repmat(w_part,2,1),2);
% end
% 
% %% SIR FILTER
% N = 1000;
% c = 0.01;
% N_thr = c*N;
% r_part = randn(2,N) * sigma_r_0;
% v_part = randn(2,N) * sigma_v_0;
% 
% w_part = ones(1,N)/N;
% 
% r_SIR = zeros(2,nmax);
% 
% set(gcf, 'DoubleBuffer', 'on')
% figure;
% for i = 1:nmax
%     % réechantillonnage
%     N_eff = 1/sum(w_part.^2);
%     if N_eff < N_thr
%         new_r_part = zeros(2,N);
%         new_v_part = zeros(2,N);
%         new_w_part = zeros(1,N);
%         [cumulative_p, cumulative_i] = sort(w_part);
%         cumulative_p = cumsum(cumulative_p);
%         for k = 1:N
%             r = rand;
%             j = 2;
%             while cumulative_p(j) <= r
%                 j = j + 1;
%             end
%             index = cumulative_i(j-1);
%             new_r_part(:,k) = r_part(:, index);  
%             new_v_part(:,k) = v_part(:, index);
%             new_w_part(:,k) = w_part(:, index);
%         end
%         r_part = new_r_part;
%         v_part = new_v_part;
%         w_part = new_w_part;
%     end
%     
%     % mutation
%     r_part = r_part + delta * v_part;
%     v_part = v_part - delta * randn(2,N) * sigma_INS;
%     
%     %correction
%     for j = 1:N
%         w_part(:,j) = w_part(:,j) * exp(-1.0/(2.0 * (sigma_BAR^2 + sigma_ALT^2)) * (h_ALT(i) - map(min(N1,max(1,floor((r_INS(1,i) + r_part(1,j) - X1MIN) / 20000.0 * N1))), ...
%                                                                                      min(N2,max(1,floor((r_INS(2,i) + r_part(2,j) - X2MIN) / 20000.0 * N2)))))^2);
%     end
%     w_part = w_part / sum(w_part);
%    
%     r_SIR(:,i) = r_INS(:,i) + sum(r_part .* repmat(w_part,2,1),2);
%     
%     % affichage
%     clf;
%     imagesc([X1MIN, X1MAX], [X2MIN, X2MAX], transpose(map));
%     hold on;
%     axis square;
%     axis off;
%     plot (r_SIR(1, 1:i), r_SIR(2, 1:i), 'y-');
%     scatter(r_INS(1,i)+r_part(1,:), r_INS(2,i)+r_part(2,:), 'k+');
%     drawnow;
%     
% end
% %% AFFICHAGE
% 
% figure
% % map
% imagesc([X1MIN, X1MAX], [X2MIN, X2MAX], transpose(map));
% hold on;
% colormap('default');
% axis square;
% axis off;
% 
% % real trajectory
% plot(rtrue(1,:), rtrue(2,:), 'r-');
% 
% % trajectory estimated with inertial data
% plot(r_INS(1,:), r_INS(2,:), 'm-');
% 
% plot(r_SIS(1,:), r_SIS(2,:), 'g-');
% 
% plot(r_SIR(1,:), r_SIR(2,:), 'y-');
% 
% % altitude
% figure
% plot(0:100,alt);
% hold on;
% scatter(0:100,h_ALT, 'x');
% 
% 
