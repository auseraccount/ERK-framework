clear;
%% This script verifies the positivity of the energy stability 
% operators obtained by first to third order ERK schemes
% Reference: Haifeng Wang, Hong Zhang, Jingwei Sun, Lingyan Tang, Global-in-time 
% energy estimate of a third-order exponential Runge--Kutta framework 
% for the phase field crystal equation

syms z un un1 un2 un3
ERK_opt = 33; % Switch ERK schemes
switch ERK_opt 
    case 11 % ERK(1, 1) based on underlying Euler scheme
        % Butcher tableau: A = [0]; b = [1];
        stage = 1;
        c0 = 0; c1 = 1;
        phi11 = (exp(-c1*z) - 1)/(-c1*z);
        a10 = phi11;
    case 22 % ERK(2, 2) based on underlying RK(2, 2) scheme
        % Butcher tableau: A = [0 0; 1/2 0]; B = [0 1]';
        stage = 2;
        c0 = 0; c1 = 1/2; c2 = 1;
        phi11 = (exp(-c1*z) - 1)/(-c1*z);
        phi12 = (exp(-c2*z) - 1)/(-c2*z);
        a10 = 1/2*phi11; a20 = 0; a21 = phi12;
    case 33 % ERK(3, 3) based on undelryin RK(3, 3) scheme
        % Butcher tableau: A = [0 0 0; 1/3 0 0; 0 2/3 0]; b = [1/4 0 3/4]';
        stage = 3;
        c0 = 0; c1 = 1/3; c2 = 2/3; c3 = 1; phi0 = exp(-z);
        phi11 = (exp(-c1*z)-1)/(-c1*z);
        phi12 = (exp(-c2*z)-1)/(-c2*z);
        phi13 = (exp(-c3*z)-1)/(-c3*z);
        a10 = c1*phi11; a20 = 0; a21 = c2*phi12;
        a30 = 1/4*phi0;
        a31 = 0;
        a32 = phi13 - 1/4*phi0;
end
switch stage
    case 1
        A = [a10];
    case 2
        A = [a10 0; a20 a21];
    case 3
        A = [a10 0 0; a20 a21 0; a30 a31 a32];
end

%% Verification of order conditions 
% Ref: H. Zhang, X. Qian, J. Xia, S. Song, Efficient inequality-preserving integrators for differential equations satisfying forward Euler
% conditions, ESAIM: Mathematical Modelling and Numerical Analysis 57 (3) (2023) 1619¨C1655.
if stage > 0
    hata10 = a10;
    hatb = hata10;
    c = [c0 c1]';
    if stage < 2
        fprintf('b * 1 - 1=\n');
        taylor(sum(hatb) - 1, z, 0, 'order', 3)
    end
end
if stage > 1
    hata10 = a10;
    hata20 = a20 + z*a21*hata10;
    hata21 = a21;
    hatb = [hata20 hata21 ];
    hatA = [0 0; hata10 0 ];
    c = [c0 c1]';
    if stage < 3
        fprintf('b * 1 - 1=\n');
        taylor(sum(hatb) - 1, z, 0, 'order', 3)
        fprintf('b * c -1/2 =\n');
        taylor(hatb*c - 1/2, z, 0, 'order', 3)
    end
end
if stage > 2
    hata30 = a30 + z*a31*hata10 + z*a32*hata20;
    hata31 = a31 + z*a32*hata21;
    hata32 = a32;
    if stage < 4
        hatb = [hata30 hata31 hata32];
        hatA = [0 0 0; hata10 0 0; hata20 hata21 0];
        c = [c0 c1 c2]';
        fprintf('b * 1 - 1=\n');
        taylor(sum(hatb) - 1, z, 0, 'order', 4)
        fprintf('b * c -1/2 =\n');
        taylor(hatb*c - 1/2, z, 0, 'order', 3)
        fprintf('b * c^2 -1/3 =\n');
        taylor(hatb*c.^2 - 1/3, z, 0, 'order', 3)
        fprintf('b * A * c -1/6=\n');
        taylor(hatb*hatA*c - 1/6, z, 0, 'order', 4)
    end
end

xv =  0.001:0.1:100;

%% Verification of the positive definiteness of energy stability operators
% z := \tau * L_\kapap in the ERK framework
%% Stage 1
phi01 = exp( c1*(-z) );
tau_Nkappa0 = ( a10 )^(-1)*( un1 - phi01*un ); % tau * N_kappa( u_{n, 0} )
tmp = coeffs(collect( -z*un1 + tau_Nkappa0 ), un1);
Delta_11 = inline(tmp(end));
subplot(2,2,1); hold on; title('One stage');
plot(xv, Delta_11(xv), 'linewidth', 1.);
xlabel('$z$','interpreter', 'latex'); ylabel('Energy stability function');
set(gca, 'fontname', 'Times New Roman');
legend('$\Delta_{1,1}$', 'interpreter', 'latex');
fprintf('min(Delta_11) = %d\n', min(Delta_11(xv)));

%% stage 2
if stage > 1
    phi02 = exp( c2*(-z) );
    tau_Nkappa1 = ( a21 )^(-1)*( un2 - phi02*un - a20*tau_Nkappa0 );
    tmp = coeffs(collect( -z*un2 + tau_Nkappa1 ), un2);
    Delta_22 = tmp(end);
    tmp = coeffs(collect( -z*un2 + tau_Nkappa1 - Delta_22*(un2 - un1) ), un1);
    Delta_21 = inline(tmp(end)); Delta_22 = inline(Delta_22);
    subplot(2,2,2); hold on; title('Two stage');
    plot(xv, Delta_11(xv) - 1/2*Delta_21(xv), 'linewidth', 1);
    plot(xv, Delta_21(xv), 'linewidth', 1);
    plot(xv, Delta_22(xv)- 1/2*Delta_21(xv), 'linewidth', 1);
    xlabel('$z$','interpreter', 'latex'); ylabel('Energy stability function');
    set(gca, 'fontname', 'Times New Roman');
    hl = legend('$\Delta_{1,1} - \frac12 \Delta_{2,1}$', '$\Delta_{2,1}$', '$\Delta_{2,2} - \frac12 \Delta_{2,1}$', 'location', 'northeast');
    set(hl, 'interpreter', 'latex');
    fprintf('min(Delta_22) = %d\n', min(Delta_22(xv)));
    fprintf('min(Delta_21) = %d\n', min(Delta_21(xv)))
    fprintf('min(Delta_11 - 1/2*Delta_21) = %d\n', min(Delta_11(xv) - 1/2*Delta_21(xv)))
    fprintf('min(Delta_22 - 1/2*Delta_21) = %d\n', min(Delta_22(xv) - 1/2*Delta_21(xv)))
end
if stage > 2
    %% stage 3
    phi03 = exp(c3*(-z));
    tau_Nkappa2 = (a32)^(-1)*(un3 - phi03*un - a30*tau_Nkappa0 - a31*tau_Nkappa1);
    tmp = coeffs(collect(-z*un3 + tau_Nkappa2), un3);
    Delta_33 = tmp(end);
    tmp = coeffs(collect(-z*un3 + tau_Nkappa2 - Delta_33*(un3 - un2)), un2);
    Delta_32 = tmp(end);
    tmp = coeffs(collect(-z*un3 + tau_Nkappa2 - Delta_33*(un3 - un2) -  Delta_32*(un2-un1)), un1);
    Delta_31 = tmp(end); Delta_31 = inline(Delta_31); % May be zero, this one is safe
    tmp = coeffs(collect(-z*un3 + tau_Nkappa2), un);
    D310_v = inline(-tmp(end)); Delta_33 = inline(Delta_33); Delta_32 = inline(Delta_32);
    subplot(2,2,3); hold on; title('Three stage');
    plot(xv, Delta_11(xv) - 1/2*Delta_21(xv) - 1/2*Delta_31(xv), 'linewidth', 1); hold on;
    plot(xv, Delta_21(xv), 'linewidth', 1);
    plot(xv, Delta_22(xv) - 1/2*Delta_21(xv) - 1/2*Delta_32(xv), 'linewidth', 1);
    plot(xv, Delta_31(xv), 'linewidth', 1);
    plot(xv, Delta_32(xv), 'linewidth', 1);
    plot(xv, Delta_33(xv) - 1/2*Delta_31(xv) - 1/2*Delta_32(xv), 'linewidth', 1);
    xlabel('$z$','interpreter', 'latex'); ylabel('Energy stability function');
    set(gca, 'fontname', 'Times New Roman');
    hl = legend('$\Delta_{1,1} - \frac12 \Delta_{2,1} - \frac12 \Delta_{3,1}$', ...
        '$\Delta_{2,1}$', '$\Delta_{2,2} - \frac12 \Delta_{2,1} - \frac12 \Delta_{3,2}$', ...
        '$\Delta_{3,1}$', '$\Delta_{3,2}$', '$\Delta_{3,3} - \frac12 \Delta_{3,1} - \frac12 \Delta_{3,2}$', ...
        'location', 'northeast');
    set(hl, 'interpreter', 'latex');
    
    fprintf('min(Delta_11 - 1/2*Delta_21 - 1/2*Delta_31) = %d\n', ...
        min(Delta_11(xv) - 1/2*Delta_21(xv) - 1/2*Delta_31(xv)));
    fprintf('min(Delta_21) = %d\n', min(Delta_21(xv)));
    fprintf('min(Delta_22 - 1/2*Delta_21 - 1/2*Delta_32) = %d\n', ...
        min(Delta_22(xv) - 1/2*Delta_21(xv) - 1/2*Delta_32(xv)));
    fprintf('min(Delta_31) = %d\n', min(Delta_31(xv)));
    fprintf('min(Delta_32) = %d\n', min(Delta_32(xv)));
    fprintf('min(Delta_33 - 1/2*Delta_31 - 1/2*Delta_32)= %d\n', ...
        min(Delta_33(xv) - 1/2*Delta_31(xv) - 1/2*Delta_32(xv)));
end