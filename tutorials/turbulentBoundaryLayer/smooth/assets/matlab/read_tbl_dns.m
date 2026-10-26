%% read_tbl_dns.m -- read and plot the smooth-wall ZPG boundary-layer DNS
%
% Data: tbl_dns_smooth.mat (one file, put it next to this script).
% Run:  >> read_tbl_dns
%
% THE CASE
%   Spatially developing, zero-pressure-gradient, incompressible turbulent
%   boundary layer over a smooth flat plate (DNS, solver mobydiff, 2nd-order
%   finite differences). Same domain, grid and trip as the CaNS reference DNS:
%     domain   Lx x Ly x Lz = 650 x 100 x 26 (in delta*_in), periodic in z
%     grid     3200 x 384 x 136; uniform in x and z (dx = 0.203, dz = 0.191),
%              tanh-stretched in y (first cell centre y = 0.0052, y+ ~ 0.25)
%     inflow   laminar Blasius profile at Re_delta* = 450 (Re_theta = 174)
%     trip     transition is FORCED, not natural (see TRIP below). There is
%              NO free-stream turbulence (Tu = 0).
%     top      y = Ly: an outflow (fixed-pressure) face through which the
%              displacement entrainment leaves; far outside the layer
%              (delta99 <= 15 here). Outlet at x = Lx: the same.
%   Statistics: averages over z and time (1250 time units, 65000 samples).
%   Stored for x <= 550 only: the last 100 units ahead of the outlet are an
%   adjustment zone and are not usable.
%
% TRIP (how transition is triggered)
%   The laminar Blasius inflow is tripped by a random volume force added to
%   the WALL-NORMAL momentum equation only (Schlatter & Orlu 2012, JFM 710,
%   "Turbulent boundary layers at moderate Reynolds numbers: inflow length
%   and tripping effects"; the same implementation as the CaNS and SIMSON
%   reference DNS). In the units of this file:
%
%     dv/dt + ... = f_v(x,y,z,t)
%     f_v = exp( -((x - x0)/lx)^2 - (y/ly)^2 ) * (A/nmodes) * g(z,t)
%
%   - Envelope: a Gaussian bump at the wall, x0 = 10, lx = 4, ly = 1 (so it
%     acts over roughly 2 < x < 18 and y < 2, i.e. inside the laminar layer,
%     whose delta99 is ~3 there). Fluid upstream of x ~ 2 is unforced
%     Blasius flow.
%   - Spanwise signal: a Fourier series over the periodic span Lz = 26 with
%     random phases phi_m (uniform in [0, 2 pi]) and fixed amplitudes,
%       g_k(z) = sqrt(2) cos(phi_0) + sum_{m=1}^{nmodes/2} 2 cos(2 pi m z/Lz + phi_m)
%     with nmodes = 16, i.e. 8 harmonics, spanwise wavelengths 26/m =
%     26 ... 3.25 delta*_in. Its spanwise rms is about sqrt(nmodes), so the
%     effective force amplitude is A/sqrt(nmodes) = 0.18854/4 = 0.047
%     (in U_inf^2/delta*_in).
%   - Time dependence: a new random realisation g_k is drawn every
%     ts = 4 time units (delta*_in/U_inf); between draws the force blends
%     smoothly from one to the next,
%       g(z,t) = (1 - b) g_k(z) + b g_{k+1}(z),  b = 3p^2 - 2p^3,
%       p = t/ts - k,  k = floor(t/ts),
%     so the forcing is continuous in time (a C^1 blend). The steady part of
%     the original formulation (amplitude A_s) is switched OFF here (A_s = 0).
%   - Effect: the layer becomes turbulent within ~60 delta*_in of the trip
%     (c_f overshoot, 0.0062, at x ~ 70, Re_theta ~ 355; H falls from 2.59
%     below 1.55 at the same station). The transition location is therefore set by
%     the trip, not by a receptivity/instability process: a transition model
%     should be compared from the trip onward with a prescribed transition
%     onset, or on the fully turbulent part (Re_theta > ~400).
%   Parameters as in the solver input (input_ini in the file): trip_amp =
%   0.18854, trip_amp_s = 0, trip_x0 = 10, trip_lx = 4, trip_ly = 1,
%   trip_ts = 4, trip_nmodes = 16.
%
% NONDIMENSIONALISATION
%   Lengths  by delta*_in, the displacement thickness of the inflow profile
%   Velocity by U_inf, the free-stream velocity
%   Pressure by rho*U_inf^2 (kinematic, rho = 1); time by delta*_in/U_inf
%   => nu = 1/Re_dstar_in = 1/450, and every stored quantity is in these
%      "outer" units (e.g. uu is <u'u'>/U_inf^2, tauw is tau_w/(rho U_inf^2)).
%   Wall ("plus") units, used in the profile plots below, are local:
%     u_tau(x) = sqrt(tau_w/rho),  y+ = y*u_tau/nu,  U+ = U/u_tau,
%     <u'u'>+ = <u'u'>/u_tau^2.
%   x is measured from the INFLOW plane, not from the leading edge. The
%   Blasius inflow has a virtual leading edge x_v upstream of it:
%     Re_x = U_inf (x + x_v)/nu = Re_x_in + x*Re_dstar_in,
%     Re_x_in = (Re_dstar_in/1.72079)^2 = 68 388  (Blasius delta*).
%
% VARIABLES IN THE FILE (2D: ny x nx, rows = y, columns = x)
%   x [1 x nx], y [ny x 1]         cell centres (the wall is y = 0)
%   U, V, P                        mean velocity components, mean pressure
%   uu, vv, ww, uv                 Reynolds stresses <u'u'> ... <u'v'>
%   k                              turbulent kinetic energy (uu+vv+ww)/2
%   dUdy                           dU/dy (2nd-order, through U(0) = 0)
%   nut_over_nu                    eddy viscosity from the DNS (Boussinesq):
%                                  nu_t = -<u'v'>/(dU/dy + dV/dx), over nu;
%                                  NaN above delta99 where it is undefined
% 1D (functions of x):
%   Ue                             edge velocity U(x, y = Ly)
%   delta99, dstar, theta, H       99 % thickness, displacement, momentum
%                                  thickness, shape factor dstar/theta
%   tauw, utau, cf                 wall shear (2nd-order one-sided gradient),
%                                  friction velocity, c_f = 2 tau_w/Ue^2
%   Re_x, Re_theta, Re_dstar       Ue*(x+x_v)/nu ... (see above)
%   Re_tau                         utau*delta99/nu
%   p_wall, beta_clauser           wall pressure; Clauser parameter
%                                  (dstar/tau_w) dp_e/dx (~0 here: ZPG)
%   urms_max                       max over y of sqrt(<u'u'>) (transition)
% Scalars: nu, Re_dstar_in, Re_x_in, time_window, n_samples,
%   input_ini (the solver input file, verbatim), source.

clear; close all;
set(groot, 'defaultTextInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', ...
          'defaultAxesTickLabelInterpreter', 'latex');
d = load('tbl_dns_smooth.mat');
x = d.x(:)'; y = d.y(:);
ymax = 1.3*max(d.delta99);           % plot window in y

%% 2D: mean velocity and Reynolds stresses versus (x, y)
maps(x, y, d, {'U','V'}, '$x/\delta^*_{in}$', '$y/\delta^*_{in}$', ...
     x, d.delta99, ymax, 'Mean velocity vs (x, y)');
stress = {'uu','vv','ww','k','uv'};
maps(x, y, d, stress, '$x/\delta^*_{in}$', '$y/\delta^*_{in}$', ...
     x, d.delta99, ymax, 'Reynolds stresses vs (x, y)');

%% 2D: the same versus (Re_theta, y)
rt = d.Re_theta(:)';
maps(rt, y, d, {'U','V',stress{:}}, '$Re_\theta$', '$y/\delta^*_{in}$', ...
     rt, d.delta99, ymax, 'Mean velocity and Reynolds stresses vs $(Re_\theta, y)$');

%% 2D: eddy viscosity nu_t/nu versus (x, y) and (x/theta, y/theta)
maps(x, y, d, {'nut_over_nu'}, '$x/\delta^*_{in}$', '$y/\delta^*_{in}$', ...
     x, d.delta99, ymax, '$\nu_t/\nu$ vs $(x, y)$');
th = d.theta(:)';
[X, Y] = meshgrid(x, y);
Xt = X./th; Yt = Y./th;              % theta(x) varies along x: 2D coordinates
maps(Xt, Yt, d, {'nut_over_nu'}, '$x/\theta$', '$y/\theta$', ...
     x./th, d.delta99(:)'./th, 1.3*max(d.delta99(:)'./th), '$\nu_t/\nu$ vs $(x/\theta, y/\theta)$');

%% 1D: integral quantities
figure('Name', 'Integral quantities', 'Position', [50 50 1300 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact');
Rt = d.Re_theta; Rx = d.Re_x;
rtl = linspace(min(Rt), max(Rt), 200); rxl = linspace(min(Rx), max(Rx), 200);

nexttile;                                            % i) c_f vs Re_theta
plot(Rt, d.cf, 'k', rtl, 0.441./rtl, 'b--', rtl, 0.024*rtl.^(-1/4), 'r--', ...
     rtl, 2*(log(rtl)/0.384 + 4.127).^(-2), 'r:');
xlabel('$Re_\theta$'); ylabel('$c_f$'); ylim([0 0.008]);
legend('DNS', 'laminar (Blasius) $0.441/Re_\theta$', '$0.024\,Re_\theta^{-1/4}$', ...
       'Coles-Fernholz ($\kappa=0.384$, $C=4.127$)', 'Location', 'northeast');

nexttile;                                            % ii) Re_theta vs Re_x
plot(Rx, Rt, 'k', rxl, 0.664*sqrt(rxl), 'b--');
xlabel('$Re_x$'); ylabel('$Re_\theta$');
legend('DNS', 'laminar (Blasius) $0.664\,Re_x^{1/2}$', 'Location', 'northwest');

nexttile;                                            % iii) Re_tau vs Re_theta
plot(Rt, d.Re_tau, 'k', rtl, 1.13*rtl.^0.843, 'r--');
xlabel('$Re_\theta$'); ylabel('$Re_\tau = u_\tau\delta_{99}/\nu$');
legend('DNS', '$1.13\,Re_\theta^{0.843}$ (Schlatter \& Orlu 2010)', 'Location', 'northwest');

nexttile;                                            % c_f vs Re_x: transition location
plot(Rx, d.cf, 'k', rxl, 0.664./sqrt(rxl), 'b--', rxl, 0.0592*rxl.^(-1/5), 'r--');
xlabel('$Re_x$'); ylabel('$c_f$'); ylim([0 0.008]);
legend('DNS', 'laminar $0.664\,Re_x^{-1/2}$', 'turbulent $0.0592\,Re_x^{-1/5}$');

nexttile;                                            % shape factor
plot(Rt, d.H, 'k'); yline(2.59, 'b--', 'Blasius'); yline(1.4, 'r--', 'turbulent');
ylim([1.3 2.7]);
xlabel('$Re_\theta$'); ylabel('$H = \delta^*/\theta$');

nexttile;                                            % peak streamwise fluctuation
plot(Rx, d.urms_max, 'k');
xlabel('$Re_x$'); ylabel('$\max_y\, u''_{rms}/U_\infty$');

%% Profiles at a few stations, wall units
targets = [300 500 690 850];                         % Re_theta of the stations
figure('Name', 'Profiles in wall units', 'Position', [80 80 1300 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact');
for n = 1:6, ax(n) = nexttile; hold on; end         %#ok<SAGROW>
for rtn = targets
    [~, i] = min(abs(d.Re_theta - rtn));
    ut = d.utau(i); yp = y*ut/d.nu;
    lab = sprintf('$Re_\\theta = %.0f$', d.Re_theta(i));
    semilogx(ax(1), yp, d.U(:,i)/ut, 'DisplayName', lab);
    semilogx(ax(2), yp, sqrt(d.uu(:,i))/ut, 'DisplayName', lab);
    semilogx(ax(3), yp, -d.uv(:,i)/ut^2, 'DisplayName', lab);
    semilogx(ax(4), yp, sqrt(d.vv(:,i))/ut, 'DisplayName', lab);
    semilogx(ax(5), yp, sqrt(d.ww(:,i))/ut, 'DisplayName', lab);
    plot(ax(6), y/d.delta99(i), d.nut_over_nu(:,i), 'DisplayName', lab);
end
yv = logspace(-1, log10(15), 50); yl = logspace(1, 3, 50);    % viscous sublayer, log law
plot(ax(1), yv, yv, 'k:', yl, log(yl)/0.41 + 5.0, 'k--', 'HandleVisibility', 'off');
ylim(ax(1), [0 30]);
names = {'$U^+$', '$u''^+_{rms}$', '$-\langle u''v''\rangle^+$', ...
         '$v''^+_{rms}$', '$w''^+_{rms}$', '$\nu_t/\nu$'};
for n = 1:6
    set(ax(n), 'XScale', 'log'); xlabel(ax(n), '$y^+$'); ylabel(ax(n), names{n});
    box(ax(n), 'on');
end
set(ax(6), 'XScale', 'lin'); xlabel(ax(6), '$y/\delta_{99}$'); xlim(ax(6), [0 1]);
legend(ax(1), 'Location', 'northwest');

%% ------------------------------------------------------------------------
function maps(X, Y, d, names, xl, yl, xe, de, ymax, ttl)
% Filled contour maps of the fields NAMES of struct D on the grid (X, Y);
% the delta99 line (xe, de) is drawn on top.
figure('Name', regexprep(ttl, '[$\\]', ''), 'Position', [30 30 1300 300*ceil(numel(names)/2)]);
tiledlayout(ceil(numel(names)/2), min(2, numel(names)), 'TileSpacing', 'compact');
for n = 1:numel(names)
    nexttile;
    contourf(X, Y, d.(names{n}), 40, 'LineColor', 'none'); hold on;
    plot(xe, de, 'w-', 'LineWidth', 1);
    ylim([0 ymax]); colorbar;
    xlabel(xl, 'Interpreter', 'latex'); ylabel(yl, 'Interpreter', 'latex');
    title(names{n}, 'Interpreter', 'none');
end
sgtitle(ttl, 'Interpreter', 'latex');
end
