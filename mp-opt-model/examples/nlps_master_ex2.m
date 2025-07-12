function nlps_master_ex2(alg)
% nlps_master_ex2 - Example of constrained nonlinear optimization.
% ::
%
%   nlps_master_ex2()
%   nlps_master_ex2('MIPS')
%   nlps_master_ex2('FMINCON')
%   nlps_master_ex2('IPOPT')
%   nlps_master_ex2('KNITRO')
%
% Example of solving the following problem using nlps_master.
%
% .. math:: \min_{x_1,x_2,x_3} -x_1 x_2 - x_2 x_3
%
% subject to
%
% .. math:: x_1^2 - x_2^2 + x_3^2 - 2 \le 0
% .. math:: x_1^2 + x_2^2 + x_3^2 - 10 \le 0

if nargin < 1
    alg = 'DEFAULT';
end
problem = struct( ...
    'f_fcn',    @(x)f2(x), ...
    'gh_fcn',   @(x)gh2(x), ...
    'hess_fcn', @(x, lam, cost_mult)hess2(x, lam, cost_mult), ...
    'x0',       [1; 1; 0], ...
    'opt',      struct('verbose', 2, 'alg', alg) ...
);
[x, f, exitflag, output, lambda] = nlps_master(problem);
fprintf('\nf = %g   exitflag = %d\n', f, exitflag);
fprintf('\nx = \n');
fprintf('   %g\n', x);
fprintf('\nlambda.ineqnonlin =\n');
fprintf('   %g\n', lambda.ineqnonlin);

function [f, df, d2f] = f2(x)
f = -x(1)*x(2) - x(2)*x(3);
if nargout > 1          %% gradient is required
    df = -[x(2); x(1)+x(3); x(2)];
    if nargout > 2      %% Hessian is required
        d2f = -[0 1 0; 1 0 1; 0 1 0];   %% actually not used since
    end                                 %% 'hess_fcn' is provided
end

function [h, g, dh, dg] = gh2(x)
h = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
dh = 2 * [x(1) x(1); -x(2) x(2); x(3) x(3)];
g = []; dg = [];

function Lxx = hess2(x, lam, cost_mult)
if nargin < 3, cost_mult = 1; end   %% allows to be used with 'fmincon'
mu = lam.ineqnonlin;
Lxx = cost_mult * [0 -1 0; -1 0 -1; 0 -1 0] + ...
        [2*[1 1]*mu 0 0; 0 2*[-1 1]*mu 0; 0 0 2*[1 1]*mu];
