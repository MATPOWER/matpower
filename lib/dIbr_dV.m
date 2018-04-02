function [dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart)
%DIBR_DV   Computes partial derivatives of branch currents w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 5th argument.
%
%   [DIF_DVA, DIF_DVM, DIT_DVA, DIT_DVM, IF, IT] = DIBR_DV(BRANCH, YF, YT, V)
%   [DIF_DVA, DIF_DVM, DIT_DVA, DIT_DVM, IF, IT] = DIBR_DV(BRANCH, YF, YT, V, 0)
%
%   Returns four matrices containing partial derivatives of the complex
%   branch currents at "from" and "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle, respectively (for all buses).
%
%   [DIF_DVR, DIF_DVI, DIT_DVR, DIT_DVI, IF, IT] = DIBR_DV(BRANCH, YF, YT, V, 1)
%
%   Returns four matrices containing partial derivatives of the complex
%   branch currents at "from" and "to" ends of each branch w.r.t real and
%   imaginary parts of voltage, respectively (for all buses).
%
%   If YF is a sparse matrix, the partial derivative matrices will be as well.
%   Optionally returns vectors containing the currents themselves. The
%   following explains the expressions used to form the matrices:
%
%   If = Yf * V;
%
%   Polar coordinates:
%     Partials of V, Vf & If w.r.t. voltage angles
%       dV/dVa  = j * diag(V)
%       dVf/dVa = sparse(1:nl, f, j * V(f)) = j * sparse(1:nl, f, V(f))
%       dIf/dVa = Yf * dV/dVa = Yf * j * diag(V)
%
%     Partials of V, Vf & If w.r.t. voltage magnitudes
%       dV/dVm  = diag(V./abs(V))
%       dVf/dVm = sparse(1:nl, f, V(f)./abs(V(f))
%       dIf/dVm = Yf * dV/dVm = Yf * diag(V./abs(V))
%
%   Cartesian coordinates:
%     Partials of V, Vf & If w.r.t. real part of complex voltage
%       dV/dVr  = diag(ones(n,1))
%       dVf/dVr = Cf
%       dIf/dVr = Yf
%     where Cf is the connection matrix for line & from buses
%
%     Partials of V, Vf & If w.r.t. imaginary part of complex voltage
%       dV/dVi  = j * diag(ones(n,1))
%       dVf/dVi = j * Cf
%       dIf/dVi = j * Yf
%
%   Derivations for "to" bus are similar.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = ...
%           dIbr_dV(branch, Yf, Yt, V);
%       [dIf_dVr, dIf_dVi, dIt_dVr, dIt_dVi, If, It] = ...
%           dIbr_dV(branch, Yf, Yt, V, 1);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

%% define
nb = length(V);

if vcart
    dIf_dV1 = Yf;                   %% dIf_dVr
    dIf_dV2 = 1j * Yf;              %% dIf_dVi
    dIt_dV1 = Yt;                   %% dIt_dVr
    dIt_dV2 = 1j * Yt;              %% dIt_dVi
else
    Vnorm = V ./ abs(V);
    if issparse(Yf)             %% sparse version (if Yf is sparse)
        diagV       = sparse(1:nb, 1:nb, V, nb, nb);
        diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
    else                        %% dense version
        diagV       = diag(V);
        diagVnorm   = diag(Vnorm);
    end
    dIf_dV1 = Yf * 1j * diagV;      %% dIf_dVa
    dIf_dV2 = Yf * diagVnorm;       %% dIf_dVm
    dIt_dV1 = Yt * 1j * diagV;      %% dIt_dVa
    dIt_dV2 = Yt * diagVnorm;       %% dIt_dVm
end

%% compute currents
if nargout > 4
    If = Yf * V;
    It = Yt * V;
end
