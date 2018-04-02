function [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart)
%DSBR_DV   Computes partial derivatives of branch power flows w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 5th argument.
%
%   [DSF_DVA, DSF_DVM, DST_DVA, DST_DVM, SF, ST] = DSBR_DV(BRANCH, YF, YT, V)
%   [DSF_DVA, DSF_DVM, DST_DVA, DST_DVM, SF, ST] = DSBR_DV(BRANCH, YF, YT, V, 0)
%
%   Returns four matrices containing partial derivatives of the complex
%   branch power flows at "from" and "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle, respectively (for all buses).
%
%   [DSF_DVR, DSF_DVI, DST_DVR, DST_DVI, SF, ST] = DSBR_DV(BRANCH, YF, YT, V, 1)
%
%   Returns four matrices containing partial derivatives of the complex
%   branch power flows at "from" and "to" ends of each branch w.r.t real and
%   imaginary parts of voltage, respectively (for all buses).
%
%   If YF is a sparse matrix, the partial derivative matrices will be as well.
%   Optionally returns vectors containing the power flows themselves. The
%   following explains the expressions used to form the matrices:
%
%   If = Yf * V;
%   Sf = diag(Vf) * conj(If) = diag(conj(If)) * Vf
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
%     Partials of Sf w.r.t. voltage angles
%       dSf/dVa = diag(Vf) * conj(dIf/dVa)
%                       + diag(conj(If)) * dVf/dVa
%               = diag(Vf) * conj(Yf * j * diag(V))
%                       + conj(diag(If)) * j * sparse(1:nl, f, V(f))
%               = -j * diag(Vf) * conj(Yf * diag(V))
%                       + j * conj(diag(If)) * sparse(1:nl, f, V(f))
%               = j * (conj(diag(If)) * sparse(1:nl, f, V(f))
%                       - diag(Vf) * conj(Yf * diag(V)))
%
%     Partials of Sf w.r.t. voltage magnitudes
%       dSf/dVm = diag(Vf) * conj(dIf/dVm)
%                       + diag(conj(If)) * dVf/dVm
%               = diag(Vf) * conj(Yf * diag(V./abs(V)))
%                       + conj(diag(If)) * sparse(1:nl, f, V(f)./abs(V(f)))
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
%     Partials of Sf w.r.t. real part of complex voltage
%       dSf/dVr = conj(diag(If)) * Cf + diag(Vf) * conj(Yf)
%
%     Partials of Sf w.r.t. imaginary part of complex voltage
%       dSf/dVi = j * (conj(diag(If)) * Cf - diag(Vf) * conj(Yf))
%
%   Derivations for "to" bus are similar.
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = ...
%           dSbr_dV(branch, Yf, Yt, V);
%       [dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = ...
%           dSbr_dV(branch, Yf, Yt, V, 1);
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
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default input args
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

%% define
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);

%% compute intermediate values
Yfc = conj(Yf);
Ytc = conj(Yt);
Vc = conj(V);
Ifc = Yfc * Vc;     %% conjugate of "from" current
Itc = Ytc * Vc;     %% conjugate of "to" current

if issparse(Yf)             %% sparse version (if Yf is sparse)
    diagVf  = sparse(1:nl, 1:nl, V(f), nl, nl);
    diagVt  = sparse(1:nl, 1:nl, V(t), nl, nl);
    diagIfc = sparse(1:nl, 1:nl, Ifc, nl, nl);
    diagItc = sparse(1:nl, 1:nl, Itc, nl, nl);
    if ~vcart
        Vnorm       = V ./ abs(V);
        diagVc      = sparse(1:nb, 1:nb, Vc, nb, nb);
        diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
        CVf  = sparse(1:nl, f, V(f), nl, nb);
        CVnf = sparse(1:nl, f, Vnorm(f), nl, nb);
        CVt  = sparse(1:nl, t, V(t), nl, nb);
        CVnt = sparse(1:nl, t, Vnorm(t), nl, nb);
    end
else                        %% dense version
    diagVf  = diag(V(f));
    diagVt  = diag(V(t));
    diagIfc = diag(Ifc);
    diagItc = diag(Itc);
    if ~vcart
        Vnorm       = V ./ abs(V);
        diagVc      = diag(Vc);
        diagVnorm   = diag(Vnorm);
%         CVf        = zeros(nl, nb);    CVf(sub2ind([nl,nb], (1:nl)', f)) = V(f);
%         CVnf       = zeros(nl, nb);    CVnf(sub2ind([nl,nb], (1:nl)', f)) = Vnorm(f);
%         CVt        = zeros(nl, nb);    CVt(sub2ind([nl,nb], (1:nl)', t)) = V(t);
%         CVnt       = zeros(nl, nb);    CVnt(sub2ind([nl,nb], (1:nl)', t)) = Vnorm(t);
        CVf  = full(sparse(1:nl, f, V(f), nl, nb));
        CVnf = full(sparse(1:nl, f, Vnorm(f), nl, nb));
        CVt  = full(sparse(1:nl, t, V(t), nl, nb));
        CVnt = full(sparse(1:nl, t, Vnorm(t), nl, nb));
    end
end
if vcart
    Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
    Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
    Af = diagIfc * Cf;
    Bf = diagVf * Yfc;
    At = diagItc * Ct;
    Bt = diagVt * Ytc;

    dSf_dV1 = Af + Bf;          %% dSf_dVr
    dSf_dV2 = 1j * (Af - Bf);   %% dSf_dVi
    dSt_dV1 = At + Bt;          %% dSt_dVr
    dSt_dV2 = 1j * (At - Bt);   %% dSt_dVi
else
    dSf_dV1 = 1j * (diagIfc * CVf - diagVf * Yfc * diagVc);     %% dSf_dVa
    dSf_dV2 = diagVf * conj(Yf * diagVnorm) + diagIfc * CVnf;   %% dSf_dVm
    dSt_dV1 = 1j * (diagItc * CVt - diagVt * Ytc * diagVc);     %% dSt_dVa
    dSt_dV2 = diagVt * conj(Yt * diagVnorm) + diagItc * CVnt;   %% dSt_dVm
end

if nargout > 4
    Sf = V(f) .* Ifc;
    St = V(t) .* Itc;
end
