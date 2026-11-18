function [XB, info] = bbEmbed(Z, P)
%BBEMBED  Bi-Lipschitz embedding of the latent cube into backbone space.
%
%   XB = BBEMBED(Z,P) maps the n-by-d latent matrix Z to an n-by-b matrix of
%   backbone coordinates.
%
%   [XB,INFO] = BBEMBED(...) also returns INFO.lipLo and INFO.lipHi, the
%   theoretical bi-Lipschitz constants
%
%       lipLo * ||z - z'||  <=  ||phi(z) - phi(z')||  <=  lipHi * ||z - z'||.
%
%   Because Unif[0,1]^d is Ahlfors d-regular and phi is bi-Lipschitz onto its
%   image, the pushforward P_B satisfies Cond_G with d_B = d exactly.  The
%   constant c_B is degraded by at most lipHi^{-d}; the exponent is not.
%   This is the identity that experiment E1 is designed to confirm.
%
%   See also BBPARAMS, BBSAMPLE.

n = size(Z,1);
if size(Z,2) ~= P.d
    error('bbEmbed:badInput','Z must have %d columns.',P.d);
end

switch P.embedding
    case 'identity'
        XB = Z;
        lipLo = 1; lipHi = 1;

    case 'linear'
        XB = Z * P.Q.';                     % Q has orthonormal columns
        lipLo = 1; lipHi = 1;

    case 'graph'
        G  = bbGraphMap(Z,P);               % n x k, Lip(G) <= P.Lg
        XB = [Z, G] * P.Q.';
        lipLo = 1; lipHi = sqrt(1 + P.Lg^2);

    otherwise
        error('bbEmbed:badEmbedding','Unknown embedding "%s".',P.embedding);
end

if nargout > 1
    info = struct('lipLo',lipLo,'lipHi',lipHi,'n',n);
end
end

% -------------------------------------------------------------------------
function G = bbGraphMap(Z,P)
%BBGRAPHMAP  Lipschitz map g:R^d -> R^k with Lip(g) <= P.Lg.
%
%   g_j(z) = (Lg / (2*pi*sqrt(k))) * sin(2*pi * w_j' z),   ||w_j|| = 1.
%
%   Each component has gradient of norm at most Lg/sqrt(k), so
%   ||g(z)-g(z')||^2 <= sum_j (Lg^2/k) ||z-z'||^2 = Lg^2 ||z-z'||^2.
k = P.kGraph;
if k == 0
    G = zeros(size(Z,1),0);
    return
end
amp = P.Lg / (2*pi*sqrt(k));
G   = amp * sin(2*pi * (Z * P.Wg.'));
end
