function P = bbParams(varargin)
%BBPARAMS  Parameter struct for the block-design population generator.
%
%   P = BBPARAMS('name',value,...) returns a validated struct describing a
%   population on R^p together with a block measurement design:
%
%       Z    ~ Unif[0,1]^d                       (latent coordinates)
%       X_B  = phi(Z)                in R^b      (backbone / common block)
%       X_Qm = r_m(Z) + s_m(Z).*E_m  in R^pm     (exclusive block m)
%
%   The map phi is bi-Lipschitz onto its image, so the geometric dimension
%   of supp(P_B) is exactly d.  This is what makes Cond_G hold with the
%   intended exponent d_B = d, and it is the property a referee will check.
%
%   OPTIONS
%     d          intrinsic dimension of the backbone support        [2]
%     b          ambient dimension of the backbone                  [2]
%     embedding  'identity' | 'linear' | 'graph'                    ['identity']
%     geometry   'cube' | 'torus'  (torus needs 'identity')         ['cube']
%     kGraph     extra coordinates used by the 'graph' embedding    [2]
%     Lg         Lipschitz constant of the graph map g              [1.0]
%     M          number of measurement blocks                       [2]
%     pm         number of exclusive attributes per block           [3]
%     meanType   'abspow' (exact Holder-beta) | 'smooth'            ['abspow']
%     beta       Holder exponent of the conditional mean, in (0,1]  [1.0]
%     Lmean      scale of the anchoring constant                    [1.0]
%     sigma0     baseline conditional standard deviation            [0.5]
%     hetero     heteroscedasticity amplitude in [0,1)              [0.0]
%     rhoWithin  within-block correlation of the innovations        [0.3]
%     theta      cross-block correlation; theta = 0 <=> Cond_CI     [0.0]
%     seed       RNG seed for the fixed structural constants        [0]
%
%   The structural constants (rotation Q, directions w, offsets t) are drawn
%   once from a dedicated stream so that a given seed always yields the same
%   population, independently of any sampling done later.
%
%   See also BBSAMPLE, BBEMBED, BBNNDIST.

opts = struct( ...
    'd',         2,          ...
    'b',         2,          ...
    'embedding', 'identity', ...
    'geometry',  'cube',     ...
    'kGraph',    2,          ...
    'Lg',        1.0,        ...
    'M',         2,          ...
    'pm',        3,          ...
    'meanType',  'abspow',   ...
    'beta',      1.0,        ...
    'Lmean',     1.0,        ...
    'sigma0',    0.5,        ...
    'hetero',    0.0,        ...
    'rhoWithin', 0.3,        ...
    'theta',     0.0,        ...
    'seed',      0           );

if mod(numel(varargin),2) ~= 0
    error('bbParams:badArgs','Arguments must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    nm = varargin{k};
    if ~ischar(nm) && ~isstring(nm)
        error('bbParams:badArgs','Option names must be character vectors.');
    end
    nm = char(nm);
    if ~isfield(opts,nm)
        error('bbParams:unknownOption','Unknown option "%s".',nm);
    end
    opts.(nm) = varargin{k+1};
end
P = opts;

% ---- validation -------------------------------------------------------
mustBePos = @(x,nm) assert(isscalar(x) && isnumeric(x) && x > 0, ...
    'bbParams:badValue','%s must be a positive scalar.',nm);
mustBePos(P.d,'d');  mustBePos(P.b,'b');  mustBePos(P.M,'M');
assert(P.d == round(P.d) && P.b == round(P.b), ...
    'bbParams:badValue','d and b must be integers.');
assert(P.beta > 0 && P.beta <= 1, ...
    'bbParams:badValue','beta must lie in (0,1].');
assert(P.hetero >= 0 && P.hetero < 1, ...
    'bbParams:badValue','hetero must lie in [0,1).');

switch P.embedding
    case 'identity'
        assert(P.b == P.d, 'bbParams:badValue', ...
            'embedding ''identity'' requires b == d (got b=%d, d=%d).',P.b,P.d);
        P.kGraph = 0;
    case 'linear'
        assert(P.b >= P.d, 'bbParams:badValue', ...
            'embedding ''linear'' requires b >= d.');
        P.kGraph = 0;
    case 'graph'
        assert(P.b >= P.d + P.kGraph, 'bbParams:badValue', ...
            'embedding ''graph'' requires b >= d + kGraph (got b=%d, d+k=%d).', ...
            P.b, P.d + P.kGraph);
    otherwise
        error('bbParams:badValue','Unknown embedding "%s".',P.embedding);
end

if strcmp(P.geometry,'torus')
    assert(strcmp(P.embedding,'identity'), 'bbParams:badValue', ...
        'geometry ''torus'' is only defined for the identity embedding.');
elseif ~strcmp(P.geometry,'cube')
    error('bbParams:badValue','Unknown geometry "%s".',P.geometry);
end

% ---- structural constants (fixed by seed) -----------------------------
rs = RandStream('mt19937ar','Seed',P.seed);

% Orthonormal frame used by the embedding.  Orthonormal columns make phi an
% isometry on the (d+k)-dimensional coordinate block, which is what pins the
% bi-Lipschitz constants down to [1, sqrt(1+Lg^2)].
switch P.embedding
    case 'identity'
        P.Q = eye(P.d);
    case 'linear'
        [P.Q,~] = qr(randn(rs,P.b,P.d),0);          % b x d
    case 'graph'
        [P.Q,~] = qr(randn(rs,P.b,P.d+P.kGraph),0); % b x (d+k)
end

% Directions of the graph map g: R^d -> R^k, each of unit norm.  Setting
% g_j(z) = (Lg/(2*pi*sqrt(k))) * sin(2*pi*w_j'z) gives Lip(g_j) <= Lg/sqrt(k)
% and hence Lip(g) <= Lg.
if P.kGraph > 0
    W = randn(rs,P.kGraph,P.d);
    P.Wg = W ./ vecnorm(W,2,2);
else
    P.Wg = zeros(0,P.d);
end

% Conditional-mean structure: one direction, offset and scale per exclusive
% attribute.  With meanType 'abspow',
%     r(z) = Lmean * s * |w'z - t|^beta,   ||w|| = 1,
% which is exactly Holder-beta in z with constant Lmean*|s|.  Because the
% lower bi-Lipschitz constant of phi is 1, the Holder constant of r viewed as
% a function of x_B is bounded by the same value.
q = P.M * P.pm;
W = randn(rs,q,P.d);
P.Wr    = W ./ vecnorm(W,2,2);
P.tr    = 0.5 + 0.15*randn(rs,q,1);
P.sr    = 1 + 0.25*randn(rs,q,1);
P.Wsig  = randn(rs,q,P.d);
P.Wsig  = P.Wsig ./ vecnorm(P.Wsig,2,2);

% ---- innovation covariance -------------------------------------------
% Sigma(i,j) = 1 if i == j; rhoWithin if same block; theta otherwise.
blockId = repelem((1:P.M)',P.pm,1);
same    = blockId == blockId';
Sig     = P.theta*ones(q) + (P.rhoWithin - P.theta)*same;
Sig(1:q+1:end) = 1;
[Lc,flag] = chol(Sig,'lower');
if flag ~= 0
    error('bbParams:notPSD', ...
        ['The innovation covariance is not positive definite for ' ...
         'rhoWithin=%.3f, theta=%.3f, M=%d, pm=%d. Shrink theta.'], ...
        P.rhoWithin,P.theta,P.M,P.pm);
end
P.SigmaE  = Sig;
P.cholE   = Lc;
P.blockId = blockId;
P.q       = q;
P.p       = P.b + q;
end
