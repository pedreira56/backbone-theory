function out = E6_enumeration(cfg)
%E6_ENUMERATION  Experiment E6: exhaustive enumeration of backbone candidates.
%
%   OUT = E6_ENUMERATION(CFG) evaluates the predictive selection criterion of
%   Section 8 at EVERY admissible backbone, by enumeration rather than by
%   search.  No relaxation, no gradient, no heuristic: with a small universe
%   of eligible attributes the optimum is available exactly, so the criterion
%   can be examined rather than optimised.
%
%   WHAT IS CHECKED
%   (1) That the criterion selects the backbone that anchors the exclusive
%       attributes.  The generator has a known latent dimension and a known
%       set of attributes that read it cleanly; if the formulation is sound,
%       the minimiser at the corresponding cardinality is that set.
%   (2) That backbones of equal cardinality are not interchangeable.  Section
%       8 asserts this; the distribution of the criterion over all subsets of
%       a given size measures how far apart they are.
%   (3) The behaviour of the criterion as the backbone grows, and in
%       particular whether the degeneracy identified in Section 8.3 -- every
%       added coordinate removes a term from the sum, which favours large
%       backbones unconditionally -- is visible.
%   (4) That per-attribute costs remove that degeneracy and yield a
%       well-posed problem without an artificial cardinality constraint.
%
%   THE GENERATOR
%   Jointly Gaussian, which makes the two population quantities entering the
%   criterion available in closed form rather than estimated: the residual
%   conditional variance is a Schur complement, and the anchoring constant is
%   the norm of the regression coefficient vector, since the conditional mean
%   is linear.  Only the nearest-neighbour moment is simulated.  All
%   attributes are standardised to unit variance, so that Euclidean distance
%   in backbone space does not depend on an arbitrary choice of units.
%
%   CFG fields (optional): d, p, bMax, n, nRecv, repDelta, costs, seed,
%   verbose.

if nargin < 1, cfg = struct(); end
dflt = struct('d',3,'p',14,'bMax',6,'n',3000,'nRecv',400,'repDelta',4, ...
              'costs',[0 0.10 0.25],'seed',24601,'verbose',true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end
d = cfg.d; p = cfg.p;

% ---- population -------------------------------------------------------
% Attribute j is a_j' Z + tau_j e_j with Z standard normal in d dimensions.
%   1..d        clean readings of the latent coordinates
%   d+1..2d     noisy duplicates of the same coordinates
%   2d+1..p     targets, partially predictable from the latent vector
% The intended optimum at cardinality d is therefore the first d attributes:
% they pin the latent vector down with the least added noise.
rs = RandStream('mt19937ar','Seed',cfg.seed);
A = zeros(p,d); tau = zeros(p,1);
A(1:d,:)       = eye(d);      tau(1:d)       = 0.25;
A(d+1:2*d,:)   = eye(d);      tau(d+1:2*d)   = 0.90;
V = randn(rs,p-2*d,d); V = V./vecnorm(V,2,2);
A(2*d+1:p,:)   = V;           tau(2*d+1:p)   = 0.60;

S  = A*A.' + diag(tau.^2);
sc = sqrt(diag(S));
S  = S ./ (sc*sc.');                      % unit variances
out.Sigma = S;  out.trueB = 1:d;

% ---- enumerate --------------------------------------------------------
subsets = {}; card = [];
for b = 0:cfg.bMax
    C = nchoosek(1:p,b);
    if b == 0, C = zeros(1,0); end
    for i = 1:size(C,1)
        subsets{end+1} = C(i,:);                                      %#ok<AGROW>
        card(end+1) = b;                                              %#ok<AGROW>
    end
end
nS = numel(subsets);
if cfg.verbose
    fprintf('[E6] enumerating %d backbones with cardinality up to %d\n',nS,cfg.bMax);
end

Ubar = nan(nS,1);    % the criterion
Sres = nan(nS,1);    % summed residual conditional variance
Uloc = nan(nS,1);    % local part, per remaining attribute
Del2 = nan(nS,1);    % E[R^2] in backbone space

for i = 1:nS
    B    = subsets{i};
    rest = setdiff(1:p,B);
    if isempty(rest), Ubar(i)=0; Sres(i)=0; Uloc(i)=0; Del2(i)=0; continue, end
    if isempty(B)
        Sres(i) = sum(diag(S(rest,rest)));
        Uloc(i) = 0; Del2(i) = 0; Ubar(i) = Sres(i);
        continue
    end
    SBB = S(B,B);
    % In high dimension the nearest-neighbour search dominates the cost of
    % the whole enumeration, while the criterion there is driven by the
    % vanishing number of remaining attributes rather than by the geometry.
    % Trimming the receiver sample above a moderate cardinality keeps the
    % enumeration to the full universe affordable without affecting the
    % comparisons that matter, which occur at small b.
    nq = cfg.nRecv; if numel(B) > 8, nq = max(120, round(cfg.nRecv/3)); end
    Del2(i) = nnMoment(SBB, cfg.n, nq, cfg.repDelta, cfg.seed + i);
    coef = SBB \ S(B,rest);                       % |B| x |rest|
    s2   = diag(S(rest,rest)).' - sum(S(B,rest).*coef,1);
    La2  = sum(coef.^2,1);
    Sres(i) = sum(s2);
    Uloc(i) = sum(La2)*Del2(i)/numel(rest);
    Ubar(i) = Sres(i) + sum(La2)*Del2(i);
    if cfg.verbose && mod(i,500)==0
        fprintf('   %5d / %5d\n',i,nS);
    end
end

% ---- summaries --------------------------------------------------------
% best at each cardinality, with and without per-attribute costs
nb   = cfg.bMax + 1;
Phi  = nan(nb,1); PhiLoc = nan(nb,1); bestB = cell(nb,1);
Utot = nan(nb,numel(cfg.costs)); bestBc = cell(nb,numel(cfg.costs));
for b = 0:cfg.bMax
    idx = find(card == b);
    [Phi(b+1),k] = min(Ubar(idx));
    bestB{b+1}   = subsets{idx(k)};
    PhiLoc(b+1)  = min(Uloc(idx));
    for c = 1:numel(cfg.costs)
        [Utot(b+1,c),k2] = min(Ubar(idx) + cfg.costs(c)*b);
        bestBc{b+1,c}    = subsets{idx(k2)};
    end
end

% ranking at the latent cardinality
idx3 = find(card == d);
[v3,ord] = sort(Ubar(idx3));
rank3 = cell(numel(ord),1);
for i = 1:numel(ord), rank3{i} = subsets{idx3(ord(i))}; end
isTrue = isequal(sort(rank3{1}), 1:d);
% Value and rank of the INTENDED backbone, so the figure can mark it without
% assuming it is the minimiser.  Marking the minimum and labelling it as the
% intended set would be a claim rather than a measurement.
valTrue = NaN; rankTrue = NaN;
for i = 1:numel(rank3)
    if isequal(sort(rank3{i}), 1:d), valTrue = v3(i); rankTrue = i; break, end
end

% dispersion at fixed cardinality
disp_tbl = table((1:cfg.bMax)', 'VariableNames',{'b'});
mn = nan(cfg.bMax,1); md = mn; mx = mn; ct = mn;
for b = 1:cfg.bMax
    v = Ubar(card == b);
    mn(b)=min(v); md(b)=median(v); mx(b)=max(v); ct(b)=numel(v);
end
disp_tbl.count = ct; disp_tbl.min = mn; disp_tbl.median = md;
disp_tbl.max = mx;   disp_tbl.max_over_min = mx./mn;

out.subsets=subsets; out.card=card(:); out.Ubar=Ubar; out.Sres=Sres;
out.Uloc=Uloc; out.Delta2=Del2; out.Phi=Phi; out.PhiLoc=PhiLoc;
out.bestB=bestB; out.Utot=Utot; out.bestBcost=bestBc; out.costs=cfg.costs;
out.rank3=rank3; out.val3=v3; out.trueIsBest=isTrue; out.dispersion=disp_tbl;
out.valTrue=valTrue; out.rankTrue=rankTrue;
out.cfg=cfg;
end

% =======================================================================
function m = nnMoment(SBB, n, nRecv, rep, seed)
%NNMOMENT  E[R^2] for a Gaussian backbone with covariance SBB.
%   Depends on the candidate only through SBB, which is why the criterion is
%   comparable across backbones of the same cardinality.
b = size(SBB,1);
L = chol(SBB + 1e-12*eye(b),'lower');
acc = zeros(rep,1);
for t = 1:rep
    rs = RandStream('mt19937ar','Seed',seed + 977*t);
    Xd = randn(rs,n,b)     * L.';
    Xq = randn(rs,nRecv,b) * L.';
    [~,D] = knnsearch(createns(Xd,'NSMethod','kdtree'), Xq, 'K',1);
    acc(t) = mean(D.^2);
end
m = mean(acc);
end
