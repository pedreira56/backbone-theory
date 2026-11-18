function out = E5_coupling(regimes, cfg)
%E5_COUPLING  Experiment E5: the limiting coupling with several blocks (T4).
%
%   OUT = E5_COUPLING(REGIMES,CFG) measures the Wasserstein discrepancy
%   between the law the reconstruction induces and the product coupling, and
%   checks the two things Theorem T4 asserts about it: the rate, and the way
%   the constant assembles from the per-block anchoring constants.
%
%   WHY THE LEFT-HAND SIDE IS EXACT
%   Theorem T4 bounds the discrepancy by tensorisation.  For product measures
%   that step is an equality, not merely an inequality: any coupling of the
%   products projects to couplings of the components, so the quadratic
%   Wasserstein distance between two products equals the sum of the distances
%   between corresponding components.  With homoscedastic blocks the
%   conditional covariance does not depend on the backbone, the Bures term
%   vanishes, and the per-block distance reduces to the squared distance
%   between conditional means.  The quantity on the left of T4 is therefore
%   computed exactly, with no simulation of residuals and no estimation of a
%   Wasserstein distance -- which matters, because the empirical estimator of
%   W_2 in dimension D converges at n^{-1/D} and would mask the very rate
%   under test.
%
%   THE TWO PREDICTIONS
%     rate      the discrepancy is O(n^{-2 beta / d_B}), by T1 applied to the
%               moment of order 2 beta.
%     constant  it is controlled by sum over donor blocks of Lambda_l^2.  The
%               prediction is not that the error grows with the number of
%               blocks, which would be unsurprising, but that it depends on
%               the per-block constants ONLY through that sum: two different
%               arrangements with the same sum must give the same error.
%
%   CFG fields (optional): nGrid, nRecv, repMin, repMax, fitFrom, seed,
%   verbose.

if nargin < 2, cfg = struct(); end
dflt = struct('nGrid',round(logspace(2,5,10)),'nRecv',2000, ...
              'repMin',8,'repMax',40,'fitFrom',0.5, ...
              'seed',86420,'verbose',true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

nGrid = cfg.nGrid(:).';
nN    = numel(nGrid);
nReg  = numel(regimes);

W2sum  = nan(nReg,nN);   % exact discrepancy, summed over donor blocks
R2b    = nan(nReg,nN);   % E[R^{2 beta}], the factor T1 controls
lam2    = nan(nReg,1);   % sum over donor blocks of Lambda_l^2
Mv = nan(nReg,1); betaV = nan(nReg,1); dBV = nan(nReg,1);

for g = 1:nReg
    P = bbParams(regimes(g).params{:});
    assert(P.hetero == 0, 'E5_coupling:needsHomoscedastic', ...
        ['Regime "%s" must set hetero to zero: the exact per-block ' ...
         'Wasserstein distance assumes a backbone-free conditional ' ...
         'covariance.'], regimes(g).name);
    Mv(g) = P.M; betaV(g) = P.beta; dBV(g) = P.d;

    % Lambda_l^2 for block l: the squared Lipschitz constant of r_l, which
    % for the smooth mean is the sum over the block's attributes of
    % (Lmean * sr_a)^2.  Block 1 is the receiver's block and is excluded.
    lam2(g) = 0;
    for m = 2:P.M
        idx = (P.blockId == m);
        lam2(g) = lam2(g) + sum((P.Lmean * P.sr(idx)).^2);
    end

    if cfg.verbose
        fprintf('[E5] regime %d/%d: %-26s M=%d  sum Lambda^2 = %8.3f\n', ...
            g,nReg,regimes(g).name,P.M,lam2(g));
    end

    for j = 1:nN
        n   = nGrid(j);
        rep = min(cfg.repMax, max(cfg.repMin, ceil(cfg.repMax*sqrt(1e3/n))));
        accW = zeros(rep,1); accR = zeros(rep,1);
        for t = 1:rep
            rs = RandStream('mt19937ar','Seed', cfg.seed + 1e6*g + 1e3*j + t);
            Sr = bbSample(cfg.nRecv,P,rs);       % receivers
            tot = zeros(cfg.nRecv,1);
            rsum = zeros(cfg.nRecv,1);
            for m = 2:P.M
                Sd = bbSample(n,P,rs);           % an independent donor block
                [Rv,J] = bbNNDist(Sr.XB, Sd.XB, P);
                dm  = Sd.r{m}(J,:) - Sr.r{m};
                tot = tot + sum(dm.^2,2);
                rsum = rsum + Rv.^(2*P.beta);
            end
            accW(t) = mean(tot);
            accR(t) = mean(rsum)/(P.M-1);        % average over donor blocks
        end
        W2sum(g,j) = mean(accW);
        R2b(g,j)   = mean(accR);
    end
end

% ---- rate ------------------------------------------------------------
jFit = max(1,ceil((1-cfg.fitFrom)*nN)) : nN;
X    = [ones(numel(jFit),1), log(nGrid(jFit)).'];
slope = nan(nReg,1);
for g = 1:nReg
    bh = X\log(reshape(W2sum(g,jFit),[],1));
    slope(g) = bh(2);
end

% ---- the two diagnostics of the constant -----------------------------
% ratio    the discrepancy per unit of sum Lambda^2.  T4's constant structure
%          predicts this is the same across regimes that differ only in how
%          the anchoring constants are distributed among blocks.
% tight    the discrepancy divided by the majorant of T4.  Must not exceed
%          one; how far below one it falls measures the slack introduced by
%          Cond_A*, not by the tensorisation step, which is exact here.
ratio = W2sum(:,end) ./ lam2;
tight = W2sum(:,end) ./ (lam2 .* R2b(:,end));

tbl = table(string({regimes.name}).', Mv, betaV, dBV, lam2, ...
            -2*betaV./dBV, slope, ratio, tight, ...
    'VariableNames',{'regime','M','beta','d_B','sumLambda2', ...
                     'theory','slope','discr_per_Lambda2','ratio_to_bound'});

out = struct('tbl',tbl,'nGrid',nGrid,'W2sum',W2sum,'R2b',R2b, ...
             'lam2',lam2,'slope',slope,'ratio',ratio,'tight',tight, ...
             'M',Mv,'beta',betaV,'dB',dBV,'cfg',cfg,'regimes',regimes, ...
             'jFit',jFit);
end
