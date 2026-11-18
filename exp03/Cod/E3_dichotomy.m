function out = E3_dichotomy(regimes, cfg)
%E3_DICHOTOMY  Experiment E3: the pointwise / distributional dichotomy.
%
%   OUT = E3_DICHOTOMY(REGIMES,CFG) evaluates a family of transfer rules
%   under two criteria that Propositions P3-P4 and Theorems T2-T3 predict to
%   rank them in OPPOSITE orders.
%
%   THE RULES
%   Averaging the exclusive vectors of the k nearest donors, for k on a grid,
%   plus the oracle rule that transfers the true conditional mean of the
%   single nearest donor.  Increasing k is the natural way to interpolate:
%   k = 1 transfers one observed vector and keeps its full dispersion, large
%   k averages the dispersion away, and the oracle mean removes it entirely.
%
%   THE TWO CRITERIA
%     individual   E|| X_Q - imputed ||^2, the quadratic risk of Theorem T2.
%     distributional  W_2^2 between the CONDITIONAL LAW of the imputed vector
%                  given the backbones and the true conditional law of the
%                  receiver.  This is the object of Proposition P4 and
%                  Theorem T3.
%
%   WHY NO RESIDUALS ARE DRAWN
%   Both criteria are available in closed form once the backbones are fixed.
%   Conditionally on the backbones, the imputed vector of the k-average is
%   Gaussian with mean equal to the average of the donor conditional means
%   and covariance equal to k^{-2} times the sum of the donor conditional
%   covariances.  Under sdShared every conditional covariance is a scalar
%   multiple of one fixed matrix, so the Wasserstein distance collapses to
%       W_2^2 = || mean difference ||^2 + tr(Sigma_m) (a - b)^2,
%   and the individual risk is tau^2(x_B) + || mean difference ||^2 plus the
%   trace of the imputed covariance.  Everything is therefore exact, and the
%   only sampling error left is that of averaging over receivers.
%
%   CFG fields (optional): nGrid, nRecv, repMin, repMax, kList, fitFrom,
%   seed, nOracle, verbose.

if nargin < 2, cfg = struct(); end
dflt = struct('nGrid',   round(logspace(2,5,10)), ...
              'nRecv',   500, ...
              'repMin',  20, ...
              'repMax',  200, ...
              'kList',   [1 2 3 5 10 20 50 100], ...
              'fitFrom', 0.5, ...
              'seed',    13579, ...
              'nOracle', 2e6, ...
              'verbose', true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

nGrid = cfg.nGrid(:).';
nN    = numel(nGrid);
nK    = numel(cfg.kList);
nReg  = numel(regimes);

riskK = nan(nReg,nN,nK);   % individual risk, k-average rule
w2K   = nan(nReg,nN,nK);   % distributional criterion, k-average rule
riskM = nan(nReg,nN);      % individual risk, oracle mean rule
w2M   = nan(nReg,nN);      % distributional criterion, oracle mean rule
tau2E = nan(nReg,1);       % E[tau_l^2(X_B)], the floor of Proposition P4
betaV = nan(nReg,1);  dBV = nan(nReg,1);

for g = 1:nReg
    P = bbParams(regimes(g).params{:});
    assert(P.sdShared, 'E3_dichotomy:needsShared', ...
        ['Regime "%s" must set sdShared true: the closed-form Wasserstein ' ...
         'distance requires proportional conditional covariances.'], ...
         regimes(g).name);
    betaV(g) = P.beta;  dBV(g) = P.d;

    idx1  = find(P.blockId == 1);
    trSig = trace(P.SigmaE(idx1,idx1));     % tr(Sigma_m), the block trace

    if cfg.verbose
        fprintf('[E3] regime %d/%d: %-24s beta=%.2f d=%d\n', ...
            g,nReg,regimes(g).name,P.beta,P.d);
    end

    % Floor E[tau^2(X_B)] = sigma0^2 tr(Sigma_m) E[g(X_B)^2].
    rsO = RandStream('mt19937ar','Seed', cfg.seed + 5555 + g);
    So  = bbSample(cfg.nOracle, P, rsO);
    gO  = So.sd{1}(:,1)/P.sigma0;
    tau2E(g) = P.sigma0^2 * trSig * mean(gO.^2);
    clear So gO

    for j = 1:nN
        n   = nGrid(j);
        rep = min(cfg.repMax, max(cfg.repMin, ceil(cfg.repMax*sqrt(1e3/n))));
        Kmx = min(max(cfg.kList), n);

        accR = zeros(rep,nK); accW = zeros(rep,nK);
        accRM = zeros(rep,1); accWM = zeros(rep,1);

        for t = 1:rep
            rs = RandStream('mt19937ar','Seed', cfg.seed + 1e6*g + 1e3*j + t);
            Sd = bbSample(n,        P, rs);
            Sr = bbSample(cfg.nRecv,P, rs);

            ns = createns(Sd.XB,'NSMethod','kdtree');
            J  = knnsearch(ns, Sr.XB, 'K', Kmx);   % nRecv x Kmx

            rRec = Sr.r{1};                        % nRecv x pm
            gRec = Sr.sd{1}(:,1)/P.sigma0;         % nRecv x 1
            tauR = P.sigma0^2 * trSig * gRec.^2;   % tau^2 at the receiver
            rDon = Sd.r{1};                        % n x pm
            gDon = Sd.sd{1}(:,1)/P.sigma0;         % n x 1

            for kk = 1:nK
                k  = min(cfg.kList(kk), Kmx);
                Jk = J(:,1:k);

                % mean of the k donor conditional means
                mBar = zeros(size(rRec));
                for c = 1:size(rRec,2)
                    col = rDon(:,c);
                    mBar(:,c) = mean(reshape(col(Jk),[],k),2);
                end
                dm2 = sum((mBar - rRec).^2, 2);

                % conditional covariance of the k-average is
                % (sigma0^2 / k^2) * sum_i g(z_i)^2 * Sigma_m
                g2s = sum(reshape(gDon(Jk),[],k).^2, 2);
                aV  = P.sigma0 * sqrt(g2s) / k;    % scale of the imputed law
                bV  = P.sigma0 * gRec;             % scale of the target law

                accR(t,kk) = mean(tauR + dm2 + trSig * aV.^2);
                accW(t,kk) = mean(dm2 + trSig * (aV - bV).^2);
            end

            % oracle mean rule: a point mass at the nearest donor's mean
            mOne = rDon(J(:,1),:);
            dm2  = sum((mOne - rRec).^2, 2);
            accRM(t) = mean(tauR + dm2);
            accWM(t) = mean(tauR + dm2);
        end

        riskK(g,j,:) = mean(accR,1);   w2K(g,j,:) = mean(accW,1);
        riskM(g,j)   = mean(accRM);    w2M(g,j)   = mean(accWM);
    end
end

% ---- rate of the distributional criterion for the value rule (T3) ----
jFit  = max(1,ceil((1-cfg.fitFrom)*nN)) : nN;
X     = [ones(numel(jFit),1), log(nGrid(jFit)).'];
slopeW = nan(nReg,1);
for g = 1:nReg
    bh = X\log(reshape(w2K(g,jFit,1),[],1));
    slopeW(g) = bh(2);
end

% ---- consistency with T2: the value/mean risk ratio must tend to 2 ---
ratio2 = riskK(:,end,1) ./ riskM(:,end);
% ---- consistency with P4: the mean rule's W2^2 must tend to the floor -
floorRatio = w2M(:,end) ./ tau2E;

tbl = table(string({regimes.name}).', betaV, dBV, ...
            -2*betaV./dBV, slopeW, tau2E, ratio2(:), floorRatio(:), ...
    'VariableNames',{'regime','beta','d_B','theoryW','slopeW', ...
                     'floor_tau2','riskRatio_val_med','w2med_over_floor'});

out = struct('tbl',tbl,'nGrid',nGrid,'kList',cfg.kList, ...
             'riskK',riskK,'w2K',w2K,'riskM',riskM,'w2M',w2M, ...
             'tau2E',tau2E,'slopeW',slopeW,'beta',betaV,'dB',dBV, ...
             'cfg',cfg,'regimes',regimes,'jFit',jFit);
end
