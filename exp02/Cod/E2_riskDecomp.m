function out = E2_riskDecomp(regimes, cfg)
%E2_RISKDECOMP  Experiment E2: risk decomposition of local transfer (T2).
%
%   OUT = E2_RISKDECOMP(REGIMES,CFG) measures, for a grid of donor block
%   sizes, the individual quadratic risk of the two transfer rules of
%   Theorem T2 and checks the three statements the theorem makes about them.
%
%   THE TWO RULES
%     mean transfer   the donor's true conditional mean r_a evaluated at the
%                     donor's backbone.  An oracle predictor: it isolates the
%                     neighbourhood approximation error from the statistical
%                     error of estimating r_a.
%     value transfer  the donor's raw observed value.
%
%   WHAT IS CHECKED
%   (1) The two conditional identities.  Given the receiver backbone and the
%       donor sample, T2 states the conditional risks exactly:
%           mean :  sigma_a^2(x_B)             + [r_a(z) - r_a(x_B)]^2
%           value:  sigma_a^2(x_B) + s^2_a(z)  + [r_a(z) - r_a(x_B)]^2
%       Every term on the right is a population quantity the generator knows
%       in closed form, so the check compares a noisy empirical average
%       against an EXACT reference rather than against another estimate.
%   (2) The difference identity, risk_val - risk_mean = E[sigma_a^2(Z)].
%       This is exact at every n, not asymptotic, and is the sharpest of the
%       three: it is the statement that transferring a value rather than a
%       mean costs precisely one conditional variance.
%   (3) The rate.  The excess of the mean-transfer risk over its floor,
%       E[(r_a(Z) - r_a(x_B))^2], is O(n^{-2*beta/d_B}) by Cond_A and T1.
%
%   CFG fields (optional): nGrid, nRecv, repMin, repMax, fitFrom, seed,
%   nOracle (sample size used for the closed-form floor), verbose.

if nargin < 2, cfg = struct(); end
dflt = struct('nGrid',   round(logspace(2,5,10)), ...
              'nRecv',   500, ...
              'repMin',  20, ...
              'repMax',  200, ...
              'fitFrom', 0.5, ...
              'seed',    24680, ...
              'nOracle', 2e6, ...
              'verbose', true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

nGrid = cfg.nGrid(:).';
nN    = numel(nGrid);
nReg  = numel(regimes);

riskMed = nan(nReg,nN);   % empirical risk, mean transfer
riskVal = nan(nReg,nN);   % empirical risk, value transfer
oracMed = nan(nReg,nN);   % same, from the exact conditional expression
oracVal = nan(nReg,nN);
excess  = nan(nReg,nN);   % E[(r(Z) - r(x_B))^2], computed exactly
sigDon  = nan(nReg,nN);   % E[sigma^2(Z)], computed exactly
floorS2 = nan(nReg,1);    % s_a^2(B) = E[sigma_a^2(X_B)]
relMed  = nan(nReg,nN); seMed  = nan(nReg,nN);   % identity discrepancy + s.e.
relVal  = nan(nReg,nN); seVal  = nan(nReg,nN);
relDiff = nan(nReg,nN); seDiff = nan(nReg,nN);
betaV   = nan(nReg,1);
dBV     = nan(nReg,1);

for g = 1:nReg
    P = bbParams(regimes(g).params{:});
    betaV(g) = P.beta;  dBV(g) = P.d;
    if cfg.verbose
        fprintf('[E2] regime %d/%d: %-26s beta=%.2f d=%d %s\n', ...
            g,nReg,regimes(g).name,P.beta,P.d,P.meanType);
    end

    % Floor s_a^2(B): a population constant, estimated once from a large
    % independent draw so its Monte-Carlo error is negligible against the
    % quantities it is subtracted from.
    rsO = RandStream('mt19937ar','Seed', cfg.seed + 7777 + g);
    So  = bbSample(cfg.nOracle, P, rsO);
    floorS2(g) = mean(So.sd{1}(:,1).^2);
    clear So

    for j = 1:nN
        n   = nGrid(j);
        rep = min(cfg.repMax, max(cfg.repMin, ceil(cfg.repMax*sqrt(1e3/n))));

        % Accumulate PER REPLICATION rather than pooling the receivers.
        % Receivers within a replication share the donor block, so they are
        % not independent, and the three identities have different variances
        % from one another.  A theoretical standard error of the form
        % sqrt(2/N) gets both facts wrong.  Working replication by
        % replication and taking the spread across replications measures the
        % actual sampling error of each quantity, dependence included.
        rM = zeros(rep,1); rV = zeros(rep,1);
        oM = zeros(rep,1); oV = zeros(rep,1);
        eX = zeros(rep,1); sG = zeros(rep,1);
        for t = 1:rep
            rs = RandStream('mt19937ar','Seed', cfg.seed + 1e6*g + 1e3*j + t);
            Sd = bbSample(n,        P, rs);
            Sr = bbSample(cfg.nRecv,P, rs);

            [~,J] = bbNNDist(Sr.XB, Sd.XB, P);

            % attribute a = first exclusive attribute of block 1
            xTrue  = Sr.XQ{1}(:,1);
            rDon   = Sd.r{1}(J,1);        % r_a at the donor backbone
            vDon   = Sd.XQ{1}(J,1);       % raw observed value at the donor
            rRec   = Sr.r{1}(:,1);        % r_a at the receiver backbone
            sdRec2 = Sr.sd{1}(:,1).^2;
            sdDon2 = Sd.sd{1}(J,1).^2;

            rM(t) = mean((xTrue - rDon).^2);
            rV(t) = mean((xTrue - vDon).^2);
            oM(t) = mean(sdRec2 + (rDon - rRec).^2);
            oV(t) = mean(sdRec2 + sdDon2 + (rDon - rRec).^2);
            eX(t) = mean((rDon - rRec).^2);
            sG(t) = mean(sdDon2);
        end

        riskMed(g,j) = mean(rM);   riskVal(g,j) = mean(rV);
        oracMed(g,j) = mean(oM);   oracVal(g,j) = mean(oV);
        excess(g,j)  = mean(eX);   sigDon(g,j)  = mean(sG);

        % Relative discrepancy of each identity, with its own standard error
        % estimated from the spread across replications.
        dM = (rM - oM)./oM;
        dV = (rV - oV)./oV;
        dD = ((rV - rM) - sG)./sG;
        relMed(g,j)  = mean(dM);  seMed(g,j)  = std(dM)/sqrt(rep);
        relVal(g,j)  = mean(dV);  seVal(g,j)  = std(dV)/sqrt(rep);
        relDiff(g,j) = mean(dD);  seDiff(g,j) = std(dD)/sqrt(rep);
    end
end

% ---- fit the rate of the excess --------------------------------------
jFit  = max(1,ceil((1-cfg.fitFrom)*nN)) : nN;
x     = log(nGrid(jFit)).';
X     = [ones(numel(jFit),1), x];
slope = nan(nReg,1);
for g = 1:nReg
    y = log(excess(g,jFit)).';
    bh = X\y;
    slope(g) = bh(2);
end

% ---- identity checks -------------------------------------------------
% Reported as the largest absolute discrepancy over the grid, and separately
% as the largest discrepancy measured in units of its own standard error.
% The second is the one that matters: an identity that holds exactly should
% produce discrepancies of order one standard error regardless of how large
% that standard error happens to be.
errMed  = max(abs(relMed), [],2);
errVal  = max(abs(relVal), [],2);
errDiff = max(abs(relDiff),[],2);
zMed    = max(abs(relMed)./seMed, [],2);
zVal    = max(abs(relVal)./seVal, [],2);
zDiff   = max(abs(relDiff)./seDiff,[],2);

tbl = table(string({regimes.name}).', betaV, dBV, ...
            -2*betaV./dBV, slope, floorS2, ...
            100*errMed, 100*errVal, 100*errDiff, zMed, zVal, zDiff, ...
    'VariableNames',{'regime','beta','d_B','theory','slope','floor_s2', ...
                     'errMed_pct','errVal_pct','errDiff_pct', ...
                     'zMed','zVal','zDiff'});

out = struct('tbl',tbl,'nGrid',nGrid,'riskMed',riskMed,'riskVal',riskVal, ...
             'oracMed',oracMed,'oracVal',oracVal,'excess',excess, ...
             'sigDon',sigDon,'floorS2',floorS2,'slope',slope, ...
             'relMed',relMed,'relVal',relVal,'relDiff',relDiff, ...
             'seMed',seMed,'seVal',seVal,'seDiff',seDiff, ...
             'beta',betaV,'dB',dBV,'cfg',cfg,'regimes',regimes,'jFit',jFit);
end
