function out = E1_nnRate(regimes, cfg)
%E1_NNRATE  Experiment E1: nearest-neighbour distance rate (verifies T1).
%
%   OUT = E1_NNRATE(REGIMES,CFG) estimates E[R_n^q] for a grid of donor
%   sample sizes n and fits
%
%       log E[R_n^q]  =  a  -  (q/d_B) * log n
%
%   by ordinary least squares on the upper part of the grid.  Theorem T1
%   predicts slope -q/d_B, where d_B is the GEOMETRIC dimension of the
%   backbone support -- not its cardinality b.  The whole point of the
%   experiment is the regime with b = 10 and d = 2, where the fitted slope
%   must come out near -1/2 and not near -1/10.
%
%   REGIMES is a struct array with fields 'name' and 'params' (a cell array
%   of name/value pairs passed to BBPARAMS).
%
%   CFG fields (all optional):
%     nGrid     vector of donor sample sizes        [round(logspace(2,5,10))]
%     nRecv     receivers drawn per replication     [500]
%     repMin    minimum replications per (regime,n) [20]
%     repMax    maximum replications per (regime,n) [200]
%     qList     moments to fit                      [1 2]
%     fitFrom   fraction of the grid used to fit    [0.5]
%     seed      base RNG seed for sampling          [12345]
%     verbose   print progress                      [true]
%
%   OUT is a struct with a summary table OUT.tbl and the raw curves.

if nargin < 2, cfg = struct(); end
dflt = struct('nGrid',   round(logspace(2,5,10)), ...
              'nRecv',   500, ...
              'repMin',  20, ...
              'repMax',  200, ...
              'qList',   [1 2], ...
              'fitFrom', 0.5, ...
              'seed',    12345, ...
              'verbose', true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

nGrid = cfg.nGrid(:).';
nN    = numel(nGrid);
nQ    = numel(cfg.qList);
nReg  = numel(regimes);

MOM = nan(nReg,nN,nQ);      % estimated E[R_n^q]
SE  = nan(nReg,nN,nQ);      % standard error of that estimate
REP = nan(nReg,nN);
DB  = nan(nReg,1);

allR = cell(nReg,nN);

for g = 1:nReg
    P     = bbParams(regimes(g).params{:});
    DB(g) = P.d;
    if cfg.verbose
        fprintf('[E1] regime %d/%d: %-22s d=%d b=%d %s/%s\n', ...
            g,nReg,regimes(g).name,P.d,P.b,P.embedding,P.geometry);
    end

    for j = 1:nN
        n = nGrid(j);
        % Keep the total number of (receiver,replication) pairs roughly
        % constant across n, so every point of the curve carries comparable
        % Monte-Carlo precision at roughly constant cost.
        rep = min(cfg.repMax, max(cfg.repMin, ceil(cfg.repMax*sqrt(1e3/n))));
        REP(g,j) = rep;

        Rall = zeros(rep*cfg.nRecv,1);
        for t = 1:rep
            rs = RandStream('mt19937ar','Seed', cfg.seed + 1e6*g + 1e3*j + t);
            Sd = bbSample(n,        P, rs);      % donors
            Sq = bbSample(cfg.nRecv,P, rs);      % receivers, independent
            Rv = bbNNDist(Sq.XB, Sd.XB, P);
            Rall((t-1)*cfg.nRecv + (1:cfg.nRecv)) = Rv;
        end

        allR{g,j} = Rall;
        for k = 1:nQ
            v            = Rall.^cfg.qList(k);
            MOM(g,j,k)   = mean(v);
            SE(g,j,k)    = std(v)/sqrt(numel(v));
        end
    end
end

% ---- fit the slopes ---------------------------------------------------
jFit  = max(1,ceil((1-cfg.fitFrom)*nN)) : nN;
x     = log(nGrid(jFit)).';
X     = [ones(numel(jFit),1), x];

slope = nan(nReg,nQ); slopeCI = nan(nReg,nQ,2); dHat = nan(nReg,nQ);
for g = 1:nReg
    for k = 1:nQ
        y   = log(squeeze(MOM(g,jFit,k))).';
        bhat = X\y;
        res  = y - X*bhat;
        dfe  = numel(y) - 2;
        s2   = (res.'*res)/max(dfe,1);
        C    = s2 * inv(X.'*X);                                    %#ok<MINV>
        se   = sqrt(C(2,2));
        tcrit = 2.0;                       % ~95% for the grid sizes used here
        slope(g,k)      = bhat(2);
        slopeCI(g,k,:)  = bhat(2) + tcrit*se*[-1 1];
        dHat(g,k)       = -cfg.qList(k)/bhat(2);
    end
end

% ---- summary table ----------------------------------------------------
name = string({regimes.name}).';
tbl  = table(name, DB, 'VariableNames',{'regime','d_true'});
for k = 1:nQ
    q = cfg.qList(k);
    tbl.(sprintf('slope_q%d',q))     = slope(:,k);
    tbl.(sprintf('theory_q%d',q))    = -q./DB;
    tbl.(sprintf('dHat_q%d',q))      = dHat(:,k);
end

out = struct('tbl',tbl,'nGrid',nGrid,'moment',MOM,'se',SE,'rep',REP, ...
             'slope',slope,'slopeCI',slopeCI,'dHat',dHat,'dTrue',DB, ...
             'qList',cfg.qList,'cfg',cfg,'regimes',regimes,'R',{allR}, ...
             'jFit',jFit);
end
