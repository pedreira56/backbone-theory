function out = E7_completedBase(cfg)
%E7_COMPLETEDBASE  Experiment E7: Theorem T5, the completed dataset.
%
%   OUT = E7_COMPLETEDBASE(CFG) measures how far the empirical measure of a
%   completed dataset is from the product coupling, as a function of the
%   number N of receivers and the size n of the donor block.
%
%   WHAT IS CHECKED
%   Theorem T5 bounds, for f bounded by one and 1-Lipschitz,
%
%       E[ ( int f dPhat_N  -  int f dPperp )^2 ]  <=  4/N + Delta_n^2 + 8 pi_n,
%
%   and concludes that the completed dataset converges weakly to the product
%   coupling whenever N and n both grow, WITH NO RELATION REQUIRED BETWEEN
%   THEM.  Three things follow, and each is measured here.
%
%   (1) Growing N alone is not enough.  With n fixed, the mean over the
%       completed rows converges to a quantity that differs from the limit by
%       the approximation error of the matching, so the squared error settles
%       on a plateau set by n and does not fall further.
%   (2) Growing n alone is not enough either: the Monte-Carlo term 4/N stays.
%   (3) With both growing, in any relation, the error vanishes -- including
%       along N = 100 n, where receivers outnumber donors by two orders of
%       magnitude and every donor serves many receivers.
%
%   The quantity pi_n, the probability that two distinct receivers select the
%   same donor, is the one that carries the dependence between completed rows
%   and is the reason the theorem is not a plain law of large numbers.  It is
%   computed exactly here, from the multiplicities of the selected donors,
%   rather than by enumerating pairs.
%
%   THE POPULATION
%   The dependent population of experiment E4, in which Cond_CI fails.  The
%   limit of the completed dataset is then the PRODUCT coupling and not the
%   population itself, so the experiment also exhibits the gap of Section 9.6
%   as the distance between two horizontal lines.
%
%   THE TEST FUNCTIONS
%   Three functions of the pair formed by the first exclusive attribute of
%   each block, each bounded by one in absolute value and 1-Lipschitz for the
%   Euclidean norm, as Theorem T5 requires.  Each depends on attributes of
%   BOTH blocks, which is what makes the cross-block coupling visible: a
%   function of one block alone would be reproduced exactly by any coupling.
%
%       f1(u,v) = g(u) g(v) / sqrt(2),   g(x) = min(max(x - t,0),1)
%                 a smoothed joint upper tail, the quantity of Section 9.6
%       f2(u,v) = cos(u - v) / sqrt(2)
%       f3(u,v) = 1 - min(|u - v|,2)/2
%
%   CFG fields (optional): nGrid, NGrid, repBase, repMin, repMax, nOracle,
%   quantile, theta, seed, verbose.

if nargin < 1, cfg = struct(); end
dflt = struct('nGrid',round(logspace(2,5,7)), 'NGrid',round(logspace(2,5,7)), ...
              'repBase',400,'repMin',8,'repMax',60,'nOracle',2e6, ...
              'quantile',0.80,'theta',0.45,'seed',56789,'verbose',true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

% hetero = 0 keeps the conditional dispersion free of the backbone, which is
% what puts the anchoring constant of the donor block in closed form:
%     Lambda_2^2 = sum over the block's attributes of (Lmean * sr_a)^2,
% exactly as in experiment E5.  The bound of Theorem T5 can then be evaluated
% numerically and compared with the measured error, rather than only its
% rate being checked.
base = {'d',2,'b',2,'M',2,'pm',3,'meanType','smooth','beta',1, ...
        'sigma0',0.5,'hetero',0,'rhoWithin',0.3,'sdShared',true};
Pperp  = bbParams(base{:},'theta',0);           % Cond_CI holds
Pprime = bbParams(base{:},'theta',cfg.theta);   % Cond_CI fails: the truth

% ---- thresholds and the three limits -----------------------------------
rs = RandStream('mt19937ar','Seed',cfg.seed);
S0 = bbSample(cfg.nOracle, Pperp, rs);
t1 = quantile(S0.XQ{1}(:,1), cfg.quantile);
t2 = quantile(S0.XQ{2}(:,1), cfg.quantile);
F  = @(u,v) testFuns(u,v,t1,t2);

thetaPerp = mean(F(S0.XQ{1}(:,1), S0.XQ{2}(:,1)), 1);   % the limit of T5
clear S0
rs = RandStream('mt19937ar','Seed',cfg.seed);
S1 = bbSample(cfg.nOracle, Pprime, rs);
thetaTrue = mean(F(S1.XQ{1}(:,1), S1.XQ{2}(:,1)), 1);   % the population value
clear S1
nF = numel(thetaPerp);

if cfg.verbose
    fprintf('[E7] limits of the three test functions\n');
    for k = 1:nF
        fprintf('     f%d: product coupling %.6f   population %.6f   gap %+.6f\n', ...
                k, thetaPerp(k), thetaTrue(k), thetaTrue(k)-thetaPerp(k));
    end
end

% ---- the anchoring constant of the donor block, in closed form ---------
idx2  = (Pprime.blockId == 2);
lam2  = sum((Pprime.Lmean * Pprime.sr(idx2)).^2);
if cfg.verbose, fprintf('[E7] Lambda_2^2 = %.4f\n',lam2); end

% ---- the grid ----------------------------------------------------------
nG = cfg.nGrid(:).';  NG = cfg.NGrid(:).';
MSE = nan(numel(nG),numel(NG),nF);
BIAS= nan(numel(nG),numel(NG),nF);
PI  = nan(numel(nG),numel(NG));
R2  = nan(numel(nG),numel(NG));     % E[R^2], the matching moment
REP = nan(numel(nG),numel(NG));

for a = 1:numel(nG)
    n = nG(a);
    for b = 1:numel(NG)
        N = NG(b);
        rep = min(cfg.repMax, max(cfg.repMin, ...
                  round(cfg.repBase*sqrt(1e3/max(n,N)))));
        REP(a,b) = rep;
        acc  = zeros(rep,nF);
        accP = zeros(rep,1);
        accR = zeros(rep,1);
        for r = 1:rep
            rsr = RandStream('mt19937ar','Seed', cfg.seed + 1e6*a + 1e3*b + r);
            % Receivers: their backbone and their own block is observed.
            Sr = bbSample(N, Pprime, rsr);
            % Donor block, an independent sample of the same population,
            % carrying only its own exclusive attributes (Condition D1).
            Sd = bbSample(n, Pprime, rsr);
            [j,dist] = knnsearch(createns(Sd.XB,'NSMethod','kdtree'), Sr.XB, 'K',1);
            accR(r) = mean(dist.^2);
            v  = Sd.XQ{2}(j,1);                 % the transferred value
            acc(r,:) = mean(F(Sr.XQ{1}(:,1), v), 1);
            % pi_n, exactly: the probability that two distinct receivers
            % selected the same donor.  With c_j receivers assigned to donor
            % j, the number of ordered pairs sharing a donor is sum c_j(c_j-1).
            c = accumarray(j, 1, [n 1]);
            accP(r) = sum(c.*(c-1)) / (N*(N-1));
        end
        d2 = (acc - thetaPerp).^2;
        MSE(a,b,:) = mean(d2,1);
        BIAS(a,b,:) = mean(acc,1) - thetaPerp;
        PI(a,b)    = mean(accP);
        R2(a,b)    = mean(accR);
        if cfg.verbose
            fprintf('[E7] n = %6d  N = %6d  rep = %3d   rmse(f1) = %.3e   pi = %.3e\n', ...
                    n,N,rep,sqrt(MSE(a,b,1)),PI(a,b));
        end
    end
end

% The bound of Theorem T5, term by term.
D2    = lam2 * R2;                       % Delta_n^2
BOUND = 4./repmat(NG,numel(nG),1) + D2 + 8*PI;

out = struct('nGrid',nG,'NGrid',NG,'mse',MSE,'bias',BIAS,'pi',PI,'rep',REP, ...
             'R2',R2,'lam2',lam2,'Delta2',D2,'bound',BOUND, ...
             'thetaPerp',thetaPerp,'thetaTrue',thetaTrue,'thr',[t1 t2], ...
             'Pperp',Pperp,'Pprime',Pprime,'cfg',cfg, ...
             'fNames',{{'f1 joint tail','f2 cosine','f3 distance'}});
end

% =======================================================================
function V = testFuns(u,v,t1,t2)
%TESTFUNS  Three functions bounded by one and 1-Lipschitz on the pair (u,v).
g1 = min(max(u - t1,0),1);
g2 = min(max(v - t2,0),1);
V  = [ g1.*g2/sqrt(2), ...
       cos(u - v)/sqrt(2), ...
       1 - min(abs(u - v),2)/2 ];
end
