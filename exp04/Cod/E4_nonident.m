function out = E4_nonident(cfg)
%E4_NONIDENT  Experiment E4: Proposition P5 made operational.
%
%   OUT = E4_NONIDENT(CFG) builds the two populations of Proposition P5 --
%   one satisfying Cond_CI and one violating it, with identical block-wise
%   marginals -- and measures three things.
%
%   (1) THAT NO PROCEDURE SEPARATES THEM.  A two-sample energy test with a
%       permutation null is applied to the observable data of each block.
%       P5 says the two populations induce the same law on anything the
%       design produces, so the test must have no power: its p-values must
%       be uniform and its rejection rate must equal the nominal level.
%       Reporting a rejection rate turns "no procedure can distinguish them"
%       from an assertion into a measurement.
%
%   (2) WHAT ACTUALLY DIFFERS.  The conditional joint law of the exclusive
%       vectors of two blocks, at a fixed backbone value.  This is the object
%       the design never observes: no observation carries coordinates of two
%       blocks at once.
%
%   (3) WHAT IT COSTS.  A joint tail probability -- the chance that two
%       attributes measured in different blocks both exceed their upper
%       quintile.  The reconstruction converges to the value this quantity
%       takes under the product coupling, whatever the population that
%       generated the data.  When the truth is the dependent population, the
%       procedure therefore converges, with growing confidence, to a wrong
%       answer, and every block-wise diagnostic remains impeccable while it
%       does so.
%
%   CFG fields (optional): d, pm, theta, nGrid, nRecv, nTest, nPerm, nRep,
%   nOracle, quantile, seed, verbose.

if nargin < 1, cfg = struct(); end
dflt = struct('d',2,'pm',3,'theta',0.45, ...
              'nGrid',round(logspace(2,5,8)), 'nRecv',20000, ...
              'nTest',250,'nPerm',200,'nRep',400, ...
              'nOracle',2e6,'quantile',0.80,'seed',31415,'verbose',true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

base = {'d',cfg.d,'b',cfg.d,'M',2,'pm',cfg.pm,'meanType','smooth', ...
        'beta',1,'sigma0',0.5,'hetero',0.4,'rhoWithin',0.3,'sdShared',true};
Pperp = bbParams(base{:},'theta',0);            % Cond_CI holds
Pprime = bbParams(base{:},'theta',cfg.theta);   % Cond_CI fails

% The cross-block correlation is bounded by positive definiteness:
% theta < (1 + (pm-1) rhoWithin)/pm.  bbParams raises an error beyond it.
if cfg.verbose
    fprintf('[E4] two populations built: theta = 0 and theta = %.2f\n',cfg.theta);
end

% ---- block-wise summaries, which must agree --------------------------
rs = RandStream('mt19937ar','Seed',cfg.seed);
Sp = bbSample(cfg.nOracle,Pperp, rs);
rs = RandStream('mt19937ar','Seed',cfg.seed);
Sq = bbSample(cfg.nOracle,Pprime,rs);

% The first five rows are computable from the observed data; the sixth is
% not, because no observation carries attributes of two blocks at once.  The
% table is built so that the two columns agree on the first five and differ
% only on the sixth -- which is the whole content of Proposition P5.
qty = ["mean, block 1 attribute 1"; "sd, block 1 attribute 1"; ...
       "corr within block 1";       "mean, block 2 attribute 1"; ...
       "sd, block 2 attribute 1";   "corr ACROSS blocks (unobservable)"];
val = @(S) [mean(S.XQ{1}(:,1)); std(S.XQ{1}(:,1)); ...
            corr(S.XQ{1}(:,1),S.XQ{1}(:,2)); ...
            mean(S.XQ{2}(:,1)); std(S.XQ{2}(:,1)); ...
            corr(S.XQ{1}(:,1),S.XQ{2}(:,1))];
summ = table(qty, val(Sp), val(Sq), ...
    'VariableNames',{'quantity','under_Pperp','under_Pprime'});

% ---- the joint tail probability, under each population ---------------
t1 = quantile(Sp.XQ{1}(:,1), cfg.quantile);
t2 = quantile(Sp.XQ{2}(:,1), cfg.quantile);
tailPerp  = mean(Sp.XQ{1}(:,1) > t1 & Sp.XQ{2}(:,1) > t2);
tailPrime = mean(Sq.XQ{1}(:,1) > t1 & Sq.XQ{2}(:,1) > t2);

% conditional joint at one backbone value, for the scatter panel
z0  = 0.5*ones(1,cfg.d);
n0  = 4000;
rs0 = RandStream('mt19937ar','Seed',cfg.seed+99);
cond = struct();
for kk = 1:2
    if kk == 1, Pk = Pperp; else, Pk = Pprime; end
    Zc = repmat(z0,n0,1);
    Rc = bbOracleMean(Zc,Pk);
    Gc = bbOracleSd(Zc,Pk);
    Ec = randn(rs0,n0,Pk.q) * Pk.cholE.';
    Xc = Rc + Gc.*Ec;
    cond(kk).q1 = Xc(:,1);
    cond(kk).q2 = Xc(:,cfg.pm+1);
    % Exact conditional moments of the pair, for the probability ellipses.
    % Overlaid scatters of two nearly identical clouds read as one cloud;
    % the ellipses are what make the difference between the two populations
    % legible, and they are available in closed form.
    % The two attributes plotted belong to DIFFERENT blocks, and sdShared
    % equalises the dispersion direction only WITHIN a block.  Their
    % conditional standard deviations therefore differ, and the covariance
    % is D * C * D with D = diag(sd1,sd2) rather than a scalar multiple of
    % C.  Using a common sd would draw a circle for the independent law and
    % misreport the shape of both.
    sd1 = Gc(1,1);  sd2 = Gc(1,cfg.pm+1);
    Cc  = [1, Pk.theta; Pk.theta, 1];
    D   = diag([sd1 sd2]);
    cond(kk).mu  = [Rc(1,1), Rc(1,cfg.pm+1)];
    cond(kk).Sig = D * Cc * D;
    cond(kk).sd  = [sd1 sd2];
end
clear Sp Sq

% ---- (1) the two-sample test ------------------------------------------
if cfg.verbose, fprintf('[E4] two-sample tests: %d replications\n',cfg.nRep); end
pvals = nan(cfg.nRep,2);
for r = 1:cfg.nRep
    rsA = RandStream('mt19937ar','Seed',cfg.seed+2e5+r);
    rsB = RandStream('mt19937ar','Seed',cfg.seed+4e5+r);
    Sa = bbSample(cfg.nTest,Pperp, rsA);
    Sb = bbSample(cfg.nTest,Pprime,rsB);
    for m = 1:2
        A = [Sa.XB, Sa.XQ{m}];
        B = [Sb.XB, Sb.XQ{m}];
        pvals(r,m) = energyPermTest(A,B,cfg.nPerm,rsA);
    end
end

lev  = [0.05 0.10 0.20];
rej  = nan(numel(lev),2);
for m = 1:2
    for i = 1:numel(lev), rej(i,m) = mean(pvals(:,m) < lev(i)); end
end

% ---- (3) the reconstruction -------------------------------------------
if cfg.verbose, fprintf('[E4] reconstruction over the sample-size grid\n'); end
nGrid = cfg.nGrid(:).';
tailRec = nan(1,numel(nGrid));
tailSE  = nan(1,numel(nGrid));   % binomial standard error of the estimate
corrRec = nan(1,numel(nGrid));   % cross-correlation of the imputed pair
for j = 1:numel(nGrid)
    n  = nGrid(j);
    rs = RandStream('mt19937ar','Seed',cfg.seed+7e5+j);
    % receivers: only their backbones are used
    Sr = bbSample(cfg.nRecv,Pprime,rs);
    % two donor blocks, drawn from the DEPENDENT population, each carrying
    % only its own exclusive attributes, as Condition D1 requires
    S1 = bbSample(n,Pprime,rs);
    S2 = bbSample(n,Pprime,rs);
    j1 = knnsearch(createns(S1.XB,'NSMethod','kdtree'), Sr.XB, 'K',1);
    j2 = knnsearch(createns(S2.XB,'NSMethod','kdtree'), Sr.XB, 'K',1);
    v1 = S1.XQ{1}(j1,1);   v2 = S2.XQ{2}(j2,1);
    hit = (v1 > t1) & (v2 > t2);
    tailRec(j) = mean(hit);
    tailSE(j)  = sqrt(tailRec(j)*(1-tailRec(j))/numel(hit));
    corrRec(j) = corr(v1,v2);
end

out = struct('summary',summ,'pvals',pvals,'levels',lev,'rejection',rej, ...
             'tailPerp',tailPerp,'tailPrime',tailPrime,'tailRec',tailRec, ...
             'tailSE',tailSE,'corrRec',corrRec, ...
             'corrPerp',summ.under_Pperp(6),'corrPrime',summ.under_Pprime(6), ...
             'thetaCond',cfg.theta, ...
             'nGrid',nGrid,'thr',[t1 t2],'cond',cond,'z0',z0, ...
             'Pperp',Pperp,'Pprime',Pprime,'cfg',cfg);
end

% =======================================================================
function p = energyPermTest(A,B,nPerm,rs)
%ENERGYPERMTEST  Two-sample energy test with a permutation null.
%
%   The pooled distance matrix is formed once and every permutation is then
%   an indexing operation on it, rather than a fresh distance computation.
%   Without this the test dominates the cost of the whole experiment.
nA = size(A,1);  nB = size(B,1);  N = nA + nB;
D  = pdist2([A;B],[A;B]);
idx0 = [true(nA,1); false(nB,1)];
obs  = estat(D,idx0);
cnt  = 0;
for t = 1:nPerm
    pm = randperm(rs,N);
    id = false(N,1); id(pm(1:nA)) = true;
    if estat(D,id) >= obs, cnt = cnt + 1; end
end
p = (cnt + 1)/(nPerm + 1);
end

function e = estat(D,id)
a = id; b = ~id;
e = 2*mean(mean(D(a,b))) - mean(mean(D(a,a))) - mean(mean(D(b,b)));
end

% =======================================================================
function R = bbOracleMean(Z,P)
%BBORACLEMEAN  Conditional mean r(z), replicating the map used by bbSample.
U = Z * P.Wr.' - P.tr.';
switch P.meanType
    case 'weier'
        R = zeros(size(U));
        for k = 0:P.weierTerms
            R = R + P.weierBase^(-k*P.beta) * ...
                cos(2*pi * P.weierBase^k * U + P.phir.');
        end
        R = (P.Lmean * P.sr.') .* R;
    case 'abspow'
        R = (P.Lmean * P.sr.') .* (abs(U).^P.beta);
    case 'smooth'
        R = (P.Lmean * P.sr.') .* sin(2*pi*U) / (2*pi);
end
end

function SD = bbOracleSd(Z,P)
%BBORACLESD  Conditional standard deviation s(z), as used by bbSample.
if P.hetero == 0
    SD = P.sigma0 * ones(size(Z,1),P.q);
else
    SD = P.sigma0 * (1 + P.hetero * sin(2*pi * (Z * P.Wsig.')));
end
end
