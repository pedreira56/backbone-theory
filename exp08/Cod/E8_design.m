function out = E8_design(cfg)
%E8_DESIGN  Experiment E8: a worked design problem with costs and capacity.
%
%   OUT = E8_DESIGN(CFG) carries the formulation of Section 8 through a
%   concrete design question, by exhaustive enumeration of the candidates.
%
%   THE SETTING
%   Specimens of a material receive a subset of six non-destructive
%   measurements and are then destroyed in one of three mechanical tests.
%   Each specimen yields the measurements it received and the outcomes of ONE
%   test; no specimen yields the outcomes of two tests.  The three tests are
%   the blocks, the non-destructive measurements are the candidates for the
%   backbone, and the question is which of the sixty-four subsets to measure.
%
%   WHAT THIS ADDS TO E6
%   (1) The population is neither linear nor Gaussian: the conditional means
%       are non-linear in a two-dimensional latent material state and the
%       conditional dispersion is heteroscedastic with heavy-tailed
%       innovations.  No quantity is available in closed form.
%   (2) Bench time is finite and shared.  Time spent measuring is time not
%       spent testing, so enlarging the backbone REDUCES the number of
%       specimens.  A larger backbone is more informative per specimen and
%       leaves fewer specimens, and the criterion sees both effects.
%   (3) Monetary cost per measurement enters separately, on a scale that must
%       be commensurate with the loss, which the experiment makes explicit by
%       varying it rather than fixing it.
%
%   THE CRITERION
%   The implementable one: the risk of transferring an OBSERVED value, which
%   is what a practitioner can compute and what Section 8.1 identifies as the
%   rule actually used.  For a candidate B,
%
%       Risk(B) = sum over outcomes of
%                 E[ ( Y^{receiver} - Y^{donor nearest in B} )^2 ] / Var(Y),
%
%   evaluated at the donor block sizes the capacity allows for that B.  Its
%   large-sample limit is twice the residual conditional variance, the part
%   no amount of bench time removes; the excess over that limit is the price
%   of imperfect matching.  Both are reported, and their sum is the exhibit.
%
%   VARIANCE REDUCTION
%   Every candidate, every capacity and the large-sample limit are evaluated
%   on the SAME receiver samples, drawn from streams seeded independently of
%   the candidate.  The comparison between candidates is then paired, and the
%   difference between a risk and its limit is not swamped by the Monte-Carlo
%   error of two independent estimates of nearly equal quantities.
%
%   CFG fields (optional): nPilot, nRecv, nFloor, repRisk, repFloor, caps,
%   lamGrid, seed, verbose.

if nargin < 1, cfg = struct(); end
dflt = struct('nPilot',3e5,'nRecv',4000,'nFloor',20000,'repRisk',5, ...
              'repFloor',4,'caps',[15 60 300],'lamGrid',linspace(0,0.4,201), ...
              'seed',8191,'verbose',true);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(cfg,fn{i}), cfg.(fn{i}) = dflt.(fn{i}); end
end

% ---- the six candidate measurements ------------------------------------
% name, monetary cost, bench minutes per specimen, measurement noise
M.name = {'dimensions','density','hardness','ultrasonic','roughness','conductivity'};
M.cost = [1    2    3    5    4    3   ];
M.time = [5    10   15   40   30   15  ];
M.sig  = [0.10 0.08 0.12 0.05 0.14 0.09];
p = numel(M.cost);

% ---- the three destructive tests ---------------------------------------
B3.name = {'tensile','impact','fatigue'};
B3.time = [20 10 60];                      % bench minutes per specimen
nOut = 2;  M3 = numel(B3.name);  q = M3*nOut;
blockOf = repelem(1:M3, nOut);

% ---- pilot sample: standardisation --------------------------------------
% Standardising the measurements is not cosmetic: Euclidean distance in
% backbone space would otherwise depend on the units each measurement
% happens to be reported in.  Standardising the outcomes puts the criterion
% in units of outcome variance, which is what makes a monetary cost per
% measurement commensurable with it at all.
rsP = RandStream('mt19937ar','Seed',cfg.seed);
[Xp,Yp] = drawPop(cfg.nPilot, M, rsP);
std_ = struct('mx',mean(Xp,1),'sx',std(Xp,0,1),'sy',std(Yp,0,1));
clear Xp Yp

% ---- enumerate the 2^p candidate backbones ------------------------------
nS = 2^p;
subsets = cell(nS,1);  card = zeros(nS,1);  costB = zeros(nS,1);  timeB = zeros(nS,1);
for s = 0:nS-1
    b = find(bitget(s,1:p));
    subsets{s+1} = b;  card(s+1) = numel(b);
    costB(s+1) = sum(M.cost(b));  timeB(s+1) = sum(M.time(b));
end

seedQ = cfg.seed + 777777;                 % the shared receiver streams
nCap  = numel(cfg.caps);
nDon  = zeros(nS,M3,nCap);
risk  = nan(nS,nCap);
floorR= nan(nS,1);
if cfg.verbose
    fprintf('[E8] enumerating %d candidate backbones, %d capacity levels\n',nS,nCap);
end

for s = 1:nS
    b = subsets{s};
    for c = 1:nCap
        nDon(s,:,c) = floor(60*cfg.caps(c) ./ (B3.time + sum(M.time(b))));
    end
    if isempty(b)
        % With no backbone there is nothing to match on: the transferred
        % value is an independent draw from the marginal, whose squared
        % error is twice the marginal variance, hence two per outcome.
        floorR(s) = 2*q;  risk(s,:) = 2*q;
        continue
    end
    floorR(s) = sumRisk(b, repmat(cfg.nFloor,1,M3), cfg.nRecv, cfg.repFloor, ...
                        M, B3, blockOf, std_, cfg.seed + 7919*s, seedQ);
    for c = 1:nCap
        risk(s,c) = sumRisk(b, nDon(s,:,c), cfg.nRecv, cfg.repRisk, ...
                            M, B3, blockOf, std_, ...
                            cfg.seed + 7919*s + 131071*c, seedQ);
    end
    if cfg.verbose && mod(s,8) == 0
        fprintf('   %3d / %3d\n',s,nS);
    end
end

% ---- the recommendation, as the price of a measurement rises -----------
lam = cfg.lamGrid(:).';
bestIdx = zeros(numel(lam),nCap);
for c = 1:nCap
    for i = 1:numel(lam)
        [~,bestIdx(i,c)] = min(risk(:,c) + lam(i)*costB);
    end
end

out = struct('subsets',{subsets},'card',card,'costB',costB,'timeB',timeB, ...
             'risk',risk,'floor',floorR,'nDon',nDon,'caps',cfg.caps, ...
             'lam',lam,'bestIdx',bestIdx,'bestCard',card(bestIdx), ...
             'meas',M,'tests',B3,'blockOf',blockOf,'std',std_,'cfg',cfg);
end

% =======================================================================
function r = sumRisk(b, nm, nRecv, rep, M, B3, blockOf, std_, seedD, seedQ)
%SUMRISK  Summed value-transfer risk over the six outcome attributes.
acc = zeros(rep,1);
for t = 1:rep
    tot = 0;
    for m = 1:numel(B3.time)
        % Receivers: a stream that does not depend on the candidate, so that
        % all candidates are compared on the same specimens.
        rsQ = RandStream('mt19937ar','Seed', seedQ + 1013*t + m);
        [Xq,Yq] = drawPop(nRecv, M, rsQ);
        rsD = RandStream('mt19937ar','Seed', seedD + 1013*t + m);
        [Xd,Yd] = drawPop(max(nm(m),2), M, rsD);
        Zd = (Xd(:,b) - std_.mx(b)) ./ std_.sx(b);
        Zq = (Xq(:,b) - std_.mx(b)) ./ std_.sx(b);
        j  = knnsearch(createns(Zd,'NSMethod','kdtree'), Zq, 'K',1);
        idx = find(blockOf == m);
        dy  = (Yq(:,idx) - Yd(j,idx)) ./ std_.sy(idx);
        tot = tot + sum(mean(dy.^2,1));
    end
    acc(t) = tot;
end
r = mean(acc);
end

% =======================================================================
function [X,Y] = drawPop(n, M, rs)
%DRAWPOP  A two-dimensional latent material state, six non-linear noisy
%   measurements of it, and six destructive outcomes with heteroscedastic,
%   heavy-tailed conditional laws.
Z = rand(rs,n,2);
z1 = Z(:,1);  z2 = Z(:,2);

% The measurements.  Some read one latent coordinate cleanly, some mix the
% two, and the ultrasonic reading is the most precise but is not injective,
% so it identifies the state only in company.
f = [ z1, ...
      0.5*(z1+z2) + 0.3*sin(2*pi*z1), ...
      z2, ...
      0.5 + 0.5*sin(2*pi*z1).*cos(2*pi*z2), ...
      z1.*z2, ...
      z2 + 0.2*cos(2*pi*z1) ];
X = f + randn(rs,n,numel(M.sig)) .* M.sig;

% The outcomes: two per test, with non-linear conditional means, so that no
% single measurement serves the three tests equally well.
g = [ 1.2*z1 + 0.6*sin(2*pi*z2), ...
      0.8*z2.^2 + 0.4*z1, ...
      0.9*sin(2*pi*(z1+z2))/2 + 0.5*z2, ...
      1.0*z1.*z2 + 0.3*z1, ...
      0.7*cos(2*pi*z1) + 0.9*z2, ...
      0.6*z1 + 0.6*z2 + 0.5*sin(2*pi*z1.*z2) ];

% Heteroscedastic, heavy-tailed innovations: Student's t with eight degrees
% of freedom, rescaled to unit variance, with a dispersion that varies over
% the latent state and differs between outcomes.  The chi-square is built as
% a sum of eight squared normals rather than with RANDG or CHI2RND: those
% draw from the global stream, and the whole experiment must run off the
% stream it was given in order to reproduce exactly.
k   = size(g,2);
amp = [0.22 0.28 0.25 0.30 0.24 0.26];
sd  = amp .* (1 + 0.5*sin(2*pi*z1 + (0:k-1)*pi/3));
nu  = 8;                                      % must be a positive integer
W   = zeros(n,k);
for i = 1:nu
    W = W + randn(rs,n,k).^2;                 % chi-square with nu degrees
end
Tn  = randn(rs,n,k) ./ sqrt(W/nu);            % Student's t
Y   = g + sd .* (Tn / sqrt(nu/(nu-2)));       % unit-variance innovations
end
