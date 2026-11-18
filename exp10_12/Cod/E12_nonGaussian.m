function T = E12_nonGaussian(reps, outDir)
%E12_NONGAUSSIAN  Finite-sample backbone selection in two non-Gaussian,
%                 non-linear populations, with the criterion estimated
%                 DIRECTLY by hold-out validation.
%
%   T = E12_NONGAUSSIAN() runs the experiment with 50 replications per pilot
%   size and writes E12.mat and E12_table.txt to the folder of this file.
%   T = E12_NONGAUSSIAN(REPS) sets the replications.
%   T = E12_NONGAUSSIAN(REPS,OUTDIR) writes the output to OUTDIR.  By default
%   the output goes to the sibling Result folder when BBOUTDIR is on the path,
%   and to the folder of this file otherwise.
%
%   Self-contained: needs no .mat, and no other file except BBOUTDIR, which is
%   optional.  Requires KNNSEARCH.
%
%   WHY THIS EXPERIMENT EXISTS
%   E9, E10, E10b and E11 all select inside a linear-Gaussian population, in
%   which the conditional quantities of the criterion have closed forms and
%   the empirical selector estimates those forms.  A reviewer objected that
%   this validates the selector only where its functional form is correct by
%   construction.  Here the population is neither Gaussian nor linear in two
%   substantially different ways, and NO component of the criterion is
%   estimated separately.
%
%   WHAT IS ESTIMATED, AND HOW
%   In the calibration regime of Section 8.4 the risk of value transfer is
%   MEASURABLE, because the pilot observes every attribute on the same
%   individuals.  The pilot is split into pseudo-donors and pseudo-receivers;
%   for a candidate backbone each pseudo-receiver is matched to its nearest
%   pseudo-donor in the candidate's standardized coordinates; the donor's
%   outcome values are transferred; and the squared error against the
%   receiver's observed values is averaged.  That is the criterion itself, not
%   a surrogate for it.  No conditional mean, no anchoring constant, no
%   pairing term and no geometric dimension is estimated.  Two independent
%   splits are averaged, standardization is computed from the pilot alone, and
%   all 63 non-empty candidates are enumerated.
%
%   THE TARGET
%   The population criterion is computed by large Monte Carlo at the SAME
%   donor count the pilot's hold-out provides, so the empirical quantity and
%   its target refer to the same donor scale.  As in E11, the optimum may move
%   with that scale; that is Corollary C2, not an artefact.
%
%   THE TWO POPULATIONS
%     A  the destructive-testing population of Section 9.11: six noisy
%        non-linear readings of a two-dimensional latent state, six outcomes
%        with non-linear conditional means and heteroscedastic Student-t
%        innovations.
%     B  built for this experiment and deliberately unlike A: a bimodal
%        latent state, non-monotone measurements, asymmetric lognormal
%        measurement error, and outcome laws that are MIXTURES whose branch
%        probability depends on the latent state, so the conditional
%        distributions are themselves bimodal, with Laplace innovations whose
%        dispersion varies over the state.
%   The candidate space, the cost vector and the selector are identical in the
%   two, so the population is the only thing that differs.
%
%   WHAT TO READ IN THE OUTPUT
%   At price zero the criterion is FLAT near its optimum in both populations:
%   the best candidates differ by a few tenths of a percent, so the identity
%   of the argmin is not stable even under very large Monte Carlo, and exact
%   recovery would be reporting near-ties.  The primary quantity is therefore
%   the excess criterion of the selected design, with the fraction of
%   replications within one percent of the optimum beside it.  At a positive
%   price the optima separate, because candidates of different cardinality
%   differ in cost, and exact recovery becomes meaningful again.
%
%   SIMPLIFICATION RELATIVE TO E8, STATED SO IT IS NOT MISTAKEN FOR A RESULT
%   All three outcome blocks draw from one donor set of common size; there is
%   no bench-capacity allocation across tests.  The question here is whether
%   the empirical selector recovers the right design, not how capacity is
%   split.

if nargin < 1 || isempty(reps), reps = 50; end
if nargin < 2 || isempty(outDir)
    if exist('bbOutDir','file') == 2
        outDir = bbOutDir(mfilename('fullpath'));   % sibling Result folder
    else
        outDir = fileparts(mfilename('fullpath'));
    end
end
if isempty(outDir), outDir = pwd; end

cost  = [1 2 3 5 4 3];          % measurement cost, as in Section 9.11
q     = numel(cost);
ns    = [50 100 200 500 1000 2000];
nn    = numel(ns);
nSplit   = 2;                   % hold-out splits averaged per replication
nRecvBig = 20000;               % receivers for the population target
repTar   = 6;                   % target replications
tolW     = 0.01;                % "within tolerance of the optimum" band
nStd     = 200000;              % pilot for the population standardization
seedTar  = 90001;
seedRep  = 7e6;

cand = cell(2^q-1,1);  cardC = zeros(2^q-1,1);  costC = zeros(2^q-1,1);
for i = 1:2^q-1
    b = find(bitget(i, 1:q));
    cand{i} = b;  cardC(i) = numel(b);  costC(i) = sum(cost(b));
end
nC = numel(cand);

gens = {'A', 'B'};
res  = struct();

for ig = 1:2
    g = gens{ig};
    fprintf('\n======== generator %s ========\n', g);

    % ---- standardization, from a large pilot of this population ----------
    rsS = RandStream('mt19937ar','Seed', seedTar + 7);
    [Xs, Ys] = drawPop(nStd, g, rsS);
    std_ = struct('mx', mean(Xs,1), 'sx', std(Xs,0,1), 'sy', std(Ys,0,1));
    clear Xs Ys

    % ---- population criterion, at each donor count the pilots provide ----
    nD   = floor(ns/2);
    Rpop = zeros(nC, nn);
    t0 = tic;
    for iN = 1:nn
        for t = 1:repTar
            rsQ = RandStream('mt19937ar','Seed', seedTar + 1013*t);
            [Xq, Yq] = drawPop(nRecvBig, g, rsQ);
            rsD = RandStream('mt19937ar','Seed', seedTar + 5001*t + 31*iN);
            [Xd, Yd] = drawPop(nD(iN), g, rsD);
            for i = 1:nC
                Rpop(i,iN) = Rpop(i,iN) + ...
                    transferRisk(Xd, Yd, Xq, Yq, cand{i}, std_) / repTar;
            end
        end
        fprintf('  target at nD = %4d done (%.1f min)\n', nD(iN), toc(t0)/60);
    end

    % ---- the price, fixed exogenously from the scale of the risk range ---
    lamRef = (max(Rpop(:,nn)) - min(Rpop(:,nn))) / sum(cost);
    lams   = [0 lamRef];
    nL     = numel(lams);
    kStar  = zeros(nL, nn);  Ustar = zeros(nL, nn);
    for iN = 1:nn
        for iL = 1:nL
            [Ustar(iL,iN), kStar(iL,iN)] = min(Rpop(:,iN) + lams(iL)*costC);
        end
    end
    fprintf('  price lambda = %.4f per unit of measurement cost\n', lamRef);
    for iN = 1:nn
        s0 = sort(Rpop(:,iN));
        fprintf('    nD = %4d: optimum at price 0 = [%s] (U = %.4f, gap to 2nd %.4f); at price %.4f = [%s] (cost %g)\n', ...
            nD(iN), num2str(cand{kStar(1,iN)}), s0(1), s0(2)-s0(1), ...
            lamRef, num2str(cand{kStar(2,iN)}), costC(kStar(2,iN)));
    end

    % ---- replications of the empirical selector --------------------------
    HIT = zeros(nn,nL);  EXC = zeros(nn,nL);  CRD = zeros(nn,nL);  WIN = zeros(nn,nL);
    hitR = false(reps,nn,nL);  excR = zeros(reps,nn,nL);  crdR = zeros(reps,nn,nL);
    t0 = tic;  first = true;
    for iN = 1:nn
        n = ns(iN);  nd = nD(iN);  nr = n - nd;
        for r = 1:reps
            rsP = RandStream('mt19937ar','Seed', seedRep + 1000*iN + r);
            [X, Y] = drawPop(n, g, rsP);
            sd_ = struct('mx', mean(X,1), 'sx', std(X,0,1), 'sy', std(Y,0,1));
            sd_.sx(sd_.sx <= 0) = 1;  sd_.sy(sd_.sy <= 0) = 1;

            Rhat = zeros(nC,1);
            for s = 1:nSplit
                pr = randperm(rsP, n);
                di = pr(1:nd);  qi = pr(nd+1:nd+nr);
                for i = 1:nC
                    Rhat(i) = Rhat(i) + transferRisk(X(di,:), Y(di,:), ...
                                                     X(qi,:), Y(qi,:), cand{i}, sd_);
                end
            end
            Rhat = Rhat / nSplit;

            for iL = 1:nL
                [~, kh] = min(Rhat + lams(iL)*costC);
                hitR(r,iN,iL) = (kh == kStar(iL,iN));
                excR(r,iN,iL) = (Rpop(kh,iN) + lams(iL)*costC(kh)) / Ustar(iL,iN) - 1;
                crdR(r,iN,iL) = cardC(kh);
            end
            if first
                first = false;  el = toc(t0);
                fprintf('  first replication took %.2f s; projected total %.0f min\n', ...
                        el, el*reps*nn/60);
            end
        end
        for iL = 1:nL
            HIT(iN,iL) = mean(hitR(:,iN,iL));
            EXC(iN,iL) = mean(excR(:,iN,iL));
            WIN(iN,iL) = mean(excR(:,iN,iL) <= tolW);
            CRD(iN,iL) = median(crdR(:,iN,iL));
        end
        fprintf('  n = %5d done (%.1f min elapsed)\n', n, toc(t0)/60);
    end

    rg = struct();
    rg.Rpop = Rpop;   rg.lams = lams;   rg.kStar = kStar;   rg.Ustar = Ustar;
    rg.nD = nD;       rg.HIT = HIT;     rg.EXC = EXC;       rg.CRD = CRD;
    rg.WIN = WIN;     rg.tolW = tolW;
    rg.hitR = hitR;   rg.excR = excR;   rg.crdR = crdR;
    rg.bstar0 = cardC(kStar(1,:)).';    rg.bstarL = cardC(kStar(2,:)).';
    rg.Bstar0 = cand(kStar(1,:));       rg.BstarL = cand(kStar(2,:));
    res.(g) = rg;
end

n = ns(:);
T = table(n, res.A.nD(:), ...
          res.A.EXC(:,1), res.A.WIN(:,1), res.A.EXC(:,2), res.A.HIT(:,2), res.A.CRD(:,2), ...
          res.B.EXC(:,1), res.B.WIN(:,1), res.B.EXC(:,2), res.B.HIT(:,2), res.B.CRD(:,2), ...
    'VariableNames', {'n','nD', ...
                      'excA_p0','winA_p0','excA_pL','hitA_pL','bselA_pL', ...
                      'excB_p0','winB_p0','excB_pL','hitB_pL','bselB_pL'});
disp(' ');
disp('=== E12: direct empirical selection in two non-Gaussian populations ===');
disp(T);

save(fullfile(outDir,'E12.mat'),'T','res','ns','reps','cost','cand','cardC','costC', ...
     'nSplit','nRecvBig','repTar','-v7');
fid = fopen(fullfile(outDir,'E12_table.txt'),'w');
fprintf(fid,'E12 - direct empirical selection in two non-Gaussian, non-linear populations\n\n');
fprintf(fid,'the criterion is estimated by hold-out validation of the transfer itself;\n');
fprintf(fid,'no conditional mean, anchoring constant, pairing term or dimension is estimated\n');
fprintf(fid,'candidates: the %d non-empty subsets of the six measurements\n', nC);
fprintf(fid,'measurement cost: [%s]\n', num2str(cost));
fprintf(fid,'hold-out splits averaged per replication: %d; replications: %d\n', nSplit, reps);
fprintf(fid,'target: population criterion at the same donor count, %d receivers, %d replications\n', ...
        nRecvBig, repTar);
fprintf(fid,'exc = mean excess population criterion of the selected design\n');
fprintf(fid,'win = fraction of replications within %g%% of the optimum; the criterion is flat\n', 100*tolW);
fprintf(fid,'      near its optimum in these populations, so the loss, not the argmin, is the\n');
fprintf(fid,'      quantity that matters at price zero\n');
fprintf(fid,'hit = fraction recovering the optimal subset exactly; bsel = median cardinality\n\n');
for ig = 1:2
    g = gens{ig};  r = res.(g);
    fprintf(fid,'generator %s: price lambda = %.4f\n', g, r.lams(2));
    for iN = 1:nn
        fprintf(fid,'  nD = %4d: optimum at price 0 = [%s]; at price lambda = [%s]\n', ...
            r.nD(iN), num2str(r.Bstar0{iN}), num2str(r.BstarL{iN}));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E12.mat and E12_table.txt to:\n  %s\n', outDir);
end

% =======================================================================
function r = transferRisk(Xd, Yd, Xq, Yq, b, std_)
%TRANSFERRISK  Measured risk of value transfer on the candidate B.
%   Each receiver is matched to its nearest donor in the candidate's
%   standardized coordinates; the donor's outcomes are transferred; the
%   squared error against the receiver's observed outcomes is averaged and
%   summed over outcomes, in units of outcome standard deviation.
Zd = (Xd(:,b) - std_.mx(b)) ./ std_.sx(b);
Zq = (Xq(:,b) - std_.mx(b)) ./ std_.sx(b);
j  = knnsearch(Zd, Zq, 'K', 1);
dy = (Yq - Yd(j,:)) ./ std_.sy;
r  = sum(mean(dy.^2, 1));
end

% =======================================================================
function [X, Y] = drawPop(n, g, rs)
%DRAWPOP  Draw N individuals from population G ('A' or 'B').
if strcmp(g,'A'), [X,Y] = drawA(n, rs); else, [X,Y] = drawB(n, rs); end
end

% -----------------------------------------------------------------------
function [X, Y] = drawA(n, rs)
%DRAWA  The destructive-testing population of Section 9.11.
sig = [0.10 0.08 0.12 0.05 0.14 0.09];
Z = rand(rs, n, 2);  z1 = Z(:,1);  z2 = Z(:,2);
f = [ z1, ...
      0.5*(z1+z2) + 0.3*sin(2*pi*z1), ...
      z2, ...
      0.5 + 0.5*sin(2*pi*z1).*cos(2*pi*z2), ...
      z1.*z2, ...
      z2 + 0.2*cos(2*pi*z1) ];
X = f + randn(rs, n, numel(sig)) .* sig;
g = [ 1.2*z1 + 0.6*sin(2*pi*z2), ...
      0.8*z2.^2 + 0.4*z1, ...
      0.9*sin(2*pi*(z1+z2))/2 + 0.5*z2, ...
      1.0*z1.*z2 + 0.3*z1, ...
      0.7*cos(2*pi*z1) + 0.9*z2, ...
      0.6*z1 + 0.6*z2 + 0.5*sin(2*pi*z1.*z2) ];
k   = size(g,2);
amp = [0.22 0.28 0.25 0.30 0.24 0.26];
sd  = amp .* (1 + 0.5*sin(2*pi*z1 + (0:k-1)*pi/3));
nu  = 8;
W   = zeros(n,k);
for i = 1:nu
    W = W + randn(rs, n, k).^2;
end
Tn = randn(rs, n, k) ./ sqrt(W/nu);
Y  = g + sd .* (Tn / sqrt(nu/(nu-2)));
end

% -----------------------------------------------------------------------
function [X, Y] = drawB(n, rs)
%DRAWB  A population built to be unlike A in every respect that matters:
%   bimodal latent state, non-monotone measurements, asymmetric measurement
%   error, and conditional outcome laws that are mixtures with state-dependent
%   branch probability, hence bimodal, with heteroscedastic Laplace noise.

% Latent state: one coordinate with two separated modes, one uniform.
S  = rand(rs, n, 1) < 0.45;
z1 = min(max(0.26 + 0.48*double(S) + 0.11*randn(rs, n, 1), 0), 1);
z2 = rand(rs, n, 1);

% Measurements: non-monotone and partly non-injective, so that no single
% reading identifies the state and the useful subsets are not nested.
f = [ tanh(3*(z1 - 0.5)), ...
      abs(z2 - 0.5), ...
      z1 + z2 - 2*z1.*z2, ...
      exp(-4*(z1 - 0.3).^2), ...
      sign(z2 - 0.5) .* sqrt(z1), ...
      z1.^2 + 0.5*z2 ];

% Measurement error: asymmetric and multiplicative in scale, unlike the
% additive Gaussian error of A.  Centred so that it adds no systematic shift.
sigX = [0.12 0.07 0.10 0.06 0.13 0.09];
E    = exp(0.6*randn(rs, n, numel(sigX))) - exp(0.18);
X    = f + E .* (sigX / 0.7862);          % 0.7862 = sd of the centred factor

% Outcomes: two per destructive test.  The conditional law given the latent
% state is a two-component mixture whose branch probability moves with the
% state, so the conditional distributions are bimodal, not merely skewed.
m1 = [ 1.0*z2 + 0.5*sin(2*pi*z1), ...
       0.9*z1.*z2, ...
       0.8*abs(z1 - z2), ...
       1.1*z2.^2, ...
       0.7*cos(2*pi*z2) + 0.4*z1, ...
       0.5*z1 + 0.9*z2 ];
gap = [0.85 0.70 0.95 0.60 0.80 0.75];
m2  = m1 - gap .* (0.6 + 0.8*z2);
w   = 1 ./ (1 + exp(-6*(z1 - 0.5) + 2*(z2 - 0.5)));   % branch probability
k   = size(m1,2);
J   = rand(rs, n, k) < w;
mu  = m2;  mu(J) = m1(J);

% Heteroscedastic Laplace innovations, built from exponentials so that the
% whole experiment runs off the stream it was given.
amp = [0.18 0.22 0.20 0.24 0.19 0.21];
sd  = amp .* (0.5 + 1.5*z2);
U   = rand(rs, n, k) - 0.5;
L   = -sign(U) .* log(1 - 2*abs(U)) / sqrt(2);        % unit variance
Y   = mu + sd .* L;
end
