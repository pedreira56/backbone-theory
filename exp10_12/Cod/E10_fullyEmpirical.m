function T = E10_fullyEmpirical(reps, NL, hint)
%E10_FULLYEMPIRICAL  Backbone selection with EVERY component of the criterion
%                    estimated from a finite pilot, including the geometry.
%
%   T = E10_FULLYEMPIRICAL() finds E6.mat (see BBFINDMAT) and, for a grid of
%   pilot sizes, runs the whole selection pipeline on pilot data alone: the
%   conditional quantities AND the pairing term are estimated, no population
%   quantity is consulted.  The selected backbone is then scored against the
%   true population criterion.  Writes E10.mat and E10_table.txt next to the
%   .mat.
%   T = E10_FULLYEMPIRICAL(REPS) uses REPS replications per pilot size
%   (default 50).  T = E10_FULLYEMPIRICAL(REPS,NL) sets the planned donor
%   block size that the design is being chosen for (default: read from
%   E6.cfg, else 1000).  T = E10_FULLYEMPIRICAL(REPS,NL,FOLDER) searches
%   FOLDER first.
%
%   Needs only this file, bbFindMat.m and E6.mat.
%
%   WHAT THIS ADDS TO E9
%   E9 held the pairing term Delta(B) = E||R||^2 at its population value, so
%   only the conditional quantities were estimated.  A reviewer rightly
%   objected that the criterion depends on the INTERACTION between
%   conditional quality, anchoring regularity and the geometry of the
%   donor-receiver neighbourhoods, and that holding the geometric part
%   oracular validates only part of the pipeline.  Here nothing is oracular.
%
%   HOW THE PAIRING TERM IS ESTIMATED WITHOUT POPULATION INFORMATION
%   Delta(B) is the mean squared distance from a receiver to its nearest
%   donor among NL donors.  The pilot is far smaller than NL, and resampling
%   it with replacement is useless: repeated points give zero distances.  So
%   the estimate EXTRAPOLATES ALONG THE RATE OF THEOREM T1.  The last quarter
%   of the pilot is held out as receivers; prefixes of the rest, of three
%   sizes, serve as donors; the mean squared nearest-neighbour distance is
%   measured at each size; a line is fitted in log-log and extrapolated to
%   NL.  The fitted slope estimates -2/d_B, so the geometric dimension is
%   estimated in passing and never read off the population.  Receivers are
%   disjoint from donors by construction: without that, the largest donor set
%   contains the receivers, the distances collapse to zero and the fit is
%   destroyed.
%
%   KNOWN WEAKNESS, WORTH WATCHING IN THE OUTPUT
%   In a prototype on an analogous population the extrapolation was within
%   15% of truth for d_B >= 3 even with a pilot of 50, but UNDERESTIMATED
%   Delta by an order of magnitude for d_B = 1, where the rate is n^-2 and
%   the extrapolation amplifies the fitting error.  Underestimating Delta
%   flatters low-dimensional candidates.  The table therefore reports the
%   ratio of estimated to population Delta along the greedy path, so the bias
%   is visible and not merely suspected.
%
%   GREEDY SEARCH, NOT ENUMERATION
%   With every component estimated, each candidate costs three
%   nearest-neighbour searches, so enumerating 16384 of them per replication
%   is out of reach.  The search is forward greedy: at each step add the
%   attribute that most reduces the estimated criterion, giving a path from
%   cardinality zero to p and about p(p+1)/2 evaluations.  This is also what
%   an analyst would do in a large universe, and Section 10.5 already lists
%   search procedures as an open direction.  A failure to recover the
%   population optimum therefore mixes estimation error with greedy
%   suboptimality; the column rec_path separates them, by reporting whether
%   the population optimum lies on the greedy path at all.
%
%   WHAT IT REPORTS, the five items the review asked for
%     bsel      median selected cardinality, per cost
%     pb        fraction of replications recovering the population-optimal
%               cardinality, per cost
%     hit       fraction recovering the population-optimal backbone itself
%     excess    mean excess population criterion of the selected design
%     rec_path  fraction of replications whose greedy path contains the
%               population optimum (isolates search from estimation)
%     dratio    median ratio of estimated to population pairing term, at the
%               latent cardinality

if nargin < 1 || isempty(reps), reps = 50; end
if nargin < 2, NL = []; end
if nargin < 3, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E6', hint);
E6 = load(matFile);

S = E6.Sigma;  p = size(S,1);  d = E6.cfg.d;
card = E6.card(:);  Ubar = E6.Ubar(:);  Sres = E6.Sres(:);
subsets = E6.subsets;  nC = numel(subsets);

% planned donor block size
if isempty(NL)
    for f = {'nl','NL','nDonor','nDon','n','nBlock','nm'}
        if isfield(E6.cfg, f{1}) && isscalar(E6.cfg.(f{1}))
            NL = double(E6.cfg.(f{1}));  break
        end
    end
end
if isempty(NL), NL = 1000;
    warning('E10: planned donor block size not found in E6.cfg; using NL = 1000. Pass it explicitly if that is wrong.');
end
fprintf('E10: p = %d, d = %d, planned donor block size NL = %d, %d replications\n', p, d, NL, reps);

% ---- population reference, and Delta(B) recovered from the enumeration ----
key = zeros(nC,1);  Bidx = cell(nC,1);  Ridx = cell(nC,1);  L2 = zeros(nC,1);
for i = 1:nC
    B = sort(double(subsets{i}(:)).');
    Bidx{i} = B;  Ridx{i} = setdiff(1:p, B);
    key(i) = sum(2.^(B-1));
    if ~isempty(B) && ~isempty(Ridx{i})
        be = S(B,B) \ S(B,Ridx{i});  L2(i) = sum(be(:).^2);
    end
end
Dpop = zeros(nC,1);  ok = L2 > 0;
Dpop(ok) = (Ubar(ok) - Sres(ok)) ./ L2(ok);
if any(Dpop(ok) < -1e-9)
    error('E10: negative pairing term recovered; the criterion in E6.mat is not Sres + Delta*sum(L^2).');
end
Dpop = max(Dpop, 0);
i3 = find(card == d);  [Ustar, k3] = min(Ubar(i3));
Bstar = Bidx{i3(k3)};  kstar = sum(2.^(Bstar-1));
costs = [0.40 0.50];  bstar = zeros(1,2);
for c = 1:2, [~,kc] = min(Ubar + costs(c)*card);  bstar(c) = card(kc); end
fprintf('  population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
fprintf('  population-optimal cardinality: %d (cost 0.40), %d (cost 0.50)\n', bstar(1), bstar(2));

Cf = chol(S,'lower');
ns = [50 100 200 500 1000 5000];
nn = numel(ns);
HIT = zeros(nn,1);  EXC = zeros(nn,1);  RPATH = zeros(nn,1);  DRAT = zeros(nn,1);
BSEL = zeros(nn,2);  PB = zeros(nn,2);
t0 = tic;  first = true;

for iN = 1:nn
    n = ns(iN);
    nr = max(8, min(200, floor(n/4)));          % held-out receivers
    npool = n - nr;
    fr = [0.25 0.5 1.0];
    hit = 0;  exc = 0;  rpath = 0;  drat = nan(reps,1);
    bs = nan(reps,2);
    for r = 1:reps
        rng(2e6 + 1000*iN + r, 'twister');
        X = randn(n,p) * Cf.';
        Sg = cov(X,1);  dgS = diag(Sg);
        rec = X(npool+1:end, :);  pool = X(1:npool, :);

        B = [];  path = cell(p+1,1);  Uh = nan(p+1,1);  Dh = nan(p+1,1);
        path{1} = [];  Uh(1) = sum(dgS);           % no backbone
        for b = 1:p
            R0 = setdiff(1:p, B);
            bestU = inf;  bestA = R0(1);  bestD = NaN;
            for a = R0
                Bc = sort([B a]);  Rc = setdiff(1:p, Bc);
                if isempty(Rc), Uc = 0;  Dc = NaN;
                else
                    be  = Sg(Bc,Bc) \ Sg(Bc,Rc);
                    s2  = dgS(Rc).' - sum(Sg(Bc,Rc).*be, 1);
                    l2  = sum(be(:).^2);
                    % --- pairing term from the pilot, by extrapolation ---
                    xs = zeros(1,3);  ys = zeros(1,3);
                    for q = 1:3
                        m  = max(6, round(npool*fr(q)));
                        Dn = sum(rec(:,Bc).^2,2) + sum(pool(1:m,Bc).^2,2).' ...
                             - 2*rec(:,Bc)*pool(1:m,Bc).';
                        xs(q) = log(m);  ys(q) = log(max(mean(max(min(Dn,[],2),0)), 1e-14));
                    end
                    sl = (xs - mean(xs)) * (ys - mean(ys)).' / sum((xs - mean(xs)).^2);
                    ic = mean(ys) - sl*mean(xs);
                    Dc = exp(ic) * NL^sl;
                    Uc = sum(max(s2,0)) + Dc*l2;
                end
                if Uc < bestU, bestU = Uc;  bestA = a;  bestD = Dc; end
            end
            B = sort([B bestA]);
            path{b+1} = B;  Uh(b+1) = bestU;  Dh(b+1) = bestD;
        end

        % population criterion along the path, and selection under cost
        Upop_path = nan(p+1,1);
        for b = 0:p
            kk = find(key == sum(2.^(path{b+1}-1)), 1);
            if ~isempty(kk), Upop_path(b+1) = Ubar(kk); end
        end
        hit  = hit  + isequal(path{d+1}, Bstar);
        if ~isnan(Upop_path(d+1)), exc = exc + Upop_path(d+1)/Ustar - 1; end
        rpath = rpath + any(cellfun(@(B2) isequal(B2,Bstar), path));
        if ~isnan(Dh(d+1)) && Dpop(find(key==kstar,1)) > 0
            drat(r) = Dh(d+1) / Dpop(find(key==kstar,1));
        end
        for c = 1:2
            [~, kb] = min(Uh + costs(c)*(0:p).');
            bs(r,c) = kb - 1;
        end
        if first
            first = false;  el = toc(t0);
            fprintf('  first replication took %.2f s; projected total %.0f min\n', el, el*reps*nn/60);
        end
    end
    HIT(iN) = hit/reps;  EXC(iN) = exc/reps;  RPATH(iN) = rpath/reps;
    DRAT(iN) = median(drat(~isnan(drat)));
    for c = 1:2
        BSEL(iN,c) = median(bs(:,c));  PB(iN,c) = mean(bs(:,c) == bstar(c));
    end
    fprintf('  n = %5d done (%.1f min elapsed)\n', n, toc(t0)/60);
end

n = ns(:);
T = table(n, HIT, EXC, RPATH, DRAT, BSEL(:,1), PB(:,1), BSEL(:,2), PB(:,2), ...
    'VariableNames', {'n','hit_backbone','excess_criterion','opt_on_path', ...
                      'delta_ratio','bsel_c40','pb_c40','bsel_c50','pb_c50'});
disp(' ');
disp('=== E10: fully empirical selection (geometry estimated too) ===');
disp(T);

save(fullfile(outDir,'E10.mat'),'T','ns','reps','NL','costs','bstar','Bstar','Ustar','-v7');
fid = fopen(fullfile(outDir,'E10_table.txt'),'w');
fprintf(fid,'E10 - fully empirical backbone selection: every criterion component estimated\n\n');
fprintf(fid,'population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
fprintf(fid,'planned donor block size used for extrapolation: NL = %d\n', NL);
fprintf(fid,'population-optimal cardinality: %d (cost 0.40), %d (cost 0.50)\n', bstar(1), bstar(2));
fprintf(fid,'replications per pilot size: %d\n', reps);
fprintf(fid,'search: forward greedy; opt_on_path says whether the population optimum\n');
fprintf(fid,'was reachable by the greedy path at all, separating search from estimation\n');
fprintf(fid,'delta_ratio: estimated / population pairing term at the latent cardinality\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E10.mat and E10_table.txt to:\n  %s\n', outDir);
end
