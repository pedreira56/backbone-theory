function T = E9_pilotSelection(reps, hint)
%E9_PILOTSELECTION  Backbone selection with the criterion ESTIMATED from a
%                   finite pilot sample (Section 8.4 meets Section 9).
%
%   T = E9_PILOTSELECTION() finds E6.mat (see BBFINDMAT), and for a grid of
%   pilot sizes draws pilot samples from the population of E6, estimates the
%   selection criterion of Section 8 from each pilot by three strategies,
%   selects a backbone with each, and scores it against the TRUE population
%   criterion.  Writes E9.mat and E9_table.txt next to the .mat.
%   T = E9_PILOTSELECTION(REPS) uses REPS replications per pilot size
%   (default 100).  T = E9_PILOTSELECTION(REPS,FOLDER) searches FOLDER first.
%
%   Needs only this file, bbFindMat.m and E6.mat.
%
%   WHY THIS EXPERIMENT EXISTS
%   Every other experiment evaluates candidates with population quantities
%   known in closed form.  That validates the criterion but says nothing
%   about selecting with data.  Section 8.4 recommends sample splitting or
%   cross-validation and warns against fitting and selecting on the same
%   data; Section 10.4 admits the effect of estimation error on the ordering
%   is not quantified.  This quantifies it.
%
%   THE THREE STRATEGIES, all over the full enumeration
%     in-sample  Sigma-hat from the whole pilot; residual conditional
%                variance is the in-sample Schur complement.  Its bias grows
%                with |B|, so it distorts comparisons ACROSS cardinalities
%                while leaving comparisons at FIXED cardinality nearly
%                unaffected.  That is the distinction this experiment isolates.
%     split      beta from the first half, residual variance measured on the
%                second half.  Honest but spends half the pilot; a candidate
%                is declared not evaluable when the half is too small for its
%                cardinality (|B| + 5 > n/2), which is itself a finding.
%     cv5        five-fold cross-validation of the residual variance: honest
%                and uses every observation.  Costs five solves per candidate.
%   Anchoring constants come from the full pilot in all three: they are norms
%   of regression coefficients, not out-of-sample risks.
%
%   TWO COSTS, ON PURPOSE
%   Cardinality is selected under a per-attribute cost.  With cost 0.50, well
%   above the threshold 0.384 of Proposition P8, the population answer is the
%   latent cardinality by a comfortable margin.  With cost 0.40, barely above
%   the threshold, the decision is marginal and estimation error should bite.
%   Reporting both separates a robust recommendation from a fragile one.
%
%   WHAT IS HELD AT ITS POPULATION VALUE, AND WHY
%   The matching term Delta(B) = E||R||^2 keeps the value the enumeration of
%   E6 established.  Estimating it too would need a nearest-neighbour
%   simulation inside each of the 16384 candidates of every replication.  The
%   error therefore enters through one channel, the conditional quantities,
%   which is the channel Section 8.4 discusses.  The rates below are
%   consequently optimistic, and the written table says so.

if nargin < 1 || isempty(reps), reps = 100; end
if nargin < 2, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E6', hint);
E6 = load(matFile);

S = E6.Sigma;  p = size(S,1);  d = E6.cfg.d;
card = E6.card(:);  Ubar = E6.Ubar(:);  Sres = E6.Sres(:);
subsets = E6.subsets;  nC = numel(subsets);
fprintf('E9: universe p = %d, latent dimension d = %d, %d candidates\n', p, d, nC);

% ---- index sets and population anchoring constants ---------------------
Bidx = cell(nC,1);  Ridx = cell(nC,1);  L2sum = zeros(nC,1);
for i = 1:nC
    B = sort(double(subsets{i}(:)).');
    Bidx{i} = B;  Ridx{i} = setdiff(1:p, B);
    if isempty(Ridx{i}) || isempty(B), continue, end
    be = S(B,B) \ S(B,Ridx{i});
    L2sum(i) = sum(be(:).^2);
end

% ---- recover Delta(B) from the stored enumeration ----------------------
%  U(B) = Sres(B) + Delta(B) * sum_a L_a(B)^2
Delta = zeros(nC,1);  ok = L2sum > 0;
Delta(ok) = (Ubar(ok) - Sres(ok)) ./ L2sum(ok);
bad = ok & (Delta < -1e-9);
if any(bad)
    error(['E9: recovered a negative matching term for %d of %d candidates. ' ...
           'The criterion stored in E6.mat is then not Sres + Delta*sum(L^2) ' ...
           'with unit weights, and this script cannot reconstruct it. ' ...
           'Send me the lines of your E6 script that build Ubar and Sres.'], ...
           sum(bad), nC);
end
Delta = max(Delta, 0);
fprintf('  matching term recovered from the enumeration:\n');
for b = 0:min(p,6)
    j = card == b;
    if any(j)
        fprintf('     b = %2d:  mean %.3e   sd %.1e   (%d candidates)\n', ...
                b, mean(Delta(j)), std(Delta(j)), sum(j));
    end
end

% ---- population reference ----------------------------------------------
i3 = find(card == d);
[Ustar, k3] = min(Ubar(i3));
Bstar = Bidx{i3(k3)};
costs = [0.40 0.50];
bstar = zeros(1,numel(costs));
for c = 1:numel(costs)
    [~, kc] = min(Ubar + costs(c)*card);  bstar(c) = card(kc);
end
fprintf('  population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
for c = 1:numel(costs)
    fprintf('  cost %.2f per attribute: population-optimal cardinality %d\n', ...
            costs(c), bstar(c));
end

Cf = chol(S, 'lower');
ns = [25 35 50 100 200 500 1000 5000];        % ns = [50 100 200 500 1000 5000];
nn = numel(ns);  nS = 3;
hit3 = zeros(nn,nS);  regret = zeros(nn,nS);
pb = zeros(nn,nS,numel(costs));  bmed = zeros(nn,nS,numel(costs));
skipped = zeros(nn,1);
t0 = tic;  firstDone = false;

for iN = 1:nn
    n = ns(iN);  h = floor(n/2);
    folds = 5;  cut = round(linspace(0, n, folds+1));
    evalSplit = true(nC,1);
    for i = 1:nC
        if numel(Bidx{i}) + 5 > h, evalSplit(i) = false; end
    end
    skipped(iN) = sum(~evalSplit);
    H = zeros(1,nS);  RG = zeros(1,nS);
    BSEL = nan(reps,nS,numel(costs));
    for r = 1:reps
        rng(1e6 + 1000*iN + r, 'twister');
        X  = randn(n,p) * Cf.';
        Sg = cov(X,1);  dgS = diag(Sg);
        S1 = cov(X(1:h,:),1);  X2 = X(h+1:end,:);
        Sf = cell(folds,1);  Ko = cell(folds,1);
        for f = 1:folds
            keep = true(n,1);  keep(cut(f)+1:cut(f+1)) = false;
            Sf{f} = cov(X(keep,:),1);  Ko{f} = find(~keep);
        end
        U = inf(nC,nS);
        for i = 1:nC
            B = Bidx{i};  R = Ridx{i};
            if isempty(R), U(i,:) = 0; continue, end
            if isempty(B)
                L2 = 0;
                s2a = dgS(R).';
                s2b = mean(X2(:,R).^2, 1);
                s2c = mean(X(:,R).^2, 1);
            else
                be  = Sg(B,B) \ Sg(B,R);  L2 = sum(be(:).^2);
                s2a = dgS(R).' - sum(Sg(B,R).*be, 1);
                if evalSplit(i)
                    b1  = S1(B,B) \ S1(B,R);
                    s2b = mean((X2(:,R) - X2(:,B)*b1).^2, 1);
                else
                    s2b = [];
                end
                s2c = zeros(1, numel(R));
                for f = 1:folds
                    bk = Sf{f}(B,B) \ Sf{f}(B,R);
                    Xo = X(Ko{f},:);
                    s2c = s2c + sum((Xo(:,R) - Xo(:,B)*bk).^2, 1);
                end
                s2c = s2c / n;
            end
            loc = Delta(i) * L2;
            U(i,1) = sum(max(s2a,0)) + loc;
            if ~isempty(s2b), U(i,2) = sum(max(s2b,0)) + loc; end
            U(i,3) = sum(max(s2c,0)) + loc;
        end
        for m = 1:nS
            [~, k] = min(U(i3,m));  Bh = Bidx{i3(k)};
            H(m)  = H(m)  + isequal(Bh, Bstar);
            RG(m) = RG(m) + Ubar(i3(k))/Ustar - 1;
            for c = 1:numel(costs)
                [~, kk] = min(U(:,m) + costs(c)*card);
                BSEL(r,m,c) = card(kk);
            end
        end
        if ~firstDone
            firstDone = true;  el = toc(t0);
            fprintf('  first replication took %.2f s; projected total %.0f min\n', ...
                    el, el*reps*nn/60);
        end
    end
    hit3(iN,:)   = H / reps;
    regret(iN,:) = RG / reps;
    for m = 1:nS
        for c = 1:numel(costs)
            pb(iN,m,c)   = mean(BSEL(:,m,c) == bstar(c));
            bmed(iN,m,c) = median(BSEL(:,m,c));
        end
    end
    fprintf('  n = %5d done (%d candidates not evaluable by splitting; %.1f min elapsed)\n', ...
            n, skipped(iN), toc(t0)/60);
end

n = ns(:);
T = table(n, hit3(:,1), hit3(:,2), hit3(:,3), ...
             regret(:,1), regret(:,2), regret(:,3), ...
             bmed(:,1,1), bmed(:,2,1), bmed(:,3,1), ...
             pb(:,1,1), pb(:,2,1), pb(:,3,1), ...
             bmed(:,1,2), bmed(:,3,2), pb(:,1,2), pb(:,3,2), skipped, ...
    'VariableNames', {'n','hit3_in','hit3_split','hit3_cv5', ...
        'reg_in','reg_split','reg_cv5', ...
        'bmed_in_c40','bmed_split_c40','bmed_cv5_c40', ...
        'pb_in_c40','pb_split_c40','pb_cv5_c40', ...
        'bmed_in_c50','bmed_cv5_c50','pb_in_c50','pb_cv5_c50','not_evaluable_split'});

disp(' ');
disp('=== E9: selection from a finite pilot sample ===');
disp(T);

save(fullfile(outDir,'E9.mat'), 'T','ns','reps','costs','bstar', ...
     'Bstar','Ustar','Delta','hit3','regret','pb','bmed','skipped','-v7');
fid = fopen(fullfile(outDir,'E9_table.txt'),'w');
fprintf(fid,'E9 - backbone selection with the criterion estimated from a finite pilot\n\n');
fprintf(fid,'population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
fprintf(fid,'cost threshold of Proposition P8 in this population: 0.384\n');
for c = 1:numel(costs)
    fprintf(fid,'cost %.2f per attribute: population-optimal cardinality %d\n', costs(c), bstar(c));
end
fprintf(fid,'replications per pilot size: %d\n', reps);
fprintf(fid,'the matching term is held at its population value, so these rates\n');
fprintf(fid,'are optimistic: an analyst would have to estimate it as well\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E9.mat and E9_table.txt to:\n  %s\n', outDir);
end
