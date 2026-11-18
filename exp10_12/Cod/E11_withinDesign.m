function T = E11_withinDesign(reps, B0, hint)
%E11_WITHINDESIGN  Backbone selection using ONLY data the block design itself
%                  produces: candidates are subsets of the realized backbone.
%
%   T = E11_WITHINDESIGN() finds E6.mat (see BBFINDMAT) and simulates the
%   paper's own design: a realized backbone B0 measured in every block, two
%   blocks with exclusive attributes of their own, n individuals per block,
%   and NO individual measured in everything.  Candidates are the subsets of
%   B0.  Every quantity in the criterion is estimated where the design
%   identifies it.  Writes E11.mat and E11_table.txt next to the .mat.
%   T = E11_WITHINDESIGN(REPS) sets replications per block size (default 50).
%   T = E11_WITHINDESIGN(REPS,B0) sets the realized backbone (default
%   [1 2 3 7 10 14]).  T = E11_WITHINDESIGN(REPS,B0,FOLDER) searches FOLDER.
%
%   Needs only this file, bbFindMat.m and E6.mat.
%
%   WHY THIS EXPERIMENT EXISTS
%   E9 and E10 select from a complete pilot, in which all attributes are
%   measured on the same individuals.  That is Section 8.4's SECOND regime,
%   the calibration regime, and it is the one regime the paper's central
%   design cannot produce: joint measurement is impossible under destructive
%   measurement, and where it is possible the cross-block dependence becomes
%   identified, Cond_CI becomes testable and Proposition P6 stops biting.
%   This experiment runs in the THIRD regime, candidates restricted to subsets
%   of the realized backbone, where no new conditional law is required.
%
%   WHERE EACH QUANTITY COMES FROM, SO THAT NOTHING IS HIDDEN
%     attribute exclusive to block m  estimated in block m, where it and the
%                                     backbone were measured together
%     attribute of B0 outside the     estimated on the two blocks pooled, as
%     candidate                       it is observed in both
%     pairing term Delta(B)           MEASURED DIRECTLY, not extrapolated:
%                                     mean squared distance from receivers of
%                                     one block to their nearest donor in the
%                                     other, in the candidate's coordinates,
%                                     at the real block sizes.  Section 8.4
%                                     notes this is available because the
%                                     distance uses backbone coordinates only,
%                                     observed in every block.
%   No conditional law across blocks is ever invoked, so Proposition P6 holds
%   throughout.  There is no extrapolation and no estimated dimension.
%
%   THE TARGET
%   The population optimum is taken over subsets of B0, with exact conditional
%   quantities from Sigma and the pairing term from large Monte Carlo at each
%   block size.  It may move with block size: that is Corollary C2 at work,
%   not an artefact.

if nargin < 1 || isempty(reps), reps = 50; end
if nargin < 2 || isempty(B0), B0 = [1 2 3 7 10 14]; end
if nargin < 3, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E6', hint);
E6 = load(matFile);
S = E6.Sigma;  p = size(S,1);  dS = diag(S);
B0 = sort(B0(:)).';
if max(B0) > p
    error('E11:backboneOutsideUniverse', ['The realized backbone B0 = [%s] names ' ...
        'attribute %d, but E6.mat describes a universe of %d attributes.  Pass a B0 ' ...
        'inside the universe, or rerun runE6 in full: a quick run of runE6 produces ' ...
        'a smaller universe.'], num2str(B0), max(B0), p);
end
rest = setdiff(1:p, B0);
half = floor(numel(rest)/2);
Q1 = rest(1:half);  Q2 = rest(half+1:end);
fprintf('E11: realized backbone B0 = [%s]\n', num2str(B0));
fprintf('  block 1 exclusive: [%s]   block 2 exclusive: [%s]\n', num2str(Q1), num2str(Q2));

b0 = numel(B0);  nCand = 2^b0;
cand = cell(nCand,1);  cardC = zeros(nCand,1);
for i = 1:nCand
    sel = find(bitget(i-1, 1:b0));
    cand{i} = B0(sel);  cardC(i) = numel(sel);
end

Cf = chol(S,'lower');
ns = [100 200 500 1000 2000 5000];
nn = numel(ns);
cost = 0.50;

NBIG = 120000;
rng(4242,'twister');
BIG = randn(NBIG,p) * Cf.';

HIT = zeros(nn,1); EXC = zeros(nn,1); BSEL = zeros(nn,1);
PB  = zeros(nn,1); BOPT = zeros(nn,1);
excR = zeros(reps,nn);  hitR = false(reps,nn);  bsR = zeros(reps,nn);
EXCQ = zeros(nn,3);  EXCP = zeros(nn,2);   % quartis e intervalo de 90% do excesso
t0 = tic;  first = true;

for iN = 1:nn
    n = ns(iN);

    % ---------- target: population criterion over the subsets of B0 --------
    Upop = zeros(nCand,1);
    Rtar = min(600, n);  nrep = 3;
    for i = 1:nCand
        B = cand{i};  R = setdiff(1:p, B);
        if isempty(B)
            Upop(i) = sum(dS(R));  continue
        end
        be = S(B,B) \ S(B,R);
        s2 = dS(R).' - sum(S(B,R).*be, 1);
        l2 = sum(be(:).^2);
        dd = 0;
        for q = 1:nrep
            idx = randperm(NBIG, n + Rtar);
            dd = dd + nnMeanSq(BIG(idx(n+1:end), B), BIG(idx(1:n), B));
        end
        Upop(i) = sum(max(s2,0)) + (dd/nrep)*l2;
    end
    [Ustar, kstar] = min(Upop);
    Bstar = cand{kstar};
    [~, kc] = min(Upop + cost*cardC);  BOPT(iN) = cardC(kc);
    fprintf('  n = %5d: restricted optimum [%s], U = %.4f; optimal cardinality under cost %.2f: %d\n', ...
            n, num2str(Bstar), Ustar, cost, BOPT(iN));

    % ---------- replications of the design ----------
    hit = 0;  exc = 0;  bs = zeros(reps,1);
    for r = 1:reps
        rng(3e6 + 1000*iN + r, 'twister');
        X1 = randn(n,p) * Cf.';           % block 1 observes B0 and Q1
        X2 = randn(n,p) * Cf.';           % block 2 observes B0 and Q2
        Rsub = min(500, n);
        i1 = randperm(n, Rsub);  i2 = randperm(n, Rsub);

        Uh = zeros(nCand,1);
        for i = 1:nCand
            B = cand{i};  outB = setdiff(B0, B);
            if isempty(B)
                Uh(i) = sum(var(X1(:,[Q1 outB]),1)) + sum(var(X2(:,Q2),1));
                continue
            end
            s2tot = 0;  l2tot = 0;
            % exclusive attributes, each inside its own block
            for blk = 1:2
                if blk == 1, Xb = X1;  cols = Q1;  else, Xb = X2;  cols = Q2;  end
                if isempty(cols), continue, end
                Sb  = cov(Xb(:,[B cols]), 1);  k = numel(B);  dSb = diag(Sb);
                be  = Sb(1:k,1:k) \ Sb(1:k,k+1:end);
                s2  = dSb(k+1:end).' - sum(Sb(1:k,k+1:end).*be, 1);
                s2tot = s2tot + sum(max(s2,0));  l2tot = l2tot + sum(be(:).^2);
            end
            % attributes of B0 left out of the candidate: observed in both
            if ~isempty(outB)
                Xp  = [X1(:,[B outB]); X2(:,[B outB])];
                Sp  = cov(Xp, 1);  k = numel(B);  dSp = diag(Sp);
                be  = Sp(1:k,1:k) \ Sp(1:k,k+1:end);
                s2  = dSp(k+1:end).' - sum(Sp(1:k,k+1:end).*be, 1);
                s2tot = s2tot + sum(max(s2,0));  l2tot = l2tot + sum(be(:).^2);
            end
            % pairing term, measured in both directions at the real sizes
            d12 = nnMeanSq(X1(i1,B), X2(:,B));
            d21 = nnMeanSq(X2(i2,B), X1(:,B));
            Uh(i) = s2tot + 0.5*(d12 + d21)*l2tot;
        end

        [~, kh] = min(Uh);
        hit = hit + isequal(cand{kh}, Bstar);
        exc = exc + Upop(kh)/Ustar - 1;
        [~, kb] = min(Uh + cost*cardC);
        bs(r) = cardC(kb);
        excR(r,iN) = Upop(kh)/Ustar - 1;
        hitR(r,iN) = isequal(cand{kh}, Bstar);
        bsR(r,iN)  = cardC(kb);
        if first
            first = false;  el = toc(t0);
            fprintf('  first replication took %.2f s; projected total %.0f min\n', el, el*reps*nn/60);
        end
    end
    HIT(iN) = hit/reps;  EXC(iN) = exc/reps;
    BSEL(iN) = median(bs);  PB(iN) = mean(bs == BOPT(iN));
    EXCQ(iN,:) = quantile(excR(:,iN), [0.25 0.50 0.75]);
    EXCP(iN,:) = quantile(excR(:,iN), [0.05 0.95]);
    fprintf('    excesso: média %.4f, mediana %.4f, IQ [%.4f, %.4f], 90%% [%.4f, %.4f]\n', ...
            EXC(iN), EXCQ(iN,2), EXCQ(iN,1), EXCQ(iN,3), EXCP(iN,1), EXCP(iN,2));
    fprintf('  n = %5d done (%.1f min elapsed)\n', n, toc(t0)/60);
end

n = ns(:);
T = table(n, HIT, EXC, EXCQ(:,2), EXCQ(:,1), EXCQ(:,3), EXCP(:,1), EXCP(:,2), ...
          BSEL, PB, BOPT, ...
    'VariableNames', {'n','hit_subset','excess_criterion','exc_med','exc_q25', ...
                      'exc_q75','exc_p05','exc_p95','bsel','pb','bopt'});
disp(' ');
disp('=== E11: selection within the realized design (no complete pilot) ===');
disp(T);

save(fullfile(outDir,'E11.mat'),'T','ns','reps','B0','Q1','Q2','cost', ...
     'excR','hitR','bsR','EXCQ','EXCP','-v7');
fid = fopen(fullfile(outDir,'E11_table.txt'),'w');
fprintf(fid,'E11 - selection using only data the block design produces\n\n');
fprintf(fid,'realized backbone B0 = [%s]\n', num2str(B0));
fprintf(fid,'block 1 exclusive [%s], block 2 exclusive [%s]\n', num2str(Q1), num2str(Q2));
fprintf(fid,'candidates: the %d subsets of B0; cost per attribute %.2f\n', nCand, cost);
fprintf(fid,'replications per block size: %d\n', reps);
fprintf(fid,'the pairing term is MEASURED, not extrapolated; no complete pilot is used;\n');
fprintf(fid,'no conditional law across blocks is invoked, so Proposition P6 holds throughout\n');
fprintf(fid,'bopt is the population-optimal cardinality at that block size\n\n');
fprintf(fid,'exc_med, exc_q25, exc_q75: median and interquartile range of the excess\n');
fprintf(fid,'exc_p05, exc_p95: 90%% interval of the excess across replications\n');
fprintf(fid,'per-replication vectors are stored in E11.mat as excR, hitR and bsR\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E11.mat and E11_table.txt to:\n  %s\n', outDir);
end

% ----------------------------------------------------------------------
function dm = nnMeanSq(rec, don)
%NNMEANSQ  mean squared distance from each row of REC to its nearest row of DON
dm = 0;  tot = 0;  nb = 256;
sd = sum(don.^2, 2).';
for ii = 1:nb:size(rec,1)
    blk = rec(ii:min(ii+nb-1, size(rec,1)), :);
    D = sum(blk.^2, 2) + sd - 2*blk*don.';
    dm = dm + sum(max(min(D, [], 2), 0));
    tot = tot + size(blk,1);
end
dm = dm / max(tot,1);
end
