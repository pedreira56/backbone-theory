function T = E10b_greedyOnly(reps, hint)
%E10B_GREEDYONLY  Greedy search with the pairing term held at its POPULATION
%                 value: isolates the cost of the search from the cost of
%                 estimating the geometry.
%
%   T = E10B_GREEDYONLY() finds E6.mat (see BBFINDMAT) and repeats the
%   forward-greedy selection of E10, but with Delta(B) taken from the stored
%   enumeration instead of extrapolated from the pilot.  Only the conditional
%   quantities are estimated.  Writes E10b.mat and E10b_table.txt next to the
%   .mat.  T = E10B_GREEDYONLY(REPS) sets replications per pilot size
%   (default 50, matching E10).  T = E10B_GREEDYONLY(REPS,FOLDER) searches
%   FOLDER first.
%
%   Needs only this file, bbFindMat.m and E6.mat.  Runs in seconds: there are
%   no nearest-neighbour searches.
%
%   WHY IT EXISTS
%   Three experiments now differ in two dimensions:
%     E9    exhaustive enumeration, population Delta  -> recovery 1.00 from n = 100
%     E10b  greedy search,          population Delta  -> this script
%     E10   greedy search,          estimated Delta   -> recovery 0.54 at n = 100
%   The gap between E9 and E10b is the price of the search; the gap between
%   E10b and E10 is the price of estimating the geometry.  Without this middle
%   row the two are confounded, and the remedy differs: a better search in one
%   case, a larger pilot in the other.
%
%   The greedy path, the cost grid, the pilot sizes, the seeds and the
%   in-sample estimation of the conditional quantities are exactly those of
%   E10, so the rows are comparable line by line.

if nargin < 1 || isempty(reps), reps = 50; end
if nargin < 2, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E6', hint);
E6 = load(matFile);

S = E6.Sigma;  p = size(S,1);  d = E6.cfg.d;
card = E6.card(:);  Ubar = E6.Ubar(:);  Sres = E6.Sres(:);
subsets = E6.subsets;  nC = numel(subsets);

key = zeros(nC,1);  L2 = zeros(nC,1);  Bidx = cell(nC,1);
for i = 1:nC
    B = sort(double(subsets{i}(:)).');
    Bidx{i} = B;  key(i) = sum(2.^(B-1));
    R = setdiff(1:p, B);
    if ~isempty(B) && ~isempty(R)
        be = S(B,B) \ S(B,R);  L2(i) = sum(be(:).^2);
    end
end
Dpop = zeros(nC,1);  ok = L2 > 0;
Dpop(ok) = (Ubar(ok) - Sres(ok)) ./ L2(ok);
if any(Dpop(ok) < -1e-9)
    error('E10b: negative pairing term recovered; the criterion in E6.mat is not Sres + Delta*sum(L^2).');
end
Dpop = max(Dpop, 0);
Dof = containers.Map(num2cell(key), num2cell(Dpop));   % Delta by subset key

i3 = find(card == d);  [Ustar, k3] = min(Ubar(i3));
Bstar = Bidx{i3(k3)};
costs = [0.40 0.50];  bstar = zeros(1,2);
for c = 1:2, [~,kc] = min(Ubar + costs(c)*card);  bstar(c) = card(kc); end
fprintf('E10b: population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
fprintf('  population-optimal cardinality: %d (cost 0.40), %d (cost 0.50)\n', bstar(1), bstar(2));

Cf = chol(S,'lower');
ns = [50 100 200 500 1000 5000];
nn = numel(ns);
HIT = zeros(nn,1);  EXC = zeros(nn,1);  BSEL = zeros(nn,2);  PB = zeros(nn,2);
excR = nan(reps,nn);  hitR = false(reps,nn);  bsR = nan(reps,nn,2);
EXCQ = zeros(nn,3);  EXCP = zeros(nn,2);   % quartis e intervalo de 90% do excesso
t0 = tic;

for iN = 1:nn
    n = ns(iN);
    hit = 0;  exc = 0;  bs = nan(reps,2);
    for r = 1:reps
        rng(2e6 + 1000*iN + r, 'twister');       % same seeds as E10
        X  = randn(n,p) * Cf.';
        Sg = cov(X,1);  dgS = diag(Sg);

        B = [];  path = cell(p+1,1);  Uh = nan(p+1,1);
        path{1} = [];  Uh(1) = sum(dgS);
        for b = 1:p
            R0 = setdiff(1:p, B);
            bestU = inf;  bestA = R0(1);
            for a = R0
                Bc = sort([B a]);  Rc = setdiff(1:p, Bc);
                if isempty(Rc)
                    Uc = 0;
                else
                    be = Sg(Bc,Bc) \ Sg(Bc,Rc);
                    s2 = dgS(Rc).' - sum(Sg(Bc,Rc).*be, 1);
                    kc = sum(2.^(Bc-1));
                    if ~isKey(Dof, kc)
                        error('E10b:missingSubset', ['E6.mat does not contain the subset ' ...
                            '{%s}.  E10b needs the pairing term of every subset the greedy ' ...
                            'search visits, which requires the complete enumeration that ' ...
                            'runE6 produces in a full run; an E6.mat from a quick run of ' ...
                            'runE6 covers a smaller universe and fails here.'], num2str(Bc));
                    end
                    Uc = sum(max(s2,0)) + Dof(kc) * sum(be(:).^2);
                end
                if Uc < bestU, bestU = Uc;  bestA = a; end
            end
            B = sort([B bestA]);  path{b+1} = B;  Uh(b+1) = bestU;
        end

        hit = hit + isequal(path{d+1}, Bstar);
        hitR(r,iN) = isequal(path{d+1}, Bstar);
        kk = find(key == sum(2.^(path{d+1}-1)), 1);
        if ~isempty(kk)
            exc = exc + Ubar(kk)/Ustar - 1;
            excR(r,iN) = Ubar(kk)/Ustar - 1;
        end
        for c = 1:2
            [~, kb] = min(Uh + costs(c)*(0:p).');
            bs(r,c) = kb - 1;
            bsR(r,iN,c) = kb - 1;
        end
    end
    HIT(iN) = hit/reps;  EXC(iN) = exc/reps;
    for c = 1:2
        BSEL(iN,c) = median(bs(:,c));  PB(iN,c) = mean(bs(:,c) == bstar(c));
    end
    v = excR(:,iN);  v = v(~isnan(v));
    EXCQ(iN,:) = quantile(v, [0.25 0.50 0.75]);
    EXCP(iN,:) = quantile(v, [0.05 0.95]);
    fprintf('    excesso: média %.4f, mediana %.4f, IQ [%.4f, %.4f], 90%% [%.4f, %.4f]\n', ...
            EXC(iN), EXCQ(iN,2), EXCQ(iN,1), EXCQ(iN,3), EXCP(iN,1), EXCP(iN,2));
    fprintf('  n = %5d done (%.1f s elapsed)\n', n, toc(t0));
end

n = ns(:);
T = table(n, HIT, EXC, EXCQ(:,2), EXCQ(:,1), EXCQ(:,3), EXCP(:,1), EXCP(:,2), ...
          BSEL(:,1), PB(:,1), BSEL(:,2), PB(:,2), ...
    'VariableNames', {'n','hit_backbone','excess_criterion','exc_med','exc_q25', ...
                      'exc_q75','exc_p05','exc_p95','bsel_c40','pb_c40', ...
                      'bsel_c50','pb_c50'});
disp(' ');
disp('=== E10b: greedy search, population pairing term ===');
disp(T);

save(fullfile(outDir,'E10b.mat'),'T','ns','reps','costs','bstar','Bstar','Ustar', ...
     'excR','hitR','bsR','EXCQ','EXCP','-v7');
fid = fopen(fullfile(outDir,'E10b_table.txt'),'w');
fprintf(fid,'E10b - greedy search with the pairing term at its population value\n\n');
fprintf(fid,'population optimum at b = %d: [%s], U = %.4f\n', d, num2str(Bstar), Ustar);
fprintf(fid,'population-optimal cardinality: %d (cost 0.40), %d (cost 0.50)\n', bstar(1), bstar(2));
fprintf(fid,'replications per pilot size: %d (same seeds as E10)\n', reps);
fprintf(fid,'compare line by line with E10 (estimated Delta) and with E9 (exhaustive\n');
fprintf(fid,'enumeration, population Delta) to separate search cost from geometry cost\n\n');
fprintf(fid,'exc_med, exc_q25, exc_q75: median and interquartile range of the excess\n');
fprintf(fid,'exc_p05, exc_p95: 90%% interval of the excess across replications\n');
fprintf(fid,'per-replication vectors are stored in E10b.mat as excR, hitR and bsR\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E10b.mat and E10b_table.txt to:\n  %s\n', outDir);
end
