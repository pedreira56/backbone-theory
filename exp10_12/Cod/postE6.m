function T = postE6(hint)
%POSTE6  Oracle and implementable orderings of the candidates of E6.
%
%   T = POSTE6() finds E6.mat (see BBFINDMAT), re-evaluates the selection
%   criterion with weight two on the residual conditional variance and
%   compares the resulting ordering of the candidates with the original one,
%   writing E6_orderings.txt next to the .mat.  The enumeration is not
%   repeated.  T = POSTE6(FOLDER) searches FOLDER first.
%
%   Needs only this file and bbFindMat.m.
%
%   THE TWO CRITERIA
%   The enumeration stored, for every candidate backbone B,
%       Sres(B)   the summed residual conditional variance   sum_a s_a^2(B)
%       Ubar(B)   Sres(B) + local term                       (oracle, mean)
%   The implementable rule transfers an observed value rather than the
%   conditional mean, and its risk carries the conditional variance twice
%   (Theorem T2).  Its criterion is therefore
%       Uval(B) = 2*Sres(B) + local term = Ubar(B) + Sres(B).
%   Section 8.1 notes that the two can order candidates differently.  This
%   function measures by how much.

if nargin < 1, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E6', hint);
E6 = load(matFile);

card = E6.card(:);  Ubar = E6.Ubar(:);  Uval = Ubar + E6.Sres(:);
bMax = E6.cfg.bMax;  d = E6.cfg.d;  p = E6.cfg.p;

b = (1:max(1,bMax-1)).';  nb = numel(b);
count = zeros(nb,1); sameMin = false(nb,1); topK = zeros(nb,1);
rho = nan(nb,1); phiO = nan(nb,1); phiV = nan(nb,1);
argO = strings(nb,1); argV = strings(nb,1);
for i = 1:nb
    idx = find(card == b(i));
    count(i) = numel(idx);
    [vo,oo] = sort(Ubar(idx));  [vv,ov] = sort(Uval(idx));
    phiO(i) = vo(1);  phiV(i) = vv(1);
    sameMin(i) = oo(1) == ov(1);
    m = min(8,numel(idx));
    topK(i) = numel(intersect(oo(1:m),ov(1:m)));
    if numel(idx) > 2
        rho(i) = corr(Ubar(idx),Uval(idx),'Type','Spearman');
    end
    argO(i) = string(num2str(E6.subsets{idx(oo(1))}));
    argV(i) = string(num2str(E6.subsets{idx(ov(1))}));
end
T = table(b,count,phiO,phiV,sameMin,topK,rho,argO,argV, ...
    'VariableNames',{'b','count','Phi_oracle','Phi_value','same_minimiser', ...
                     'top8_in_common','spearman','argmin_oracle','argmin_value'});

% Cost per attribute above which the latent cardinality beats the full
% universe.  Available only when the enumeration ran to the whole universe.
if any(card == p)
    thr  = @(U) (min(U(card == d)) - min(U(card == p)))/(p - d);
    thrO = thr(Ubar);  thrV = thr(Uval);
else
    thrO = NaN;  thrV = NaN;
end

% Rank of the intended backbone at the latent cardinality, value criterion.
idx = find(card == d);  [~,ov] = sort(Uval(idx));  rankTrue = NaN;
for i = 1:numel(ov)
    if isequal(sort(double(E6.subsets{idx(ov(i))})), 1:d), rankTrue = i; break, end
end

disp(' ');
disp('=== E6: oracle ordering against the implementable ordering ===');
disp(T);
fprintf('  intended backbone, value criterion: rank %d of %d at b = %d\n', ...
        rankTrue,numel(idx),d);
fprintf('  cost threshold, oracle criterion : %.3f\n',thrO);
fprintf('  cost threshold, value criterion  : %.3f\n',thrV);

fid = fopen(fullfile(outDir,'E6_orderings.txt'),'w');
fprintf(fid,'E6 - oracle ordering against the implementable ordering\n\n');
fprintf(fid,'%s\n',evalc('disp(T)'));
fprintf(fid,'  intended backbone, value criterion: rank %d of %d at b = %d\n', ...
        rankTrue,numel(idx),d);
fprintf(fid,'  cost threshold, oracle criterion : %.3f\n',thrO);
fprintf(fid,'  cost threshold, value criterion  : %.3f\n',thrV);
fclose(fid);
fprintf('\nWrote E6_orderings.txt to:\n  %s\n',outDir);
end
