function T = postE12_dispersao(hint)
%POSTE12_DISPERSAO  Dispersion of the excess criterion across replications, E12.
%
%   T = POSTE12_DISPERSAO() finds E12.mat (see BBFINDMAT) and reports, for each
%   population, each pilot size and each price, the mean, the median, the
%   interquartile range and the 5th and 95th percentiles of the excess
%   population criterion over the replications.  Writes E12_dispersao.txt next
%   to the .mat.  T = POSTE12_DISPERSAO(FOLDER) searches FOLDER first.
%
%   No simulation is repeated: every number comes from the per-replication
%   vectors that E12 stored in res.A.excR and res.B.excR.  Runs in a second.
%
%   WHY IT EXISTS
%   A mean excess says what is lost on average and hides whether the loss is
%   uniformly small or small in most replications and large in a few.  In a
%   near-tie regime, where exact recovery of the optimal subset is unstable
%   while the loss is negligible, the dispersion is the quantity that settles
%   whether instability of the argmin has any consequence.

if nargin < 1, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E12', hint);
E12 = load(matFile);
res = E12.res;  ns = E12.ns(:);

gens = {'A','B'};
rows = {};
for ig = 1:2
    g = gens{ig};  r = res.(g);
    nn = numel(ns);  nL = numel(r.lams);
    for iL = 1:nL
        for iN = 1:nn
            v = r.excR(:,iN,iL);  v = v(~isnan(v));
            q = quantile(v, [0.05 0.25 0.50 0.75 0.95]);
            rows(end+1,:) = {g, r.lams(iL), ns(iN), r.nD(iN), numel(v), ...
                             mean(v), q(3), q(2), q(4), q(1), q(5)};          %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'pop','preco','n','nD','reps','media','mediana','q25','q75','p05','p95'});
disp(' ');
disp('=== E12: dispersion of the excess criterion across replications ===');
disp(T);

fid = fopen(fullfile(outDir,'E12_dispersao.txt'),'w');
fprintf(fid,'E12 - dispersion of the excess criterion across replications\n\n');
fprintf(fid,'no simulation repeated: read from res.A.excR and res.B.excR\n');
fprintf(fid,'excess = population criterion of the selected design over the optimum, minus one\n');
fprintf(fid,'q25/q75 = interquartile range; p05/p95 = 90%% interval\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E12_dispersao.txt to:\n  %s\n', outDir);
end
