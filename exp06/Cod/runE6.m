function out = runE6(quick, reuse)
%RUNE6  Driver for experiment E6 (Section 8, by exhaustive enumeration).
%
%   OUT = RUNE6()             full run   (use this for the paper)
%   OUT = RUNE6(true)         quick run  (smaller universe, smoke test)
%   OUT = RUNE6(false,true)   reuse a saved enumeration and redo only the
%                             post-processing.  The enumeration itself is the
%                             expensive part and does not depend on the cost
%                             grid, so exploring costs does not require
%                             recomputing it.
%
%   Layout: put every .m file in  ...\Exp6\Cod\  and results appear in
%   ...\Exp6\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs: E6.mat, E6_fig.pdf/.png, E6_table.txt
%
%   This experiment answers, for the formulation of Section 8, the question a
%   reader is entitled to ask of any optimisation problem stated in a paper:
%   does its solution behave as the formulation intends?  With a universe of
%   fourteen eligible attributes and cardinality capped at six, every
%   candidate can be evaluated, so the answer requires no algorithm and
%   inherits no approximation error from one.

if nargin < 1 || isempty(quick), quick = false; end
if nargin < 2 || isempty(reuse), reuse = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E6] code   : %s\n[E6] results: %s\n\n', codeDir, outDir);

cfg = struct();
if quick
    cfg.p = 10; cfg.bMax = 4; cfg.n = 1500; cfg.nRecv = 200; cfg.repDelta = 2;
else
    % bMax = p: the enumeration runs to the full universe.  Stopping short
    % of it hides the degeneracy of Section 8.3, which appears only as the
    % backbone approaches the whole set of eligible attributes: each added
    % coordinate removes a term from the sum, and at b = p the sum is empty
    % and the criterion is identically zero.  With b capped well below p the
    % criterion has an interior minimum even without costs, and the role of
    % the costs -- to make the problem well posed without an artificial
    % cardinality constraint -- cannot be seen at all.
    cfg.p = 14; cfg.bMax = 14; cfg.n = 3000; cfg.nRecv = 400; cfg.repDelta = 3;
end

% Cost grid.  A cost that is too small does not remove the degeneracy: going
% from the latent cardinality to the full universe gains the whole criterion
% and costs the number of added attributes times the unit cost, so the unit
% cost must exceed the ratio of the two for an interior optimum to exist.
% The grid brackets that threshold deliberately.
cfg.costs = [0 0.25 0.50];

matFile = fullfile(outDir,'E6.mat');
if reuse && isfile(matFile)
    fprintf('[E6] reusing the saved enumeration in %s\n',matFile);
    out = load(matFile);
    out.costs = cfg.costs;
    nb = out.cfg.bMax + 1;
    out.Utot = nan(nb,numel(cfg.costs));
    out.bestBcost = cell(nb,numel(cfg.costs));
    for b = 0:out.cfg.bMax
        idx = find(out.card == b);
        for c = 1:numel(cfg.costs)
            [out.Utot(b+1,c),k2] = min(out.Ubar(idx) + cfg.costs(c)*b);
            out.bestBcost{b+1,c}  = out.subsets{idx(k2)};
        end
    end
else
    out = E6_enumeration(cfg);
end
out.outDir = outDir;

% Where the optimum sits as a function of the unit cost: the transition from
% the degenerate solution to the latent cardinality, and the threshold at
% which it happens.
cGrid = linspace(0,1.2,121);
bStar = nan(size(cGrid));
for i = 1:numel(cGrid)
    best = inf;
    for b = 0:out.cfg.bMax
        v = min(out.Ubar(out.card == b)) + cGrid(i)*b;
        if v < best, best = v; bStar(i) = b; end
    end
end
out.costGrid = cGrid; out.bStar = bStar;

d = out.cfg.d;
fprintf('\n=== E6: ranking at the latent cardinality (b = %d) ===\n',d);
for i = 1:min(5,numel(out.rank3))
    fprintf('  %2d.  %8.4f   {%s}\n',i,out.val3(i),num2str(out.rank3{i}));
end
fprintf('  ...\n  last %8.4f   {%s}\n',out.val3(end),num2str(out.rank3{end}));
verdict = {'is NOT','IS'};
fprintf('  the intended backbone {%s} %s the minimiser  (rank %d of %d)\n', ...
        num2str(1:d), verdict{out.trueIsBest + 1}, out.rankTrue, numel(out.val3));

disp(' ');
disp('=== E6: dispersion among backbones of equal cardinality ===');
disp(out.dispersion);

fprintf('=== E6: best backbone at each cardinality ===\n');
fprintf('  %3s %12s %12s   %s\n','b','Phi(b)','Phi_loc(b)','argmin');
for b = 0:out.cfg.bMax
    fprintf('  %3d %12.4f %12.6f   {%s}\n', ...
        b,out.Phi(b+1),out.PhiLoc(b+1),num2str(out.bestB{b+1}));
end

fprintf('\n=== E6: per-attribute costs restore a well-posed problem ===\n');
fprintf('  %3s',' b');
for c = 1:numel(out.costs), fprintf(' %14s',sprintf('cost %.2f',out.costs(c))); end
fprintf('\n');
for b = 0:out.cfg.bMax
    fprintf('  %3d',b);
    for c = 1:numel(out.costs), fprintf(' %14.4f',out.Utot(b+1,c)); end
    fprintf('\n');
end
for c = 1:numel(out.costs)
    [~,k] = min(out.Utot(2:end,c));
    fprintf('  cost %.2f -> optimal b = %d, backbone {%s}\n', ...
        out.costs(c), k, num2str(out.bestBcost{k+1,c}));
end
kThr = find(out.bStar <= out.cfg.d, 1);
if ~isempty(kThr)
    fprintf('  the optimum first reaches the latent cardinality at cost %.3f\n', ...
            out.costGrid(kThr));
end

% ---- figure -----------------------------------------------------------
f = figure('Position',[100 100 1300 380]);

% (a) the distribution of the criterion at the latent cardinality
subplot(1,3,1); hold on; box on
v = out.Ubar(out.card == d);
histogram(v,30,'FaceColor',[0.65 0.72 0.82],'EdgeColor',[.35 .35 .35], ...
          'DisplayName',sprintf('all %d backbones with b = %d',numel(v),d));
yl = ylim;
plot(out.valTrue*[1 1], yl,'-','Color',[0.85 0.33 0.10],'LineWidth',1.8, ...
     'DisplayName',sprintf('intended backbone (rank %d of %d)', ...
                   out.rankTrue,numel(v)));
ylim(yl);
set(gca,'FontSize',9);
xlabel('criterion value');
ylabel('number of backbones');
text(-0.02,1.05,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','northeast','FontSize',8,'Box','off');

% (b) the criterion against cardinality, with and without costs
subplot(1,3,2); hold on; box on
% b = 0 is the criterion with no backbone at all, equal to the total
% variance; it exceeds every other value by enough to compress the range
% where the comparison actually happens.  It is reported in the table and
% omitted from the plot.
bs   = 1:out.cfg.bMax;
cols = [0.25 0.25 0.25; 0.00 0.45 0.74; 0.85 0.33 0.10];
for c = 1:numel(out.costs)
    plot(bs, out.Utot(2:end,c),'o-','Color',cols(min(c,3),:),'MarkerSize',5, ...
         'LineWidth',1.3,'MarkerFaceColor','w', ...
         'DisplayName',sprintf('cost per attribute = %.2f',out.costs(c)));
end
for c = 1:numel(out.costs)
    [mv,k] = min(out.Utot(2:end,c));
    plot(bs(k), mv,'p','Color',cols(min(c,3),:),'MarkerSize',12, ...
         'MarkerFaceColor',cols(min(c,3),:),'HandleVisibility','off');
end
set(gca,'FontSize',9,'YScale','log');
xlabel('b   (backbone cardinality)');
ylabel('best value of the criterion');
xlim([0.7, out.cfg.bMax+0.3]);
% The no-cost curve falls towards zero as the backbone absorbs the universe,
% so a linear ordinate would compress everything else against the axis.
yy = out.Utot(2:end,:); yy = yy(yy > 0);
ylim([min(yy)/2, max(out.Utot(2:end,:),[],'all')*1.5]);
text(-0.02,1.05,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','northeast','FontSize',8,'Box','off');

% ---- (c) where the optimum sits as the unit cost varies ---------------
subplot(1,3,3); hold on; box on
stairs(out.costGrid, out.bStar,'-','Color',[0.25 0.25 0.25],'LineWidth',1.4, ...
       'DisplayName','optimal cardinality');
yline(out.cfg.d,'--','Color',[0.85 0.33 0.10],'LineWidth',1.3, ...
      'DisplayName',sprintf('latent dimension (%d)',out.cfg.d));
yline(out.cfg.p,':','Color',[0.00 0.45 0.74],'LineWidth',1.3, ...
      'DisplayName',sprintf('full universe (%d)',out.cfg.p));
set(gca,'FontSize',9);
xlabel('cost per attribute');
ylabel('optimal backbone cardinality');
ylim([-0.5, out.cfg.p+0.5]);
text(-0.02,1.05,'(c)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','east','FontSize',8,'Box','off');

% A destination that is open in a viewer, or being synchronised by a cloud
% client, refuses to be overwritten on Windows.  Losing a completed run to
% that is avoidable: report the obstacle, fall back to a timestamped name,
% and carry on to save the results.
saveFigure(f, outDir, 'E6_fig');

try
    save(fullfile(outDir,'E6.mat'),'-struct','out');
catch ME
    fprintf(2,'[E6] could not write E6.mat (%s); results remain in the output struct.\n',ME.identifier);
end
fid = fopen(fullfile(outDir,'E6_table.txt'),'w');
if fid < 0
    fprintf(2,'[E6] could not open E6_table.txt for writing; skipping the table file.\n');
    fprintf('\nDone (figure and .mat only).\n');
    return
end
fprintf(fid,'E6 - exhaustive enumeration of the selection criterion\n\n');
fprintf(fid,'Ranking at the latent cardinality (b = %d)\n',d);
for i = 1:min(8,numel(out.rank3))
    fprintf(fid,'  %2d. %8.4f  {%s}\n',i,out.val3(i),num2str(out.rank3{i}));
end
fprintf(fid,'  last %8.4f  {%s}\n',out.val3(end),num2str(out.rank3{end}));
fprintf(fid,'  intended backbone is the minimiser: %d\n\n',out.trueIsBest);
fprintf(fid,'Dispersion among backbones of equal cardinality\n');
fprintf(fid,'%s\n',evalc('disp(out.dispersion)'));
fprintf(fid,'\nBest backbone at each cardinality\n');
fprintf(fid,'  %3s %12s %12s   %s\n','b','Phi(b)','Phi_loc(b)','argmin');
for b = 0:out.cfg.bMax
    fprintf(fid,'  %3d %12.4f %12.6f   {%s}\n', ...
        b,out.Phi(b+1),out.PhiLoc(b+1),num2str(out.bestB{b+1}));
end
fprintf(fid,'\nWith per-attribute costs\n  %3s',' b');
for c = 1:numel(out.costs), fprintf(fid,' %14s',sprintf('cost %.2f',out.costs(c))); end
fprintf(fid,'\n');
for b = 0:out.cfg.bMax
    fprintf(fid,'  %3d',b);
    for c = 1:numel(out.costs), fprintf(fid,' %14.4f',out.Utot(b+1,c)); end
    fprintf(fid,'\n');
end
for c = 1:numel(out.costs)
    [~,k] = min(out.Utot(2:end,c));
    fprintf(fid,'  cost %.2f -> optimal b = %d, backbone {%s}\n', ...
        out.costs(c),k,num2str(out.bestBcost{k+1,c}));
end
fclose(fid);

fprintf('\nWrote E6.mat, E6_fig.pdf, E6_fig.png, E6_table.txt to:\n  %s\n', outDir);
end

% =======================================================================
function saveFigure(f, outDir, stem)
%SAVEFIGURE  Write a figure as vector PDF and raster PNG, tolerating a
%   destination that is locked by another process.
for ext = {'pdf','png'}
    e = ext{1};
    target = fullfile(outDir,[stem '.' e]);
    try
        writeOne(f,target,e);
    catch ME
        alt = fullfile(outDir, sprintf('%s_%s.%s', stem, ...
                       datestr(now,'yyyymmdd_HHMMSS'), e));                %#ok<TNOW1,DATST>
        fprintf(2,['[E6] could not write %s (%s).\n' ...
                   '     It is probably open in a viewer or being synchronised.\n' ...
                   '     Writing %s instead.\n'], target, ME.identifier, alt);
        try
            writeOne(f,alt,e);
        catch
            fprintf(2,'[E6] the fallback also failed; the figure was not saved.\n');
        end
    end
end
end

function writeOne(f,target,e)
if strcmp(e,'pdf')
    exportgraphics(f, target, 'ContentType','vector');
else
    exportgraphics(f, target, 'Resolution',200);
end
end
