function figE1(hint)
%FIGE1  Redraw the figure of experiment E1 from the saved results.
%
%   FIGE1() finds E1.mat (see BBFINDMAT), redraws the figure and writes
%   E1_fig_v2.pdf and E1_fig_v2.png next to it.  Nothing is simulated again.
%   FIGE1(FOLDER) searches FOLDER first.
%
%   Needs only this file and bbFindMat.m; it does not call bbParams.
%
%   WHY A SECOND VERSION
%   The original panel (a) carried ten curves and a legend that covered the
%   plot area.  Here the curves are split by geometry into two panels, the
%   legends sit below the axes, and the fitted dimensions are left to the
%   table, where they can be read to three figures.
%     (a) cube and flat torus, paired at d = 3, 5, 8
%     (b) curved embeddings, with the plain square as reference: the square
%         and the embedding with d = 2, b = 10 share a colour and a slope
%     (c) the tail bound, unchanged

if nargin < 1, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E1', hint);
out = load(matFile);
k1  = find(out.qList == 1, 1);

dVals = [1 2 3 4 5 8];
dCols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
         0.49 0.18 0.56; 0.00 0.60 0.60; 0.64 0.08 0.18];

% Classify the regimes from the stored option lists, so that the figure can
% be redrawn without the population code on the path.
nR = numel(out.regimes);
isGraph = false(nR,1); isTorus = false(nR,1);
dd = out.dTrue(:); bb = zeros(nR,1);
for g = 1:nR
    p = out.regimes(g).params;
    isGraph(g) = strcmp(getopt(p,'embedding','identity'),'graph');
    isTorus(g) = strcmp(getopt(p,'geometry','cube'),'torus');
    bb(g)      = getopt(p,'b',dd(g));
end

f = figure('Position',[80 80 1350 430]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');

% ---- (a) cube and torus ------------------------------------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
for g = find(~isGraph).'
    if isTorus(g), ls = '--'; nm = 'torus'; else, ls = '-'; nm = 'cube'; end
    plot(ax, out.nGrid, squeeze(out.moment(g,:,k1)), ls, ...
         'Color',dCols(dVals == dd(g),:), 'Marker','o','MarkerSize',3.5, ...
         'LineWidth',1.1, 'DisplayName',sprintf('%s, d_B = %d',nm,dd(g)));
end
set(ax,'XScale','log','YScale','log','FontSize',9);
xlabel(ax,'n   (donor block size)'); ylabel(ax,'E[ R ]');
text(ax,-0.02,1.05,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southoutside','NumColumns',3,'FontSize',8,'Box','off');

% ---- (b) curved embeddings, with the square as reference ---------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
gRef = find(~isGraph & ~isTorus & dd == 2, 1);
sel  = [gRef(:); find(isGraph)];
for g = sel.'
    if isGraph(g)
        ls = '-.'; nm = sprintf('curved, d_B = %d, b = %d',dd(g),bb(g));
    else
        ls = '-';  nm = sprintf('cube, d_B = b = %d',dd(g));
    end
    plot(ax, out.nGrid, squeeze(out.moment(g,:,k1)), ls, ...
         'Color',dCols(dVals == dd(g),:), 'Marker','o','MarkerSize',3.5, ...
         'LineWidth',1.1, 'DisplayName',nm);
end
set(ax,'XScale','log','YScale','log','FontSize',9);
xlabel(ax,'n   (donor block size)'); ylabel(ax,'E[ R ]');
text(ax,-0.02,1.05,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southoutside','NumColumns',1,'FontSize',8,'Box','off');

% ---- (c) tail bound ----------------------------------------------------
T    = out.tail;
keep = T.surv > 0 & T.surv < 1;
xx   = T.eps(keep).^T.dB;
yy   = -log(T.surv(keep))/T.n;
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
% Both axes are rescaled by the same factor, which leaves the slope -- the
% quantity being checked -- untouched and removes the 10^{-5} exponents that
% MATLAB would otherwise print in the top corners, where one of them collided
% with the panel label.
sc = 1e5;
plot(ax, sc*xx, sc*T.Vtheory*xx, '--','Color',[.45 .45 .45],'LineWidth',1.2, ...
     'DisplayName',sprintf('theoretical  V_d = %.3f',T.Vtheory));
plot(ax, sc*xx, sc*yy, 'o','Color',dCols(dVals == T.dB,:),'MarkerSize',3.5, ...
     'LineWidth',1,'DisplayName',sprintf('empirical  c_B = %.3f',T.cFit));
set(ax,'FontSize',9);
xlabel(ax,sprintf('10^{5} \\epsilon^{d_B}     (torus, d_B = %g, n = %d)',T.dB,T.n));
ylabel(ax,'-10^{5} log P(R > \epsilon) / n');
text(ax,-0.02,1.05,'(c)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southoutside','NumColumns',1,'FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E1_fig_v2.pdf'),'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E1_fig_v2.png'),'Resolution',200);
fprintf('Wrote E1_fig_v2.pdf and E1_fig_v2.png to:\n  %s\n',outDir);
end

% =======================================================================
function v = getopt(c, name, dflt)
%GETOPT  Value of an option in a name/value cell array, or a default.
v = dflt;
for i = 1:2:numel(c)-1
    if ischar(c{i}) && strcmpi(c{i},name), v = c{i+1}; end
end
end
