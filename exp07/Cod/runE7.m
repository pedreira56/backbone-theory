function out = runE7(quick)
%RUNE7  Driver for experiment E7 (Theorem T5, the completed dataset).
%
%   OUT = RUNE7()      full run   (use this for the paper)
%   OUT = RUNE7(true)  quick run  (smoke test, about one minute)
%
%   Layout: put every .m file in  ...\Exp 7\Cod\  and results appear in
%   ...\Exp 7\Result\ , created automatically.  See BBOUTDIR.
%   Copy there the five bb*.m files as well.
%
%   Outputs: E7.mat, E7_fig.pdf/.png, E7_table.txt
%
%   Panel (a) fixes the number of receivers and grows the donor block: the
%   error falls and then flattens on the Monte-Carlo floor of that N.
%   Panel (b) fixes the donor block and grows the number of receivers: the
%   error falls at the parametric rate and then flattens on a plateau set by
%   the matching error of that n.  Together they show that neither index
%   alone suffices.  The third curve of panel (a) grows both, along
%   N = 100 n, and the error vanishes: the theorem requires no relation
%   between them, and donors reused a hundred times do not impede it.
%   Panel (c) is the quantity that makes the completed rows dependent.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E7] code   : %s\n[E7] results: %s\n\n', codeDir, outDir);

cfg = struct();
if quick
    cfg.nGrid = round(logspace(2,4,5)); cfg.NGrid = round(logspace(2,4,5));
    cfg.repBase = 120; cfg.repMin = 6; cfg.repMax = 24; cfg.nOracle = 2e5;
else
    cfg.nGrid = round(logspace(2,5,7)); cfg.NGrid = round(logspace(2,5,7));
    cfg.repBase = 400; cfg.repMin = 8; cfg.repMax = 60; cfg.nOracle = 2e6;
end

tStart = tic;
out = E7_completedBase(cfg);
out.outDir  = outDir;
out.seconds = toc(tStart);
fprintf('\n[E7] elapsed: %.1f min\n', out.seconds/60);

nG = out.nGrid; NG = out.NGrid; kF = 2;     % f2, the cosine, is the one plotted
rmse = sqrt(out.mse(:,:,kF));

% indices used by the panels
[~,bLo] = min(abs(NG - 1e3));
[~,bHi] = min(abs(NG - NG(end)));
[~,aLo] = min(abs(nG - 1e2));
[~,aMd] = min(abs(nG - 1e3));
[~,aHi] = min(abs(nG - 1e4));

% the diagonal N = 100 n, taken from the cells of the grid that realise it
dn = []; dv = [];
for a = 1:numel(nG)
    [dd,b] = min(abs(NG - 100*nG(a)));
    if dd/max(NG(b),1) < 0.1
        dn(end+1) = nG(a);  dv(end+1) = rmse(a,b);  %#ok<AGROW>
    end
end

f = figure('Position',[80 80 1350 420]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

% ---- (a) growing the donor block --------------------------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
plot(ax,nG,rmse(:,bLo),'o-','Color',cols(1,:),'MarkerSize',4.5,'LineWidth',1.3, ...
     'MarkerFaceColor','w','DisplayName',sprintf('N = %d receivers',NG(bLo)));
plot(ax,nG,rmse(:,bHi),'o-','Color',cols(2,:),'MarkerSize',4.5,'LineWidth',1.3, ...
     'MarkerFaceColor','w','DisplayName',sprintf('N = %d receivers',NG(bHi)));
if ~isempty(dn)
    plot(ax,dn,dv,'s-','Color',cols(3,:),'MarkerSize',5.5,'LineWidth',1.4, ...
         'MarkerFaceColor','w','DisplayName','N = 100 n');
end
set(ax,'XScale','log','YScale','log','FontSize',9);
xlabel(ax,'n   (donor block size)');
ylabel(ax,'root mean squared error');
text(ax,-0.02,1.05,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southwest','FontSize',8,'Box','off');

% ---- (b) growing the number of receivers -------------------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
for i = [aLo aMd aHi]
    plot(ax,NG,rmse(i,:),'o-','MarkerSize',4.5,'LineWidth',1.3, ...
         'MarkerFaceColor','w','DisplayName',sprintf('n = %d donors',nG(i)));
end
ref = rmse(aHi,1)*sqrt(NG(1)./NG);
plot(ax,NG,ref,'--','Color',[.45 .45 .45],'LineWidth',1.1, ...
     'DisplayName','slope -1/2');
set(ax,'XScale','log','YScale','log','FontSize',9);
xlabel(ax,'N   (number of receivers)');
ylabel(ax,'root mean squared error');
text(ax,-0.02,1.05,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southwest','FontSize',8,'Box','off');

% ---- (c) the probability of a shared donor -----------------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
pv = out.pi(:,bHi);
plot(ax,nG,pv,'o-','Color',[0.49 0.18 0.56],'MarkerSize',4.5,'LineWidth',1.3, ...
     'MarkerFaceColor','w','DisplayName','\pi_n, measured exactly');
refp = pv(1)*(log(nG)./nG)/(log(nG(1))/nG(1));
plot(ax,nG,refp,'--','Color',[.45 .45 .45],'LineWidth',1.1, ...
     'DisplayName','log n / n');
set(ax,'XScale','log','YScale','log','FontSize',9);
xlabel(ax,'n   (donor block size)');
ylabel(ax,'probability that two receivers share a donor');
text(ax,-0.02,1.05,'(c)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','southwest','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E7_fig.pdf'),'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E7_fig.png'),'Resolution',200);

% ---- console and table --------------------------------------------------
L = {};
L{end+1} = 'E7 - Theorem T5, convergence of the completed dataset';
L{end+1} = '';
L{end+1} = 'Limits of the three test functions';
for k = 1:numel(out.thetaPerp)
    L{end+1} = sprintf('  %-16s product coupling %+.6f   population %+.6f   gap %+.6f', ...
        out.fNames{k}, out.thetaPerp(k), out.thetaTrue(k), out.thetaTrue(k)-out.thetaPerp(k));
end
L{end+1} = '';
L{end+1} = 'Root mean squared error of the completed-base mean of f1';
hdr = sprintf('  %8s','n \ N');
for b = 1:numel(NG), hdr = [hdr sprintf(' %11d',NG(b))]; end %#ok<AGROW>
L{end+1} = hdr;
for a = 1:numel(nG)
    row = sprintf('  %8d',nG(a));
    for b = 1:numel(NG), row = [row sprintf(' %11.3e',rmse(a,b))]; end %#ok<AGROW>
    L{end+1} = row;
end
L{end+1} = '';
L{end+1} = 'Probability that two receivers share a donor, and the reference log n / n';
L{end+1} = sprintf('  %8s %13s %13s %10s','n','pi_n','log n / n','ratio');
for a = 1:numel(nG)
    r = log(nG(a))/nG(a);
    L{end+1} = sprintf('  %8d %13.3e %13.3e %10.3f',nG(a),out.pi(a,bHi),r,out.pi(a,bHi)/r);
end
L{end+1} = '';
L{end+1} = sprintf('The bound of Theorem T5, term by term (f2, the cosine)');
L{end+1} = sprintf('  Lambda_2^2 = %.4f',out.lam2);
L{end+1} = sprintf('  %8s %9s %11s %11s %11s %11s %9s','n','N','measured','4/N','Delta_n^2','8 pi_n','ratio');
for a = 1:numel(nG)
    for b = [bLo bHi]
        L{end+1} = sprintf('  %8d %9d %11.3e %11.3e %11.3e %11.3e %9.4f', ...
            nG(a),NG(b),out.mse(a,b,2),4/NG(b),out.Delta2(a,b),8*out.pi(a,b), ...
            out.mse(a,b,2)/out.bound(a,b));
    end
end
L{end+1} = '';
L{end+1} = sprintf('Along N = 100 n');
for i = 1:numel(dn)
    L{end+1} = sprintf('  n = %7d  N = %9d   rmse = %.3e',dn(i),100*dn(i),dv(i));
end
L{end+1} = '';
L{end+1} = sprintf('Elapsed: %.1f min',out.seconds/60);

fprintf('\n'); fprintf('%s\n',L{:});
fid = fopen(fullfile(outDir,'E7_table.txt'),'w');
fprintf(fid,'%s\n',L{:});
fclose(fid);

save(fullfile(outDir,'E7.mat'),'-struct','out');
fprintf('\nWrote E7.mat, E7_fig.pdf, E7_fig.png, E7_table.txt to:\n  %s\n', outDir);
end
