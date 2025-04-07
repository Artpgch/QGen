clc
clear

nQ = 3; % number of queuing systems
qN = 12; % number of nodes
spL = 1; % level of sparse
dZ = true; % is the diagonal zero (are there loops)
sM = true; % is the matrix symmetric

syBw = false; % is the bandwidth matrix symmetric
smBw = false; % are the bandwidth same
minBw = 1; % (minimal) bandwidth
maxBw = 9; % (maximal) bandwidth

minFl = 3; % minimal flow
maxFl = 8; % maximal flow
spFl = 3; % sparse level of bandwidth matrix

arrowcolor = [110/255, 159/255, 1];
backcolor = [24/255, 27/255, 31/255];
bwlabelcolor = [155/255, 1, 25/255];
flowlabelcolor = [1, 120/255, 80/255];

for q = 1:nQ
    disp(strcat('Creating a queuing system', 32, num2str(q), 32, 'from', 32, num2str(nQ)))
[nM, fM] = RandQ(qN, spL, dZ, sM, syBw, smBw, minBw, maxBw, minFl, maxFl, spFl);

writematrix(nM, 'Networks.xlsx', 'Sheet', strcat('Network', 32, num2str(q)));
writematrix(fM, 'Flows.xlsx', 'Sheet', strcat('Flows', 32, num2str(q)));

sG = simpGraph(nM);
nG = netGraph(nM);
[nfl, nl, fG] = flowGraph(fM);

figure('Name', strcat('Simplified network', 32, num2str(q)),'NumberTitle','off','MenuBar','none', 'Visible','off')
plot(sG,'NodeColor', [0 1 0], 'EdgeFontSize', 8, 'EdgeColor', arrowcolor, 'LineWidth', 0.8,'EdgeAlpha', 1, 'NodeFontSize',12,'NodeLabelColor',[1 1 1], 'MarkerSize', 2, 'NodeFontSize', 8, 'EdgeFontSize', 8)
title(strcat('Simplified network', 32, num2str(q)))
set(gca, 'color', backcolor)
exportgraphics(gcf, 'SimpNetwork.pdf','Append', true)

figure('Name', strcat('Network', 32, num2str(q)),'NumberTitle','off','MenuBar','none', 'Visible','off')
plot(nG, 'Layout', 'force', 'EdgeLabel', nG.Edges.Weight, 'NodeColor', [0 1 0], 'EdgeFontSize', 6, 'EdgeColor', arrowcolor, 'ArrowSize', 8,'NodeFontSize', 8,'NodeLabelColor',[1 1 1],'LineWidth', 0.8, 'MarkerSize', 2,'EdgeLabelColor', bwlabelcolor,'EdgeAlpha',1);
title(strcat('Network', 32, num2str(q)))
set(gca, 'color', backcolor)
exportgraphics(gcf, 'Networks.pdf','Append', true)

figure('Name', strcat('Flows', 32, num2str(q)),'NumberTitle','off','MenuBar','none', 'Visible','off')
plot(fG, 'NodeLabel', nl, 'EdgeLabel', fG.Edges.Weight,'NodeColor', [0 1 0],'ArrowSize', 8, 'EdgeFontSize', 8, 'EdgeColor', arrowcolor, 'LineWidth', 0.8,'EdgeAlpha', 1, 'NodeFontSize', 8,'NodeLabelColor',[1 1 1], 'EdgeFontSize', 8, 'EdgeLabelColor', flowlabelcolor,'EdgeAlpha', 1, 'EdgeFontSize', 11, 'NodeFontSize', 8, 'EdgeFontSize', 8, 'MarkerSize', 2);
title(strcat('Flows', 32, num2str(q)))
set(gca, 'color', backcolor)
exportgraphics(gcf, 'Flows.pdf','Append', true)
clc
end

if smBw
    maxBw = minBw;
end

disp(strcat(num2str(nQ), 32,'queuing systems',  32, 'created', 32, 'with parameters:', 10,  10, ...
    'Number of nodes:', 32, num2str(qN), 10, ...
    'Level of sparsity of the adjacency matrices:', 32, num2str(spL), 10, ...
    'Minimal bandwidth:', 32, num2str(minBw), 10, ...
    'Maximal bandwidth:', 32, num2str(maxBw), 10, ...
    'Minimal flow:', 32, num2str(minFl), 10, ...
    'Maximal flow:', 32, num2str(maxFl), 10, ...
    'Level of sparsity of the flow matrices:', 32, num2str(spFl) ))
if dZ
    disp('Non-zero diagonal')
else
    disp('Zero diagonal')
end

if sM
    disp('Symmetric adjacency matrix')
else
    disp('Asymmetric adjacency matrix')
end

if syBw
    disp('Symmetric bandwidth matrix')
else
    disp('Asymmetric bandwidth matrix')
end


function sG = simpGraph(nM)
qN = length(nM);
rM = nM ./ nM;

for a = 1 : qN
    for b = 1 : qN
        if isnan(rM(a, b))
            rM(a, b) = 0;
        end
    end
end

ruM = triu(rM);
trM = rM- ruM;
rrM = zeros(qN);
for p= 1:qN
    for q = 1:qN
        rrM(q, p) = trM(p,q);
    end
end
grM = rrM + ruM;

t = zeros(1, nnz(grM));
s = t;
bW = s;
u = 1;
for p = 1 : qN
    for q = 1 : qN
        if grM(p, q) > 0
                s(u) = p;
                t(u) = q;        
                bW(u) = grM(p, q);  
                u = u + 1;
        end
     end
end
sG = graph(s,t);

end

function nG = netGraph(nM)
qN = length(nM);
t = zeros(1, nnz(nM));
s = t;
bW = s;
u = 1;
for p = 1 : qN
    for q = 1 : qN
        if nM(p, q) > 0
                s(u) = p;
                t(u) = q;       
                bW(u) = nM(p, q);   
                u = u + 1;
        end
     end
end

nG = digraph(s, t, bW);

end

function [nfl, nl, fG] = flowGraph(nM)

qN = length(nM);
nfl = nnz(nM);
nl = zeros(1,2*nfl);
flG(1,:) = 1:2:2*nfl-1;
flG(2,:) = 2:2:2*nfl;
flI = zeros(1,nfl);
c = 1;
for a = 1 : qN
    for b = 1 : qN
        if nM(a, b) > 0
            nl(2*c - 1) = a;
            nl(2*c) = b;
            flI(c) = nM(a, b);
            c = c + 1;
        end
    end
end
fG = digraph(flG(1,:), flG(2,:), flI);

end
