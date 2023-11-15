%Returns a structure with tables, matrices and neuron names for the
%forward & backwards connectome of the C. elegans
function [Celegans]=Organizer(expression_type,diagonal_symmetry,off_diagonal_symmetry)

%Creattion of name ordering for neuron based on Natt. Comm paper
Celegans.NeuronNames.NodesB_gap = {'VA02';'VA03';'VA06';'VA07';'VA08';'VA09';'VA10';'VA11';'VA12';'DA03';'DA06';'DA07';'DA05';'DA09';'DA08';'DA04';'DA01';'VA01';'DA02';'VA04';'VA05';'RIML';'RIMR';'AIBL';'AIBR';'AVEL';'AVER';'AVAL';'AVAR'};
Celegans.NeuronNames.NodesB_chem = {'VA01';'VA02';'VA03';'VA04';'VA05';'DA01';'DA02';'DA03';'DA04';'DA05';'DA08';'DA09';'VA06';'VA11';'VA08';'VA09';'VA12';'VA10';'DA06';'DA07';'AVEL';'AVER';'AVDL';'AVDR';'RIML';'RIMR';'AIBL';'AIBR';'AVAL';'AVAR'};
if expression_type=="Chem"
   Celegans.NeuronNames.NodesB = unique(vertcat( Celegans.NeuronNames.NodesB_chem, Celegans.NeuronNames.NodesB_gap),'stable'); %chem based
elseif expression_type=="Gap"
   Celegans.NeuronNames.NodesB = unique(vertcat( Celegans.NeuronNames.NodesB_gap, Celegans.NeuronNames.NodesB_chem),'stable'); %gap based
end

Celegans.NeuronNames.NodesF_gap = {'DB05';'DB06';'DB07';'VB10';'VB11';'DB04';'VB02';'DB01';'VB04';'DB03';'VB05';'DB02';'VB01';'VB06';'VB07';'VB03';'VB08';'VB09';'RIBL';'RIBR';'AVBL';'AVBR'};
Celegans.NeuronNames.NodesF_chem = {'VB06';'VB08';'VB09';'VB07';'VB03';'VB04';'VB05';'VB10';'VB11';'DB02';'DB03';'DB04';'DB06';'DB07';'PVCL';'PVCR';'AVBL';'AVBR';'RIBL';'RIBR'};
if expression_type=="Chem"
   Celegans.NeuronNames.NodesF = unique(vertcat( Celegans.NeuronNames.NodesF_chem, Celegans.NeuronNames.NodesF_gap),'stable'); %chem based
elseif expression_type=="Gap"
   Celegans.NeuronNames.NodesF = unique(vertcat( Celegans.NeuronNames.NodesF_gap, Celegans.NeuronNames.NodesF_chem),'stable'); %chem gap
end

%%%%%%%%%%%%%%% diagonal terms

if diagonal_symmetry=="Hernan Symmetric Weighted"
    files = ["FCVWLRHS.tsv","BCVWLRHS.tsv","FGVWLRHS.tsv","BGVWLRHS.tsv"];
elseif diagonal_symmetry=="Varshney Weighted"
    files = ["Forward Chemical No Symmetry Varshney Weights.tsv","Backwards Chemical No Symmetry Varshney Weights.tsv","Forward Gap No Symmetry Varshney Weights.tsv","Backwards Gap No Symmetry Varshney Weights.tsv"];
elseif diagonal_symmetry=="LR Symmetric Weighted"
    files = ["Forward Chemical LR Symmetry Varshney Weights.tsv","Backwards Chemical LR Symmetry Varshney Weights.tsv","Forward Gap LR Symmetry Varshney Weights.tsv","Backwards Gap LR Symmetry Varshney Weights.tsv"];
end

diag_tables = ["FC","BC","FG","BG"];
diag_matrices = ["W_ff_chem","W_bb_chem","W_ff_gap","W_bb_gap"];
diag_names = ["NodesF","NodesB","NodesF","NodesB"];

for i = 1:4
    T = readtable("files/"+files(i),'FileType',"text");
    G = digraph(table2array(T(:,1)),table2array(T(:,2)),table2array(T(:,3)));
    A = full(adjacency(G,"weighted"));
    Celegans.Tables.(diag_tables(i)) = array2table(A,"RowNames",table2array(G.Nodes),"VariableNames",table2array(G.Nodes));
end

for i=1:4
    Celegans.Matrices.(diag_matrices(i)) = zeros(size(Celegans.NeuronNames.(diag_names(i)),1));

    NodesJ = Celegans.NeuronNames.(diag_names(i));
    NodesK = Celegans.NeuronNames.(diag_names(i));
    Table = Celegans.Tables.(diag_tables(i));
    for j = 1:size(NodesJ,1)
        Nodej = string(NodesJ(j,1));
        for k = 1:size(NodesK,1)
            Nodek = string(NodesK(k,1));
            try
                a = table2array(Table(Nodej,Nodek));
            catch
                a = 0;
            end
            Celegans.Matrices.(diag_matrices(i))(j,k) = a;
        end
    end
end
%%%%%%%%%%%%%%%% off diagonal terms
if off_diagonal_symmetry=="LR Symmetric Weighted"
    files = ["F2B GAP Hernans symmetry neurons with Varshney weights LR symmetrized.tsv","B2F GAP Hernans symmetry neurons with Varshney weights LR symmetrized.tsv","F2B CHEM Hernans symmetry neurons with Varshney weights LR symmetrized.tsv","B2F CHEM Hernans symmetry neurons with Varshney weights LR symmetrized.tsv"];
elseif off_diagonal_symmetry=="Varshney Weighted"
    files = ["F2B GAP Hernans symmetry neurons with Varshney weights.tsv","B2F GAP Hernans symmetry neurons with Varshney weights.tsv","F2B CHEM Hernans symmetry neurons with Varshney weights.tsv","B2F CHEM Hernans symmetry neurons with Varshney weights.tsv"];
end

off_diag_tables = ["F2BG","B2FG","F2BC","B2FC"];
off_diag_matrices = ["W_fb_gap","W_bf_gap","W_fb_chem","W_bf_chem"];
off_diag_names = ["NodesF","NodesB","NodesF","NodesB"];

for i = 1:4
    T = readtable("files/"+files(i),'FileType',"text");
    G = digraph(table2array(T(:,1)),table2array(T(:,2)),table2array(T(:,3)));
    A = full(adjacency(G,"weighted"));
    Celegans.Tables.(off_diag_tables(i)) = array2table(A,"RowNames",table2array(G.Nodes),"VariableNames",table2array(G.Nodes));
end

for i=1:4
    Celegans.Matrices.(off_diag_matrices(i)) =  zeros(size(Celegans.NeuronNames.(off_diag_names(i)),1),size(Celegans.NeuronNames.(off_diag_names(5-i)),1));
    NodesJ = Celegans.NeuronNames.(off_diag_names(i));
    NodesK = Celegans.NeuronNames.(off_diag_names(5-i));%cycle backwards
    Table = Celegans.Tables.(off_diag_tables(i));
    for j = 1:size(NodesJ,1)
        Nodej = string(NodesJ(j,1));
        for k = 1:size(NodesK,1)
            Nodek = string(NodesK(k,1));
            try
                a = table2array(Table(Nodej,Nodek));
            catch
                a = 0;
            end
            Celegans.Matrices.(off_diag_matrices(i))(j,k) = a;
        end
    end
end

end 