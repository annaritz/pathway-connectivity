```
python3 count.py 
```

produces this table:

```
\begin{table}[h]
\centering
\begin{tabular}{|r||c|c|c|cc|} \hline
& Directed  & Bipartite  & Compound  & \multicolumn{2}{c|}{Hypergraph} \\ 
&  Graph &  Graph &  Graph & Full & Filtered \\ \hline
\# Nodes & 12086 &  30775 &  19650 & 19650 & 15440 \\
\# Edges/Hyperedges & 285556 &  45155 &  38218 & 11125 & 8773 \\ \hline
\end{tabular}
\caption{\textbf{Representations of the Reactome pathway database.}  The filtered hypergraph has removed all small molecules, two forms of Ubiquitinase, and the Nuclear Pore Complex from the hyperedges.}
\label{tab:representations}
\end{table}
```

```
python3 pathways.py
```
