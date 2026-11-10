# Window graphs of the 88 classes

`window_graphs.py` (needs NumPy) builds the graphs W^{(k,p)}_M of the paper
for (k,p) = (0,1), (1,1), (-1,1): the vertices are the words of length M+1,
and there is an edge from c_[0,M] to (sigma^k F^p c)_[0,M] for every value of
the cells outside the window that the image depends on. It computes the
strongly connected components.

    python3 window_graphs.py                 # M = 0,...,6, all 88 classes
    python3 window_graphs.py 0,1,2,3 88 95   # M = 0,...,3, representatives between 88 and 95

For the smallest rule of each class it prints one line with three strings,
one for each (k,p); the letter at position M is `S` if W^{(k,p)}_M is
strongly connected and a dot if not. The output for M ≤ 6 is
`window_graphs_M0-6.txt`. It gives the remark on the window graphs of the
classes:

- every class of status (G) other than those of 90, 95 and 165 has a dot in
  the string for (0,1) at some M ≤ 4;
- 90 and 165 have `SSSSSSS` for (0,1) and dots for (1,1) and (-1,1);
- the classes of status (Z) other than that of 46, the class of 95, and all
  open classes other than that of 191 have `SSSSSSS` three times;
- 46 has `S......` and 191 has `SS.....` for (0,1).

Running time: about 20 minutes for M ≤ 6 (M = 6 takes most of it); a few
seconds per class for M ≤ 4.
