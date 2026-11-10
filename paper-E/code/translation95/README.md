# Rule 95: the translation property

`t95ex.py K D` takes all pairs (P, X) of points on the background 1 whose
cores (the words between the first and the last zero, without two zeros at
distance 2) have length ≤ K, at both relative parities of placement, and
checks that in S(P, X) = (Γ, Ω) the zeros of each class of X are translated
by one common amount in Γ, and those of P in Ω. With D = 0 the scattering is
computed with the Goursat recurrence (`gs.c`, called through the module
`sim.py`); with D = 1 it is computed in
addition by direct simulation in (t, i) (`sim.direct`). Needs numpy.

    gcc -O2 -shared -fPIC -o libgs.so gs.c
    python3 t95ex.py 9 0        # K 9 pairs 3200 failures 0
    python3 t95ex.py 6 1        # K 6 pairs 200 failures 0

Running time: under 1 second each.
