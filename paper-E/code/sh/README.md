# The shortening statement (Sh) for rules 91 and 94

Computation for the remark on responses in the subsection "Responses to
probes". It is not used in a proof. A point of Φ^+ (probe) moves to the
right, a point of Φ^- (object) to the left; all scatterings are computed by
direct simulation of F in the coordinates (t, i) (`sc.h`). All running times
are for one core.

    gcc -O2 -o sh sh.c

| statement | program | command | expected output |
|---|---|---|---|
| (Sh), rule 94, cores of length ≤ 23, probes of length ≤ 24, the four tail types | `sh.c` | `./sh 94 c0 c1 d0 d1 1 23 24 0` for (c0 c1 d0 d1) = 0000, 0011, 1100, 1111 | `logs/sh_all.log`: `fail 0` in every line |
| (Sh), rule 91, cores of length ≤ 16, probes of length ≤ 20, the nine tail types | `sh.c` | `./sh 91 c0 c1 d0 d1 1 16 20 0`, c, d ∈ {01, 10, 11} | `logs/sh_all.log`: `fail 0` in every line |

`sh RULE c0 c1 d0 d1 NMIN NMAX LP 0`: for every object with left tail c,
right tail d and core of length NMIN..NMAX, searches the probes of length
≤ LP for one that shortens the object. When the two tails are different,
only probes that send the two bare kinks with these tails to bare kinks are
used, so the program tests a statement that implies (Sh). It prints one line
`n=… objects … fail … maxneed …` per core length, `maxneed` being the
largest length of the first probe found. The first line gives the number of
probes (452060 and 793310 for rule 94 on the backgrounds 00 and 11 with
LP = 24; 167914 and 203804 for rule 91 on the backgrounds 01 and 11 with
LP = 20). `./sh 94 0 0 0 0 1 10 12 0` takes about 1 second and prints the
lines n = 1..10 of the log (objects 0, 1, 1, 1, 2, 4, 7, 12, 21, 37); the
full runs take minutes to hours.
