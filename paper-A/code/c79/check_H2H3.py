# Confirmation (not used in the proof): the proofs of Propositions (H_2) and (H_3), run as recipes
# (recipe_H2.py, recipe_H3.py) at time 0.
#  (1) For every level-L state, in three-valued logic (recipe_states.A): a run that completes certifies the steps
#      for all configurations with that state. L=4 for (H_2); L=6 for (H_3), the states for which a value at
#      level 6 is unknown being completed at level 7.
#  (2) The same recipes on random whole configurations, simulated exactly (recipe_states.C), at the times |t|<=3.
# Usage: python3 check_H2H3.py seed trials
import sys, random
import recipe_H2, recipe_H3
from recipe_states import A, C, Unknown, states, vcell, Conf
seed=int(sys.argv[1]) if len(sys.argv)>1 else 11
trials=int(sys.argv[2]) if len(sys.argv)>2 else 1500
n=0
for s in states(4): recipe_H2.h2(A(s,4),0); n+=1
print('(H_2) recipe, all level-4 states:',n,flush=True)
limited=[]; n=0
for s in states(6):
    try: recipe_H3.prove(A(s,6),0); n+=1
    except Unknown: limited.append(s)
print('(H_3) recipe, level-6 states completed:',n,' limited:',len(limited),flush=True)
n7=0
for s in limited:
    for e in [(0,0),(0,1),(1,0),(1,1)]:
        recipe_H3.prove(A((s[0],(e,)+s[1]),7),0); n7+=1
print('(H_3) recipe, level-7 extensions of the limited states completed:',n7,flush=True)
rng=random.Random(seed); N=0
for trial in range(trials):
    th=rng.randrange(6)
    left=[(rng.randint(0,1),rng.randint(0,1)) for _ in range(110)]
    while left[-1]==vcell(th,0): left[-1]=(rng.randint(0,1),rng.randint(0,1))
    c=C(Conf(th,left,PADR=80))
    for t in range(-3,4): recipe_H2.h2(c,t); recipe_H3.prove(c,t); N+=1
print('cases of (H_3) Steps used:',sorted(recipe_H3.cnt.items(),key=str))
print('(H_2) and (H_3) recipes on random whole configurations: (configuration, time) pairs',N)
