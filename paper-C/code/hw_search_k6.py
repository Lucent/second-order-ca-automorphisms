"""The search of hw_search.py for windows of length 6 and closures of at most 300 words, on the 57 classes
without a fixed word of length at most 3. Output: logs/hw_search_k6.txt. Needs about 2 GB of memory.
"""
import sys, json
sys.argv=['x','6']
import hw_search as S
L=[13,14,15,24,25,26,27,28,29,30,31,41,43,45,46,47,57,58,59,60,61,62,63,75,77,79,90,91,94,95,105,107,109,111,123,126,127,133,135,141,143,151,153,155,157,159,165,167,175,179,183,187,189,191,219,223,255]
out={}
for r in L:
    fnd,sizes=S.analyse(r,6,cap=300)
    out[r]=[(int(a),int(b)) for a,b,_ in fnd]
    print(r, sorted(S.cls(r)), 'min closure',int(sizes.min()), 'hidden-involution families (size,j):', sorted(set(out[r])), flush=True)
json.dump(out,open('search_k6.json','w'))
