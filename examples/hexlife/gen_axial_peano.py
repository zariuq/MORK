def generate_peano_axial():
    # Patterns for Inc/Dec/Same
    inc = [("(P $x)", "$x"), ("Z", "(S Z)"), ("(S $x)", "(S (S $x))")]
    dec = [("(S $x)", "$x"), ("Z", "(P Z)"), ("(P $x)", "(P (P $x))")]
    same = [("$x", "$x")]

    # Axial Directions (dq, dr)
    directions = [
        (1, 0),   # inc q, same r
        (1, -1),  # inc q, dec r
        (0, -1),  # same q, dec r
        (-1, 0),  # dec q, same r
        (-1, 1),  # dec q, inc r
        (0, 1)    # same q, inc r
    ]

    rules = []
    
    for dq, dr in directions:
        pats_q = inc if dq == 1 else (dec if dq == -1 else same)
        pats_r = inc if dr == 1 else (dec if dr == -1 else same)
        
        for in_q, out_q in pats_q:
            iq = in_q.replace("$x", "$q")
            oq = out_q.replace("$x", "$q")
            
            for in_r, out_r in pats_r:
                ir = in_r.replace("$x", "$r")
                or_ = out_r.replace("$x", "$r")
                
                # FIXED: Flatten arguments to match (neighbor q r nq nr)
                rules.append(f"(neighbor {iq} {ir} {oq} {or_})")

    return rules

if __name__ == "__main__":
    for r in generate_peano_axial():
        print(r)