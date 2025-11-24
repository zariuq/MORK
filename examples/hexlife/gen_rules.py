def generate_rules():
    # Patterns
    inc_patterns = [
        ("(P $x)", "$x"),
        ("Z", "(S Z)"),
        ("(S $x)", "(S (S $x))")
    ]
    dec_patterns = [
        ("(S $x)", "$x"),
        ("Z", "(P Z)"),
        ("(P $x)", "(P (P $x))")
    ]
    nop_pattern = [("$x", "$x")]

    # Directions (index 0=q, 1=r, 2=s)
    # Changes: +1 (inc), -1 (dec), 0 (nop)
    directions = [
        (1, 0, -1), # q+, s-
        (1, -1, 0), # q+, r-
        (0, -1, 1), # r-, s+
        (-1, 0, 1), # q-, s+
        (-1, 1, 0), # q-, r+
        (0, 1, -1)  # r+, s-
    ]

    rules = []
    
    for dq, dr, ds in directions:
        # Get patterns for each coord
        pat_q = inc_patterns if dq == 1 else (dec_patterns if dq == -1 else nop_pattern)
        pat_r = inc_patterns if dr == 1 else (dec_patterns if dr == -1 else nop_pattern)
        pat_s = inc_patterns if ds == 1 else (dec_patterns if ds == -1 else nop_pattern)

        for iq, oq in pat_q:
            # Fix variable names in patterns to match coord
            iq_s, oq_s = (iq.replace("$x", "$q"), oq.replace("$x", "$q")) if dq != 0 else ("$q", "$q")
            
            for ir, or_ in pat_r:
                ir_s, or_s = (ir.replace("$x", "$r"), or_.replace("$x", "$r")) if dr != 0 else ("$r", "$r")
                
                for is_, os_ in pat_s:
                    is_s, os_s = (is_.replace("$x", "$s"), os_.replace("$x", "$s")) if ds != 0 else ("$s", "$s")
                    
                    rule = f"(neighbors ({iq_s} {ir_s} {is_s}) ({oq_s} {or_s} {os_s}))"
                    rules.append(rule)
    
    return rules

for r in generate_rules():
    print(r)
