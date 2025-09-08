from functools import partial
from multiprocessing import Pool, cpu_count
from client import MORK, ManagedMORK
import re
from collections import defaultdict
from itertools import combinations

DATASETS = (
    "lordOfTheRings",
    "adameve",
    "simpsons",
    "royal92",
)

def preprocessing(server, datasets=DATASETS):
    with server.work_at("aunt-kg") as ins:
        for dataset in datasets:
            with ins.work_at(dataset) as scope:
                with scope.work_at("src") as src:
                    src.sexpr_import_(f"https://raw.githubusercontent.com/Adam-Vandervorst/metta-examples/refs/heads/main/aunt-kg/{dataset}.metta")\
                        .block()
                    # test: see the names of the people in the dataset
                    downloaded = src.download("(Individuals $i (Fullname $name))", "$name", max_results=5)
                    print("5 names", dataset, downloaded.data)

                scope.transform(("(src (Individuals $i (Id $id)))", "(src (Individuals $i (Fullname $name)))"),
                                ("(simple (hasName $id $name))", "(simple (hasId $name $id))")).block()

                scope.transform(("(src (Individuals $i (Id $id)))", "(src (Individuals $i (Sex \"M\")))"),
                                ("(simple (male $id))",)).block()
                scope.transform(("(src (Individuals $i (Id $id)))", "(src (Individuals $i (Sex \"F\")))"),
                                ("(simple (female $id))",)).block()

                scope.transform(("(src (Relations $r (Husband $id)))", "(src (Relations $r (Children $lci $cid)))"),
                                ("(simple (parent $id $cid))",)).block()
                scope.transform(("(src (Relations $r (Wife $id)))", "(src (Relations $r (Children $lci $cid)))"),
                                ("(simple (parent $id $cid))",)).block()
                # Agh, overloading of keys in JSON, see https://github.com/trueagi-io/MORK/pull/11
                scope.transform(("(src (Relations $r (Husband $id)))", "(src (Relations $r (Children $cid)))"),
                                ("(simple (parent $id $cid))",)).block()
                scope.transform(("(src (Relations $r (Wife $id)))", "(src (Relations $r (Children $cid)))"),
                                ("(simple (parent $id $cid))",)).block()


        path = "file://" + __file__.rpartition("/")[0] + "/simple_all.metta"
        ins.sexpr_export("($dataset (simple $x))", "($dataset $x)", path).block()

    # test: see the steps performed by the server
    # for i, item in enumerate(ins.history):
    #     print("preprocessing event", i, str(item))


def is_different_hack(scope, ids, leave_out_id):
    # note that `leave_out_id` is a prefix and thread-specific, hence allowed to happen in parallel
    with scope.work_at(namespace=f'(simple (isIdDifferent "{leave_out_id}" "{{}}"))') as private_ns:
        private_ns.upload_("\n".join(id0 for id0 in ids if id0 != leave_out_id))

# def hack(server, datasets=DATASETS, parallel=True):
#     # HACK: Generate isIdDifferent relations so that sister relations aren't reflective
#     if parallel: pool = Pool(cpu_count()//2)
#     with server.work_at("aunt-kg") as ins:
#         for dataset in datasets:
#             with ins.work_at(dataset) as scope:
#                 ids_string = scope.download("(src (Individuals $i (Id $id)))", "$id").data
#                 ids = [line.strip()[1:-1] for line in ids_string.strip().splitlines() if line]
#                 if parallel: pool.map(partial(is_different_hack, scope._bare(), ids), ids)
#                 else: scope.upload_("".join(f'(simple (isIdDifferent "{leave_out_id}" "{id0}"))\n'
#                                             for id0 in ids for leave_out_id in ids if id0 != leave_out_id))
#     if parallel: pool.terminate()

def hack(server, datasets=DATASETS, parallel=False):
    with server.work_at("aunt-kg").and_time() as ins:
        for dataset in datasets:
            with ins.work_at(dataset) as scope:
                ids_text = scope.download("(src (Individuals $i (Id $id)))", "$id").data
                ids = [ln.strip()[1:-1] for ln in ids_text.strip().splitlines() if ln]
                n = len(ids)
                print(f"[hack] {dataset}: generating isIdDifferent for {n} IDs (~{n*(n-1)} pairs)")
                if parallel:
                    from multiprocessing import Pool, cpu_count
                    from functools import partial
                    pool = Pool(max(1, cpu_count() // 2))
                    try:
                        pool.map(partial(is_different_hack, scope._bare(), ids), ids)
                    finally:
                        pool.terminate()
                else:
                    K = max(1, n // 50)  # progress reporting chunk
                    for idx, left in enumerate(ids, 1):
                        scope.upload_("".join(
                            f'(simple (isIdDifferent "{left}" "{right}"))\n'
                            for right in ids if right != left
                        ))
                        if idx % K == 0 or idx == n:
                            print(f"[hack] {dataset}: {idx}/{n} ids processed")


def _chunked_upload(scope, lines_iter, chunk_size=200_000):
    """Upload in big chunks so we don't spam the server with tiny requests."""
    buf = []
    for ln in lines_iter:
        buf.append(ln)
        if len(buf) >= chunk_size:
            scope.upload_("".join(buf))
            buf.clear()
    if buf:
        scope.upload_("".join(buf))

def _parse_quoted_pairs(lines_text):
    """Parse lines like: ("ParentID" "ChildID") -> [("ParentID","ChildID"), ...]"""
    out = []
    for ln in lines_text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        parts = re.findall(r'"([^"]*)"', ln)
        if len(parts) == 2:
            out.append((parts[0], parts[1]))
    return out

def siblings_only_hack(server, datasets=("simpsons", "adameve", "lordOfTheRings", "royal92")):
    """Generate (simple (isIdDifferent A B)) only for children with the same parents."""
    with server.work_at("aunt-kg").and_time() as ins:
        for dataset in datasets:
            with ins.work_at(dataset) as scope:
                # download (simple (parent $p $c)) as pairs "(P C)"
                txt = scope.download('(simple (parent $p $c))', '($p $c)').data
                pairs = _parse_quoted_pairs(txt)

                # Build parent->children, and child->set(parents)
                parent_children = defaultdict(set)
                child_parents = defaultdict(set)
                for p, c in pairs:
                    parent_children[p].add(c)
                    child_parents[c].add(p)

                # Group children by exact parent-set (so full-siblings only)
                sibling_groups = defaultdict(set)
                for c, ps in child_parents.items():
                    key = tuple(sorted(ps))  # normalize e.g. ("Mother","Father")
                    sibling_groups[key].add(c)

                # Produce isIdDifferent only within each sibling group
                total_pairs = sum(len(g)*(len(g)-1) for g in sibling_groups.values())
                print(f"[hack-siblings-only] {dataset}: {len(sibling_groups)} sibling groups, ~{total_pairs} ordered pairs")

                def gen_lines():
                    for group in sibling_groups.values():
                        # ordered pairs within the group
                        for a in group:
                            for b in group:
                                if a != b:
                                    yield f'(simple (isIdDifferent "{a}" "{b}"))\n'

                _chunked_upload(scope, gen_lines(), chunk_size=200_000)

def processing(server, datasets=DATASETS, human_readable=True):
    with server.work_at("aunt-kg") as ins:
        for dataset in datasets:
            with ins.work_at(dataset).work_at("simple").and_time() as scope:
                scope.transform(("(parent $pid $cid)", "(female $pid)",),
                                ("(mother $pid $cid)",),).block()
                if human_readable: scope.transform(("(mother $pid $cid)", "(hasName $pid $name0)", "(hasName $cid $name1)",), ("(motherByName $name0 $name1)",),).block()

                scope.transform(("(parent $pid $cid0)", "(parent $pid $cid1)", "(isIdDifferent $cid0 $cid1)", "(female $cid0)",),
                                ("(sister $cid0 $cid1)",),).block()
                if human_readable: scope.transform(("(sister $cid0 $cid1)", "(hasName $cid0 $name0)", "(hasName $cid1 $name1)",), ("(sisterByName $name0 $name1)",),).block()

                scope.transform(("(parent $pid $cid)", "(sister $aid $pid)",),
                                ("(aunt $aid $cid)",),).block()
                if human_readable: scope.transform(("(aunt $aid $cid)", "(hasName $aid $name0)", "(hasName $cid $name1)",), ("(auntByName $name0 $name1)",),).block()

def _main():
    with ManagedMORK.connect(binary_path="../target/release/mork_server").and_terminate() as server:
        server.clear().block()
        preprocessing(server)
        # hack(server, parallel=True)
        siblings_only_hack(server, datasets=("simpsons","adameve","lordOfTheRings","royal92"))
        processing(server, human_readable=False)

if __name__ == '__main__':
    _main()
