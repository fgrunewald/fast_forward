import networkx as nx

def find_one_ismags_match(graph1, graph2, node_match):
    """
    Returns one ismags match when graphs are isomorphic
    otherwise None.
    """
    GM = nx.isomorphism.GraphMatcher(graph1, graph2, node_match=node_match)
    raw_matches = GM.subgraph_isomorphisms_iter()
    try:
        mapping = next(raw_matches)
        return mapping
    except StopIteration:
        raise IOError("no match_found")
        return None

def node_match(n1, n2):
    for attr in ["resname", "atomname"]:
        if n1[attr] != n2[attr]:
            return False
    return True
