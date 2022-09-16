from pyoptics import madlang

lhc=madlang.open('lhc_as-built_db.seq')
madlang.open('opt_inj.madx',lhc)

deps=lhc.build_dep()
deps2={ item:set(x[0] for x in dep) for item,dep in deps.items()}

import networkx as nx


DG=nx.DiGraph()

for node,edges in deps.items():
    DG.add_node(node)
    for edge in edges:
      expr=edge[-1]
      idx=edge[-2]
      objs='.'.join(edge[:-2])
      if idx is not None:
          objs='%s[%d]'%(objs,idx)
      DG.add_edge(node,objs,idx=idx,expr=expr)


def draw(DG):
  from networkx.drawing.nx_agraph import graphviz_layout
  pos= graphviz_layout(DG,prog='dot')
  nx.draw(DG,pos)
  nx.draw_networkx_labels(DG,pos)

draw(nx.bfs_tree(DG,'on_x1'))


DG=nx.DiGraph()
DG.add_nodes_from(range(6))
DG.add_edges_from([(0,1),(0,2),(2,3),(0,3),(3,4),(1,5)])

dep={1:set([0]),2:set([0]),3:set([2,0]),4:set([3]),5:set([1])}
dep2={0:set([1,2,3]),1:set([5]),3:set([4])}



`
