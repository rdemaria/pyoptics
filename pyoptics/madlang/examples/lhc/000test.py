from pyoptics import madlang

lhc=madlang.open('lhc_as-built_db.seq')
madlang.open('opt_inj.madx',lhc)
lhc.print_dep('on_sep1')
out,rest=lhc.lhcb1.expand_struct()




