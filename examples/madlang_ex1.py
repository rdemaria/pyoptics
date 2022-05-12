from pyoptics import madlang

lhc=madlang.open('hllhc_thin.seq')

print(lhc.lhcb1.mcbyyv_4l1_b1)

lhc.print_dep('on_x1')

print(lhc.acbxv1_r1)

lhc.on_x1=1
print(lhc.acbxv1_r1)
