from pymol import cmd

last_num = 14

cmd.delete('ss')
execfile('py_ss_14.py')
cmd.clip('slab', 2000)
mv = cmd.get_view()

for i in range(last_num+1):
    cmd.delete('ss')
    execfile('py_ss_%d.py' % (i))
    cmd.set_view(mv)
    cmd.ray()
    cmd.png('img_ss_%d.png' % (i))

'''
cmd.delete('ss')
run py_ss_0.py
cmd.set_view(mv)
ray
png('step_0.py')
'''
