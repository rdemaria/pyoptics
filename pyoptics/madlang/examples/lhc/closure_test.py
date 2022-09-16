
def mkclosure(env):
   def wrapper(fn):
     def newf(*args,**nargs):
       return fn(*args,**nargs)
     return newf
   return wrapper

def mkclosure3(env,expr):
   code=compile(expr,expr,'eval')
   def newf():
       return eval(code,globals(),locals())
   return newf


test3=mkclosure3(env,'env[3]')

@mkclosure(env)
def test(val):
    env[3]=val

env={}
test2=mkclosure(env)(lambda : env[3])

test(5)
test2()
test3()


