#$Id: observationSimLib.py,v 1.1.1.1 2013/02/11 21:30:45 areustle Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['observationSim'])
    env.Tool('facilitiesLib')
    env.Tool('tipLib')
    env.Tool('astroLib')
    env.Tool('fluxLib')
    env.Tool('st_facilitiesLib')
    env.Tool('celestialSourcesLib')
    env.Tool('irfsLib')
    env.Tool('dataSubselectorLib')
    env.Tool('fitsGenLib')

def exists(env):
    return 1
