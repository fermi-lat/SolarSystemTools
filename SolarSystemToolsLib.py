# $Header: $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['SolarSystemTools'])
    env.Tool('tipLib')
    env.Tool('astroLib')
    env.Tool('map_toolsLib')
    env.Tool('healpixLib')
    env.Tool('st_facilitiesLib')
    env.Tool('facilitiesLib')
    env.Tool('st_appLib')
    env.Tool('irfLoaderLib')
    env.Tool('evtbinLib')
    env.Tool('LikelihoodLib')

def exists(env):
    return 1
