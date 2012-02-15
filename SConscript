# -*- python -*-
#
# Authors: Gudlaugur Johannesson <gudlaugu@glast2.stanford.edu>
# Version: SolarSystemTools-00-00-00
# $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SConscript,v 1.2 2012/02/13 22:24:36 gudlaugu Exp $
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

#libEnv.Tool('addLinkDeps', package='SolarSystemTools', toBuild='shared')
SolarSystemToolsLib = libEnv.StaticLibrary('SolarSystemTools', 
                                          listFiles(['src/*.cxx']))

progEnv.Tool('SolarSystemToolsLib')

#testEnv = progEnv.Clone()
#testEnv.Tool('addLibrary', library = baseEnv['cppunitLibs'])
#test_SolarSystemToolsBin = testEnv.Program('test_SolarSystemTools', 
#                                          listFiles(['src/test/*.cxx']))

gtltcubesun = progEnv.Program('gtltcubesun', 
                              listFiles(['src/makeSolarExposureCube/*.cxx']))

gtexpcubesun = progEnv.Program('gtexpcubesun',
                              listFiles(['src/gtexpcubesun/*.cxx']))

gtexphpsun = progEnv.Program('gtexphpsun',
                              listFiles(['src/gtexphpsun/*.cxx']))

gtsuntemp = progEnv.Program('gtsuntemp',
                              listFiles(['src/gtsuntemp/*.cxx']))

progEnv.Tool('registerTargets', package = 'SolarSystemTools', 
             staticLibraryCxts = [[SolarSystemToolsLib, libEnv]],
             binaryCxts = [[gtltcubesun, progEnv], 
                           [gtexpcubesun, progEnv],
                           [gtexphpsun, progEnv],
                           [gtsuntemp, progEnv]],
#             testAppCxts = [[test_SolarSystemToolsBin, testEnv]],
             includes = listFiles(['SolarSystemTools/*.h']),
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))
