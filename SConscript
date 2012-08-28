# -*- python -*-
#
# Authors: Gudlaugur Johannesson <gudlaugu@glast2.stanford.edu>
# Version: SolarSystemTools-01-01-00
# $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SConscript,v 1.10 2012/06/16 11:31:10 gudlaugu Exp $
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

gtltsumsun = progEnv.Program('gtltsumsun', 
                              listFiles(['src/gtltsumsun/*.cxx']))

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
                           [gtsuntemp, progEnv],
                           [gtltsumsun, progEnv]],
#             testAppCxts = [[test_SolarSystemToolsBin, testEnv]],
             includes = listFiles(['SolarSystemTools/*.h']),
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))
