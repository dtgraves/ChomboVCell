def GetNumberOfLevels():
   return len(SILRestriction().SetsInCategory('levels'))

def TurnOnLevel(level):
   sil = SILRestriction()
   levelSetId = sil.SetsInCategory('levels')[level]
   sil.TurnOnSet(levelSetId)
   sil.TurnOffSet(sil.SetsInCategory('materials')[1])
   SetPlotSILRestriction(sil)

def TurnOffLevel(level):
   sil = SILRestriction()
   levelSetId = sil.SetsInCategory('levels')[level]
   sil.TurnOffSet(levelSetId)
   sil.TurnOffSet(sil.SetsInCategory('materials')[1])
   SetPlotSILRestriction(sil)

def TurnOnEverything():
   sil = SILRestriction()
   sil.TurnOnAll()
   SetPlotSILRestriction(sil)

def TurnOnAllLevels():
   sil = SILRestriction()
   sil.TurnOffAll()
   levelSetIds = sil.SetsInCategory('levels')
   for levelId in levelSetIds:
      sil.TurnOnSet(levelId)
   sil.TurnOffSet(sil.SetsInCategory('materials')[1])
   SetPlotSILRestriction(sil)

#DeleteAllPlots()

import sys
print "input variable number"
strvar = sys.stdin.readline()
ivar = int(strvar)
print "input time step number"
strstep = sys.stdin.readline()
step = int(strstep)
print "input the number of volumes"
strnvol = sys.stdin.readline()
nvol = int(strnvol)
print "input SpaceDim"
strsdim = sys.stdin.readline()
sdim = int(strsdim)
print 'time step number = %d, variable number = %d, number of volumes = %d'%(step, ivar,nvol)
ivol = 0
ivolmin = -1
ivolmax = -1
globalmax =  -1.0e30
globalmin =  1.0e30
while ivol < nvol:
  print "doing volume %d"%(ivol)
  strfil = "./mitochondria.00000%d.vol%d.%dd.hdf5"%(step, ivol, sdim)
  if step >= 10:
      strfil = "./mitochondria.0000%d.vol%d.%dd.hdf5"%(step, ivol, sdim)
  if step >= 100:
      strfil = "./mitochondria.000%d.vol%d.%dd.hdf5"%(step, ivol, sdim)
  if step >= 1000:
      strfil = "./mitochondria.00%d.vol%d.%dd.hdf5"%(step, ivol, sdim)
  if step >= 10000:
      strfil = "./mitochondria.0%d.vol%d.%dd.hdf5"%(step, ivol, sdim)
  if step >= 100000:
      strfil = "./mitochondria.%d.vol%d.%dd.hdf5"%(step, ivol, sdim)

  strvar = "G.var%d"%(ivar)
  print "file = "+strfil
  print "variable = "+strvar
  OpenDatabase(strfil)
  AddPlot("Pseudocolor", strvar)
  TurnOnAllLevels()
  DrawPlots()
  Query("MinMax")
  (localmin, localmax) = GetQueryOutputValue()
  if localmax > globalmax:
      globalmax = localmax
      ivolmax = ivol
  if localmin < globalmin:
      globalmin = localmin
      ivolmin = ivol
  print "local max = %e, local min = %e"%(localmax, localmin)
  ivol = ivol + 1


p = PseudocolorAttributes()
p.legendFlag = 0
p.max = globalmax
p.min = globalmin
p.maxFlag = 1
p.minFlag = 1
p.colorTableName = "hot_desaturated"
SetActivePlots(tuple(range(GetNumPlots())))
SetPlotOptions(p)


print "global max = %e, global min = %e, ivolmax = %d, ivolmin = %d"%(globalmax, globalmin, ivolmax, ivolmin)

#a = GetAnnotationAttributes()
#
#a.userInfoFlag = 0
#a.databaseInfoFlag = 0
#a.legendInfoFlag = 0
#SetAnnotationAttributes(a)

#OpenGUI() opens the graphical user interface



#for i in dir():
#...   if "Plot" in i:
#...     print i
