print "input string for help"
helpstr = raw_input()

for i in dir():
   if helpstr in i:
     print i
