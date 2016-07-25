import collections
import sys
c = collections.Counter(int(i.split()[-1]) for i in sys.stdin)
for k,v in sorted(c.iteritems()):
  print '{0}\t{1}'.format(k, v)
