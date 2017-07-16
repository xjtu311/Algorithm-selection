#!/usr/bin/env python
import sys
import time
import os
import socket


def f1(m):
  return m['x'] + m['y']
  
def f2(m):
  return 1 - m['x'] + m['y']

def f3(m):
  return m['x'] + 1-m['y']

def f4(m):
  return 1-m['x'] + 1-m['y']

def f5(m):
  return 5 - m['x'] + m['y']

def f6(m):
  return m['x'] + 5-m['y']

def f7(m):
  return 5-m['x'] + 5-m['y']


inst_func = {'x+y' : f1,
             '(1-x)+y': f2,
             'x+(1-y)': f3,
             '(1-x)+(1-y)':f4,
             '(5-x)+y': f5,
             'x+(5-y)': f6,
             '(5-x)+(5-y)':f7}


inst = sys.argv[1]



isi = sys.argv[2]
cutoff = float(sys.argv[3])
cutoff_length = sys.argv[4]
seed = sys.argv[5]
i = 6;
values = {}
while i < len(sys.argv):
  param_name = sys.argv[i][1:]
  param_value = sys.argv[i+1]
  
  values[param_name] = float(param_value)
  i = i + 2



print (values)

port = -1
sock = {}
if 'ACLIB_PORT' in os.environ:
  port = os.environ['ACLIB_PORT']
  sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  
  
runtime = inst_func[inst](values) + values['z']
crun = 0;

runtime = runtime + 0.1


timeToSleep = min(runtime, float(cutoff))
while(crun < timeToSleep):
  crun = crun+0.125
  time.sleep(0.03)
  if port > 0:
      sock.sendto(str(crun), ("127.0.0.1",int(port)))
    

#print runtime > cutoff
#print cutoff
quality=runtime
if runtime > cutoff:
  print ("Result for ParamILS: TIMEOUT, %f, 0, %f, %s"%(cutoff,cutoff,seed))
else:
  print ("Result for ParamILS: SAT, %f, 0, %f, %s"%(runtime,quality,seed))


