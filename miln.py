from __future__ import division
from pylab import *
import numpy as np
import time as tm
from random import random
import datetime
import scipy
from scipy.integrate import odeint

from pylab import *
import time
from random import random
import datetime


def yLeadUp(s):
    return s.replace("^", "**")


def yRec(s):
    s = yLeadUp(s)
    t = 'lambda x, y: ' + s
    f = eval(t)
    return f


try:
    f = open("input.txt")
except:
    print "Can`t open file"
    exit()
try:
    text = f.readlines()
except:
    print "Can`t read from file"
    exit()
try:
    funY = yRec(text[0])
    a, b, = [float(x) for x in text[3].split(',')[:-1]]
    n = int(text[3].split(',')[-1])
except:
    print "Can`r recognize data"
    exit()

startM = time.time()
X = linspace(a, b, n)
h = X[1] - X[0]
"""funY = lambda x, y: 1/x"""
trueYVals = []
trueYVals.append(0)
eps = 0.05
d = 0

""" First solution for Miln using Runge-Kutta """
for i in range(n):
    eta1 = funY(X[i], trueYVals[i])
    eta2 = funY(X[i] + h / 2, trueYVals[i] + h / 2)
    eta3 = funY(X[i] + h / 2, trueYVals[i] + h * eta2 / 2)
    eta4 = funY(X[i] + h / 2, trueYVals[i] + h * eta3 / 2)
    tmp = h * (eta1 + 2 * eta2 + 2 * eta3 + eta4) / 6
    trueYVals.append(trueYVals[i] + tmp)

yVals = trueYVals[:4]
ydVals = trueYVals[:4]
yValsM = trueYVals[:4]

""" Getting solution using Miln-method """

for i in xrange(5, n):
    if abs(d) < eps:
        ydVals.append(yVals[i - 3] + (3 * h) / 2 * (
            7 * funY(X[i - 3], yVals[i - 3]) - funY(X[i - 2], yVals[i - 2]) + 2 * funY(X[i - 1], yVals[i - 1])))

        yVals.append((yVals[i - 2] + h / 3 * (
            funY(X[i], ydVals[i]) - 4 * funY(X[i - 1], yVals[i - 1]) + funY(X[i - 2], yVals[i - 2]))))
        d = (ydVals[i] - yVals[i]) / 29
    else:
        h *= 0.25
        ydVals.append(yVals[i - 3] + (3 * h) / 2 * (
            7 * funY(X[i - 3], yVals[i - 3]) - funY(X[i - 2], yVals[i - 2]) + 2 * funY(X[i - 1], yVals[i - 1])))

        yVals.append((yVals[i - 2] + h / 3 * (
            funY(X[i], ydVals[i]) - 4 * funY(X[i - 1], yVals[i - 1]) + funY(X[i - 2], yVals[i - 2]))))
        d = (ydVals[i] - yVals[i]) / 29

finalM = time.time()

""" Getting started the Adams """
startA = time.time()
yValsA = trueYVals[:4]
yBValsA = trueYVals[:4]

for i in xrange(4, n):
    if abs(d) < eps:
        yBValsA.append(yValsA[i - 1] + h / 12 * (
            23 * funY(X[i - 1], yValsA[i - 1]) + 16 * funY(X[i - 2], yValsA[i - 2]) + 5 * funY(X[i - 3],
                                                                                               yValsA[i - 3])))
        yValsA.append(yValsA[i - 1] + h / 12 * (
            5 * funY(X[i - 1], yBValsA[i - 1]) + 8 * funY(X[i - 1], yValsA[i - 1]) + funY(X[i - 2], yValsA[i - 2])))
        yValsM.append(yValsA[i - 1] + (h / 24) * (
            9 * funY(X[i], yValsA[i]) + 19 * funY(X[i - 1], yValsA[i - 1]) - 5 * funY(X[i - 2], yValsA[i - 2]) + funY(
                X[i - 3], yValsA[i - 3]),
            yValsA[i - 3]))
        d = (19 / 720) * (yValsM[i] - yBValsA[i])
        yVals[i] *= -1
    else:
        h *= 0.25
        yBValsA.append(yValsA[i - 1] + h / 12 * (
            23 * funY(X[i - 1], yValsA[i - 1]) + 16 * funY(X[i - 2], yValsA[i - 2]) + 5 * funY(X[i - 3],
                                                                                               yValsA[i - 3])))
        yValsA.append(yValsA[i - 1] + h / 12 * (
            5 * funY(X[i - 1], yBValsA[i - 1]) + 8 * funY(X[i - 1], yValsA[i - 1]) + funY(X[i - 2], yValsA[i - 2])))
        yValsM.append(yValsA[i - 1] + (h / 24) * (
            9 * funY(X[i], yValsA[i]) + 19 * funY(X[i - 1], yValsA[i - 1]) - 5 * funY(X[i - 2], yValsA[i - 2]) + funY(
                X[i - 3], yValsA[i - 3])))
        d = (19 / 720) * (yValsM[i] - yBValsA[i])
        yValsA[i] *= -1

finishA = time.time()

y0 = 0
y = []
y = odeint(funY, y0, X)
print y
print "timeM: %0.9f s" % (finalM - startM)
print "timeA: %0.9f s" % (finishA - startA)

plt.plot(X, yVals, 'g')
plt.plot(X, yValsA, 'b')
plt.plot(X, y, 'r')

try:
    savefig(str(datetime.datetime.now()).replace(':', '-')[:-7] + ' lab3.png')
except:
    print "Can`t save image"
show()















