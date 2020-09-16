#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 23:51:31 2017

@author: rodrigodelafuente
"""

import simpy
import time

def clock(env, name, tick):
    while True:
        time.sleep(10)
        print(name, env.now)
        yield env.timeout(tick)

env = simpy.Environment()
env.process(clock(env, 'fast', 0.5))
env.process(clock(env, 'slow', 1))
env.run(until=2.00001)
