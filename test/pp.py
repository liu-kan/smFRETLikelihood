import multiprocessing
import os,sys
import random
import time

class simulation(multiprocessing.Process):
    def __init__(self, name):
        # must call this before anything else
        multiprocessing.Process.__init__(self)

        # then any other initialization
        self.name = name
        self.number = 0.0
        sys.stdout.write('[%s] created: %f\n' % (self.name, self.number))

    def run(self):
        sys.stdout.write('[%s] running ...  process id: %s\n' 
                         % (self.name, os.getpid()))

        self.number = random.uniform(0.0, 10.0)
        time.sleep( 5 )
        sys.stdout.write('[%s] completed: %f\n' % (self.name, self.number))
        

if __name__ == '__main__':
    multiprocessing.freeze_support()
    sim_list = []
    sim_list.append(simulation('foo'))
    sim_list.append(simulation('bar'))
    for sim in sim_list:
        sim.start()
    for sim in sim_list:
        sim.join()
        print("joined:",sim.pid," exitcode:",sim.exitcode)