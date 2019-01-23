import numpy as np
from deap import algorithms, base, creator, tools
from sklearn.metrics import r2_score
from data import aTpdaMpiGA as aTpdaMpi
import random,numpy

P_SIZE=int(np.ceil(np.log2(1000)))
def gaussian(x, mu, delta):
    exp = np.exp(- np.power(x - mu, 2) / (2 * np.power(delta, 2)))
    c = 1 / (delta * np.sqrt(2 * np.pi))
    return c * exp

def mixGaussian(x, mu, delta, alpha,split=False):
    if split:
        r=[]
    else:
        r=0
    
    for m,v,p in zip(mu,delta,alpha):
        t=p*gaussian(x,m,v)
        if split:
            r.append(t)
        else:
            r=r+t
    return r

class gaMixGauss:
    def __init__(self,x,y,mu,delta,snum=2):
        self.snum=snum
        self.x=x
        self.y=y
        self.mu=mu
        self.delta=delta        
        self.bitnum=P_SIZE*snum
        self.maxy=max(y)*2

    def evalMixGaussian(self,params):
        alpha=self.translateDNA(params)
        mg=mixGaussian(self.x, self.mu, self.delta, alpha)
        # print("len(mg)",len(mg)) numpy.sqrt(len(mg))*r2_score(self.y,mg)-
        return (r2_score(self.y,mg)-abs(sum(mg)-sum(self.y)),)#

    def translateDNA(self,pop):
        r=[]
        for i in range(self.snum):
            idx_bit0=P_SIZE*i
            idx_bit1=idx_bit0+P_SIZE
            # print(pop(idx_bit0,idx_bit1))
            p1=np.asarray(pop[idx_bit0:idx_bit1])
            r.append (p1.dot(2 ** np.arange(P_SIZE)[::-1]) / float(2**P_SIZE-1))#*self.maxy)    
        return r    
    def gaAlpha(self):        
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)
        toolbox = base.Toolbox()
        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("avg", numpy.mean)
        stats.register("std", numpy.std)
        stats.register("min", numpy.min)
        stats.register("max", numpy.max)        
        toolbox.register("attr_bool", random.randint, 0, 1)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=self.bitnum)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", self.evalMixGaussian)
        toolbox.register("mate", tools.cxUniform, indpb=0.5)
        toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
        # toolbox.register("select", tools.selTournament, tournsize=5)
        toolbox.register("select", tools.selRandom)
        pop = toolbox.population(n=self.bitnum*10)
        algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.05, ngen=600, verbose=False)
        res=tools.selBest(pop, k=20)
        for r in res:
            print(self.translateDNA(r))
        print(stats.compile(res))
        print(res[0].fitness.values)
        return res