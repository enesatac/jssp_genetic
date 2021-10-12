'''==========Solving job shop scheduling problem by genetic algorithm in python======='''
# importing required modules
import pandas as pd
import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import chart_studio.plotly as py
import plotly.figure_factory as ff
import datetime

#import genetic as ga
import genetic as ga

''' ================= initialization setting ======================'''

genetic = ga.Genetic()

pt_tmp=pd.read_excel("jssp_dataset.xlsx",sheet_name="Processing Time",index_col =[0])
ms_tmp=pd.read_excel("jssp_dataset.xlsx",sheet_name="Machines Sequence",index_col =[0])

genetic.initParameters(pt_tmp, ms_tmp)
    
'''==================== main code ==============================='''

start_time = time.time()

'''----- generate initial population -----'''
genetic.runGeneticAlgorithm()
        
end_time = time.time()    
    
genetic.printResult()

print('the elapsed time:%s'% (end_time - start_time))

genetic.drawFitnessPlot()
genetic.drawGanttChart()

