import pandas as pd
import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import chart_studio.plotly as py
import plotly.figure_factory as ff
import datetime

class Genetic:
    
    def __init__(self):
        # Initial基本參數
        self.__parameters = {
            "populationSize": 30,
            "crossoverRate": 0.8,
            "mutationRate": 0.2,
            "mutationSelectionRate": 0.2,
            "maxGeneration": 2000,
        }
    
    def initParameters(self, pt_temp, ms_temp):
        print("Calculating...")
        
        self.__dfshape=pt_temp.shape
        self.__num_mc=self.__dfshape[1] # number of machines
        self.__num_job=self.__dfshape[0] # number of jobs
        self.__num_gene=self.__num_mc*self.__num_job # number of genes in a chromosome
        
        self.__pt=[list(map(int, pt_temp.iloc[i])) for i in range(self.__num_job)]
        self.__ms=[list(map(int, ms_temp.iloc[i])) for i in range(self.__num_job)]    


        # raw_input is used in python 2
        self.__parameters["populationSize"]=int(input('Please input the size of population: ') or 30) # default value is 30
        self.__parameters["crossoverRate"]=float(input('Please input the size of Crossover Rate: ') or 0.8) # default value is 0.8
        self.__parameters["mutationRate"]=float(input('Please input the size of Mutation Rate: ') or 0.2) # default value is 0.2
        self.__parameters["mutationSelectionRate"]=float(input('Please input the mutation selection rate: ') or 0.2)
        self.__num_mutation_jobs=round(self.__num_gene*self.__parameters["mutationSelectionRate"])
        self.__parameters["maxGeneration"]=int(input('Please input number of iteration: ') or 2000) # default value is 2000
         
    def newPopulation(self, i):
        # Generate first population
        self.__nxm_random_num=list(np.random.permutation(self.__num_gene)) # generate a random permutation of 0 to self.__num_job*num_mc-1
        self.__population_list.append(self.__nxm_random_num)
        for j in range(self.__num_gene):
            self.__population_list[i][j]=self.__population_list[i][j]%self.__num_job # convert to job number format, every job appears m times
       
    def reproduction(self):
        # Generate new groups
        self.__Tbest_now=99999999999
        #--- two point crossover ---
        self.__parent_list=copy.deepcopy(self.__population_list)
        self.__offspring_list=copy.deepcopy(self.__population_list)
        S=list(np.random.permutation(self.__parameters["populationSize"])) # generate a random sequence to select the parent chromosome to crossover

        for m in range(int(self.__parameters["populationSize"]/2)):
            crossover_prob=np.random.rand()
            if(crossover_prob <= self.__parameters["crossoverRate"]):
                parent_1= self.__population_list[S[2*m]][:]
                parent_2= self.__population_list[S[2*m+1]][:]
                child_1=parent_1[:]
                child_2=parent_2[:]
                cutpoint=list(np.random.choice(self.__num_gene, 2, replace=False))
                cutpoint.sort()

                child_1[cutpoint[0]:cutpoint[1]]=parent_2[cutpoint[0]:cutpoint[1]]
                child_2[cutpoint[0]:cutpoint[1]]=parent_1[cutpoint[0]:cutpoint[1]]
                self.__offspring_list[S[2*m]]=child_1[:]
                self.__offspring_list[S[2*m+1]]=child_2[:]      
   
    def repair(self):
        # Repair
        for m in range(self.__parameters["populationSize"]):
            job_count={}
            larger,less=[],[] # 'larger' record jobs appear in the chromosome more than m times, and 'less' records less than m times.
            for i in range(self.__num_job):
                if i in self.__offspring_list[m]:
                    count=self.__offspring_list[m].count(i)
                    pos=self.__offspring_list[m].index(i)
                    job_count[i]=[count,pos] # store the above two values to the job_count dictionary
                else:
                    count=0
                    job_count[i]=[count,0]
                if count>self.__num_mc:
                    larger.append(i)
                elif count<self.__num_mc:
                    less.append(i)

            for k in range(len(larger)):
                chg_job=larger[k]
                while job_count[chg_job][0] > self.__num_mc:
                    for d in range(len(less)):
                        if job_count[less[d]][0] < self.__num_mc:                    
                            self.__offspring_list[m][job_count[chg_job][1]]=less[d]
                            job_count[chg_job][1]=self.__offspring_list[m].index(chg_job)
                            job_count[chg_job][0]=job_count[chg_job][0]-1
                            job_count[less[d]][0]=job_count[less[d]][0]+1                    
                        if job_count[chg_job][0]==self.__num_mc:
                            break
                        
    def mutation(self):
        # Mutation
        for m in range(len(self.__offspring_list)):
            mutation_prob=np.random.rand()
            if(mutation_prob <= self.__parameters["mutationRate"]):
                m_chg=list(np.random.choice(self.__num_gene, self.__num_mutation_jobs, replace=False)) # chooses the position to mutation
                t_value_last=self.__offspring_list[m][m_chg[0]] # save the value which is on the first mutation position
                for i in range(self.__num_mutation_jobs-1):
                    self.__offspring_list[m][m_chg[i]]=self.__offspring_list[m][m_chg[i+1]] # displacement

                self.__offspring_list[m][m_chg[self.__num_mutation_jobs-1]]=t_value_last # move the value of the first mutation position to the last mutation position  
                
    def fitness(self):
        # Calculate fitness (fitness value(calculate makespan))
        self.__total_chromosome=copy.deepcopy(self.__parent_list)+copy.deepcopy(self.__offspring_list) # parent and offspring chromosomes combination
        self.__chrom_fitness,self.__chrom_fit=[],[]
        self.__total_fitness=0
        for m in range(self.__parameters["populationSize"]*2):
            j_keys=[j for j in range(self.__num_job)]
            key_count={key:0 for key in j_keys}
            j_count={key:0 for key in j_keys}
            m_keys=[j+1 for j in range(self.__num_mc)]
            m_count={key:0 for key in m_keys}

            for i in self.__total_chromosome[m]:
                gen_t=int(self.__pt[i][key_count[i]])
                gen_m=int(self.__ms[i][key_count[i]])
                j_count[i]=j_count[i]+gen_t
                m_count[gen_m]=m_count[gen_m]+gen_t

                if m_count[gen_m]<j_count[i]:
                    m_count[gen_m]=j_count[i]
                elif m_count[gen_m]>j_count[i]:
                    j_count[i]=m_count[gen_m]

                key_count[i]=key_count[i]+1

            makespan=max(j_count.values())
            self.__chrom_fitness.append(1/makespan)
            self.__chrom_fit.append(makespan)
            self.__total_fitness=self.__total_fitness+self.__chrom_fitness[m]         
            
    def selection(self):
        # Select (roulette wheel approach)
        pk,qk=[],[]
    
        for i in range(self.__parameters["populationSize"]*2):
            pk.append(self.__chrom_fitness[i]/self.__total_fitness)
        for i in range(self.__parameters["populationSize"]*2):
            cumulative=0
            for j in range(0,i+1):
                cumulative=cumulative+pk[j]
            qk.append(cumulative)

        selection_rand=[np.random.rand() for i in range(self.__parameters["populationSize"])]

        for i in range(self.__parameters["populationSize"]):
            if selection_rand[i]<=qk[0]:
                self.__population_list[i]=copy.deepcopy(self.__total_chromosome[0])
            else:
                for j in range(0,self.__parameters["populationSize"]*2-1):
                    if selection_rand[i]>qk[j] and selection_rand[i]<=qk[j+1]:
                        self.__population_list[i]=copy.deepcopy(self.__total_chromosome[j+1])
                        break
                    
    def compare(self):
        # Compare
        for i in range(self.__parameters["populationSize"]*2):
            if self.__chrom_fit[i]<self.__Tbest_now:
                self.__Tbest_now=self.__chrom_fit[i]
                sequence_now=copy.deepcopy(self.__total_chromosome[i])
        if self.__Tbest_now<=self.__Tbest:
            self.__Tbest=self.__Tbest_now
            self.__sequence_best=copy.deepcopy(sequence_now)
        
        self.__makespan_record.append(self.__Tbest)    
    
    def runGeneticAlgorithm(self):
        self.__Tbest=999999999999999
        self.__best_list, self.__best_obj=[],[]
        self.__population_list=[]
        self.__makespan_record=[]
        
        for i in range(self.__parameters["populationSize"]):
            self.newPopulation(i)
            
        for n in range(self.__parameters["maxGeneration"]):
            print("Generation:" + str(n), end='\r')
            self.reproduction()
            self.repair()  
            self.mutation()
            self.fitness()
            self.selection()
            self.compare()
    
    def printResult(self):
        # Print the calculation results
        print("optimal sequence",self.__sequence_best)
        print("optimal value:%f"%self.__Tbest) 
        
    def drawFitnessPlot(self):
        # Fitness (makespan) Plot
        plt.plot([i for i in range(len(self.__makespan_record))],self.__makespan_record,'b')
        plt.ylabel('makespan',fontsize=15)
        plt.xlabel('generation',fontsize=15)
        plt.show() 
        
    def drawGanttChart(self):
        # Gantt Chart
        m_keys=[j+1 for j in range(self.__num_mc)]
        j_keys=[j for j in range(self.__num_job)]
        key_count={key:0 for key in j_keys}
        j_count={key:0 for key in j_keys}
        m_count={key:0 for key in m_keys}
        j_record={}
        for i in self.__sequence_best:
            gen_t=int(self.__pt[i][key_count[i]])
            gen_m=int(self.__ms[i][key_count[i]])
            j_count[i]=j_count[i]+gen_t
            m_count[gen_m]=m_count[gen_m]+gen_t

            if m_count[gen_m]<j_count[i]:
                m_count[gen_m]=j_count[i]
            elif m_count[gen_m]>j_count[i]:
                j_count[i]=m_count[gen_m]

            start_time=str(datetime.timedelta(seconds=j_count[i]-self.__pt[i][key_count[i]])) # convert seconds to hours, minutes and seconds
            end_time=str(datetime.timedelta(seconds=j_count[i]))

            j_record[(i,gen_m)]=[start_time,end_time]

            key_count[i]=key_count[i]+1


        df=[]
        for m in m_keys:
            for j in j_keys:
                df.append(dict(Task='Machine %s'%(m), Start='2018-07-14 %s'%(str(j_record[(j,m)][0])), Finish='2018-07-14 %s'%(str(j_record[(j,m)][1])),Resource='Job %s'%(j+1)))
                
        fig = ff.create_gantt(df, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True, title='Job shop Schedule')
        py.iplot(fig, filename='GA_job_shop_scheduling', world_readable=True, auto_open=True) 
        #fig.show()
        
          
                    
                      
                        