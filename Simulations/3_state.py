#!/share/apps/python-2.7.2/bin/python
#this program take an integer as input so that you can run a bunch of them with a script

# Imports
import sys
import argparse
import random, math, os 
from time import time

def main():
    for arg in sys.argv:
	print arg

    #define the landscape to simulate over
    States = ("A","B","C")
    Scores = (1,1.15,1.25)

    #Number of mutations to accept
    max_accept_mut = 100000
    #Population size
    N = 10
  
    #set staring position
    current_state=States[0]
    current_score=Scores[0]

    #Prepare data headers
    data = ['Variant,Score,ScoreDiff,Probability,Generation\n']
    #set up file to write output to
    data.append(str(current_state)+ "," + str(current_score) + ",0.0,0.0,0.0\n")

    #start evolution
    i=0
    gen=0
    while i < max_accept_mut:
            #update the number of generations that have pased
            gen+=1

	    #choose the new state to mut to
	    new_mut_key = random.randint(0,2)
	    proposed_res = States[new_mut_key]
	    while proposed_res == current_state:
		new_mut_key = random.randint(0,2)
		proposed_res = States[new_mut_key]

	    #score mutant
	    variant_score = Scores[new_mut_key]

	    #get the probability that the mutation will be accepted
	    probability = calc_prob_mh(variant_score, current_score, N, beta, thresholds)

	    #create a name for the mutant if its going to be kept 
	    variant_name = proposed_res

	    
	    #test to see if mutation is accepted
	    if random.random() < probability:
		print 'accepts:', i 
		#update the wildtype 
		current_state = proposed_res
		current_score = variant_score

	    #save name and energy change
	    data.append(current_state + "," + str(current_score) + "," + str(variant_score - current_score) + "," + str(probability) + "," + str(gen) + "\n")

            #update number of moves
	    i+=1
		

    print '\nMutations and scoring complete.'
    t1 = time()
    # Output results.
    data_filename = 'MH_3_state_'+str(arg)+'.csv'
    with open(data_filename, "w") as f:
        f.writelines(data)

    print 'Data written to:', data_filename
   


###assorted functions that have to do with scoring and prob of acceptance ####


#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):

  xi=stab_org
  xj=stab_mut

  if xi==xj:
	return(1.0)

  if xj > xi:
    return((1.0))
  else:
    ans=pow(float(xj/xi),2*float(N)-2)
    return(ans)


#score functions for met-hastings selection
def calc_prob_fix(stab_mut, stab_org, N, beta, thresholds):

  xi=stab_org
  xj=stab_mut

  if xi == xj:
    return(1.0/float(N))

  if(xj==0.0):
    return(0.0)

  try:
    p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
  except OverflowError as e:
    p = 0.0
  return (p)


#Run main program
if __name__ == '__main__':
   main()
