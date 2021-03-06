#!/share/apps/python-2.7.2/bin/python
# Calculates the probability of accepting a reversion to the ancesteral state over the next 15 accepted mutations, run this in a directory with all of the pdb files that score_mutant_pdb makes.
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 

# Imports
import sys
import argparse
import random, math, os, glob 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

from toolbox import mutate_residue
from toolbox import cleanATOM
from time import time
import re

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


def main():
    # Initialize Rosetta.
    init(extra_options='-mute basic -mute core')

    # Constants
    PACK_RADIUS = 10.0
    #Population size
    N = 37
    #Beta (temp term)
    beta = 1
    #look up what the first stored value was in the files to get the threshold
    threshold = float(-534.687360627/2)

    #Set up ScoreFunction
    sf = get_fa_scorefxn()

    #Set up MoveMap.
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    #Prepare data headers
    data = ['Generation,RevertTo,OrgScore,RevScore,Change,Prob\n']

    # Get the reversions file, the output file the score_mutant_pdb has made
    variant_scores=open('mh_rep_3_37.csv')

    #get just the mutation we want to revert to
    lines= variant_scores.readlines()
    var_line=lines[2] #gets the Nth line how ever long you want the burn to be
    var_line=var_line.split(',')[0]
  
    var_loc=int(filter(str.isdigit, var_line))
    var_rev=var_line[:1]

    gen=1
    #get all the pdb files
    sort_list=sorted(glob.glob('*.pdb'), key=numericalSort)
 
    for i in range(1,len(sort_list)-15):
     
      #calc reversion for next 15 moves
      for infile in sorted(glob.glob('*.pdb'), key=numericalSort)[i:i+15]:

	#for each mutation	
        var_line=lines[gen+1] #gets the Nth line how ever long you want the burn to be
        var_line=var_line.split(',')[0]
        var_loc=int(filter(str.isdigit, var_line))
        var_rev=var_line[:1]

      	print "Current File Being Processed is: " + infile
        initial_pose = pose_from_pdb(infile)
        initial_score = sf(initial_pose)
	print("init scored")
        mutant_pose = mutate_residue(initial_pose, var_loc , var_rev, PACK_RADIUS, sf)
        variant_score = sf(mutant_pose)
        probability = calc_prob_mh(variant_score, initial_score, N, beta, threshold)
	print(str(gen)+","+ var_line + "," +str(initial_score) + "," + str(variant_score) + "," + str(variant_score - initial_score)+ ","+ str(probability)+ "\n")
      	data.append(str(gen)+","+ var_line + "," +str(initial_score) + "," + str(variant_score) + "," + str(variant_score - initial_score)+ ","+ str(probability)+ "\n")
      gen+=1

    print '\nDONE'

    data_filename = 'rep_3_mh_37_rev_15_score.csv'
    with open(data_filename, "w") as f:
        f.writelines(data)



#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xj > xi:
    return((1.0))
  else:
    #change if you want to use exp notation
    #exponent = -2 * float(N) * (xi - xj)
    ans=pow(float(xj/xi),2*float(N)-2)
    return(ans)
    #return(safe_calc(exponent))



#score functions for met-hastings selection
def calc_prob_fix(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xi == xj:
    return(1/float(N))

  if(xj==0.0):
    return(0.0)

  try:
    p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
  except OverflowError as e:
    p = 0.0
  return (p)

#cal fitness either in logs or not
def calc_x_fix(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total = 1/(safe_calc(exponent) + 1)
  return(total)

def calc_x(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total += -math.log(safe_calc(exponent) + 1)
  return(total)

def safe_calc(exponent):
  if exponent > 700:
    print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))

#Run main program
if __name__ == '__main__':
   main()
