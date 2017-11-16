import geneticalgo
import random
import copy
import os
import json
import sys
import math
import accordionutils as accordion
from geneticalgo import algorithm
from geneticalgo import primitives
from geneticalgo import shapewriter
from sofalauncher import launcher



individualId = 0
heightTube = 25.0
radiusTube = 3.0
thickness = 1.0
number_of_cavities = 3
generate_random="OFF"
type=["frisbee","ellipsoid"]#

Vref = 1033.0

mutationType="ON"
mutationAxisX="OFF"
mutationAxisY="OFF"
mutationAxisZ="ON"
mutationRotation="OFF"

def getNextId():
    global individualId
    individualId+=1
    return individualId

def getJSONFragmentFrom(file):

    with open(file, "r") as source:
        for line in source:
            if "JSON : " in line:
                data = json.loads(line[7:])
                return data

    return None




class AccordionIndividual(algorithm.Individual):

    def __init__(self):
        self.level=None
        self.id = getNextId()
        global heightTube, radiusTube, thickness

        accordion.create(self, heightTube, radiusTube, thickness)
        
    def display(self):
        temp="heightTube="+str(self.height)+"\n"\
             +"radiusTube="+str(self.radius)+"\n"\
             +"thickness="+str(self.thickness)+"\n\n"\

        for i in range(len(self.listCavities)):
            [height,type,axisX,axisY,axisZ]=self.listCavities[i]
            temp=temp+"cavity number "+str(i)+"\n"\
                     +"type="+type+"\n"\
                     +"height="+str(height)+"\n"\
                     +"axisX="+str(axisX)+"\n"\
                     +"axisY="+str(axisY)+"\n"\
                     +"axisZ="+str(axisZ)+"\n\n"
        return temp


def newIndividualFrom(ind):
    newInd=AccordionIndividual()
    newInd.listCavities=copy.deepcopy(ind.listCavities)

    return newInd

def getShapeFromInd(ind):
    return accordion.accordionFreeDimension(ind, ind.height, ind.radius, ind.thickness, ind.listCavities)


#def copyInd (individual):
#    temp=Individual(individual.height,individual.radius,individual.thickness)
#    listCavities=copy.deepcopy(individual.listCavities)
#    temp.listCavities=listCavities
#    return temp


####################################################################################################
### MUTATION
###
def mutation_axisX(ind):

    length=len(ind.listCavities)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    index=random.randint(1,length-2)

    axisX=ind.listCavities[index][2]
    epsilon=random.uniform(-0.5,0.5)
    axisX=min(axisX+epsilon, (15.0/3.0)*ind.radius)

    axisX = max(axisX, max(2.0*ind.radius, ind.thickness+1.0))

    ind.listCavities[index][2]=axisX

def mutation_axisY(ind):

    length=len(ind.listCavities)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    index=random.randint(1,length-2)
    axisY=ind.listCavities[index][3]
    epsilon=random.uniform(-0.5,0.5)
    axisY=min(axisY+epsilon, (15.0/3.0)*ind.radius)
    axisY = max(axisY, max(2.0*ind.radius, ind.thickness+1.0))
    ind.listCavities[index][3]=axisY


def mutation_axisZ(ind):
    length=len(ind.listCavities)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    index=random.randint(1,length-2)
    axisZ=ind.listCavities[index][4]
    cavity = ind.listCavities[index]
    epsilon=random.uniform(-0.5,0.5)
    axisZ=min(3.0*(ind.height-1.5)/float(2*(number_of_cavities+1)), axisZ+epsilon)
    axisZ=min(axisZ, (ind.height-cavity[0]))
    axisZ=min(axisZ, cavity[0]-0.5)
    axisZ=max(ind.thickness+1.0, axisZ)
    ind.listCavities[index][4]=axisZ

def mutation_rotation(ind):

    length=len(ind.listCavities)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    index=random.randint(1,length-2)

    axisX = ind.listCavities[index][2]
    axisY = ind.listCavities[index][3]
    theta = random.uniform(-math.pi/4.0,math.pi/4.0)

    axisX, axisY = math.cos(theta)*axisX-math.sin(theta)*axisY, math.sin(theta)*axisX+math.cos(theta)*axisY

    ind.listCavities[index][2]=axisX
    ind.listCavities[index][3]=axisY




def mutation_type(ind):

    length=len(ind.listCavities)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    index=random.randint(1,length-2)

    type=ind.listCavities[index][1]

    if type=="ellipsoid":
        ind.listCavities[index][1]="frisbee"
    else:
        ind.listCavities[index][1]="ellipsoid"


def mutation(ind):
    if mutationAxisX=="ON":
        if random.choice([True, False]):
            mutation_axisX(ind)

    if mutationAxisY=="ON":

        if random.choice([True, False]):
            mutation_axisY(ind)

    if mutationAxisZ=="ON":

        if random.choice([True, False]):
            mutation_axisZ(ind)

    if mutationType=="ON":

        if random.choice([True, False]):
            mutation_type(ind)

    if mutationRotation=="ON":

        if random.choice([True, False]):
            mutation_rotation(ind)

def mutationFunc(pop, params):

    number_of_mutated_ind=params["nbIndMutated"]
    number_of_mutation_per_ind=params["nbMutationsPerInd"]

    if not (isinstance(number_of_mutated_ind, int) and isinstance(number_of_mutation_per_ind, int)):
        raise ValueError, "nbIndMutated and nbMutationsPerInd are not given or are not int"

    length_temp=len(pop)

    for i in range(number_of_mutated_ind):
        j=random.randint(0,length_temp-1)
        ind=newIndividualFrom(pop[j])

        for k in range(number_of_mutation_per_ind):
            mutation(ind)

        pop.append(ind)

    return pop

####################################################################################################
###CROSSING
###
def crossing_ind(individual1, individual2):

    length1=len(individual1.listCavities)

    if length1!=len(individual2.listCavities):

        raise ValueError, "they have not the same number of cavities"

    ind1=newIndividualFrom(individual1)
    ind2=newIndividualFrom(individual2)
    index=random.randint(0,length1-1)

    temp1 = ind1.listCavities[index]
    temp2 = ind2.listCavities[index]

    ind1.listCavities[index]=temp2
    ind2.listCavities[index]=temp1

    return (ind1, ind2)


def crossFunc(pop, params):

    number_of_crossing=params["crossTx"]

    if not isinstance(number_of_crossing, int):
        raise ValueError, "crossTx is not given or is not an int"
    print("crossFunc "+str(pop))

    length_temp=len(pop)

    newpop = algorithm.Population()
    for i in range(number_of_crossing):
        j=random.randint(0, length_temp-1)
        k=random.randint(0, length_temp-1)
        (a,b)=crossing_ind(pop[j], pop[k])
        newpop.append(a)
        newpop.append(b)
    return newpop

####################################################################################################
### Generate
###
def generateIndividual(aType):
        individual=AccordionIndividual()

        for i in range(1,number_of_cavities+1):
            height=0.5+i*(individual.height-1.0)/float(number_of_cavities+1)

            if aType=="ellipsoid":
                if generate_random=="ON":
                    axisX=max(individual.thickness+1.0, random.uniform((4.0/3.0)*individual.radius,(10.0/3.0)*individual.radius))
                    axisY=max(individual.thickness+1.0, random.uniform(2.0*individual.thickness,(10.0/3.0)*individual.radius))
                    axisZ=max(individual.thickness+1.0, (individual.height-1.0)/float((2*(number_of_cavities+1))))

                else:
                    axisX=max(individual.height/2.5 , max(individual.thickness+1.0, (7.0/3.0)*individual.radius))
                    axisY=max(individual.height/2.5 , max(individual.thickness+1.0, (7.0/3.0)*individual.radius))
                    axisZ=max(individual.thickness+1.0, (individual.height-1.0)/float((2*(number_of_cavities+1))))

            else:
                if generate_random=="ON":
                    axisX=max(individual.thickness+1.0, random.uniform((4.0/3.0)*individual.radius,(10.0/3.0)*individual.radius))
                    axisY=max(individual.thickness+1.0, random.uniform(2.0*individual.thickness,(10.0/3.0)*individual.radius))
                    axisZ=max(individual.thickness+1.0, 1.5*(individual.height-1.0)/float((2*(number_of_cavities+1))))
                else:
                    axisX=max(individual.height/2.5 , max(individual.thickness+1.0, (7.0/3.0)*individual.radius))
                    axisY=max(individual.height/2.5 , max(individual.thickness+1.0, (7.0/3.0)*individual.radius))
                    axisZ=max(individual.thickness+1.0, 1.5*(individual.height-1.0)/float((2*(number_of_cavities+1))))

            cavity=[height,aType,axisX,axisY,axisZ]
            accordion.addCavity(individual,cavity)

        cavityBottom=[0.5,"cylinder",5.0*individual.radius/3.0,5.0*individual.radius/3.0,0.75]
        cavityTop=[individual.height,"cylinder",5.0*individual.radius/3.0,5.0*individual.radius/3.0,0.75]

        accordion.addCavity(individual,cavityBottom)
        accordion.addCavity(individual,cavityTop)

        individual.level=None

        return individual


def generateFunc(numgen, params):
    nbInd = params["nbInd"]
    pop = algorithm.Population()

    for i in range(nbInd):
        ltype=random.choice(type)
        pop.append(generateIndividual(ltype))

    return pop

####################################################################################################
### Eval
###
def evaluationFunc(pop):
    global thickness
    print("Evaluation Function "+str(len(pop)))
    basedir=os.path.dirname(__file__)
    bestscore = -float(sys.maxint)

    filename=[]
    for ind in pop:
        shape, shapeMinus = getShapeFromInd(ind)

        ### return (shape,shapeMinus)
        ###
        fend =  "def evalField(x,y,z): \n\treturn shape.eval(primitives.Point(x,y,z))"
        f1 = shapewriter.toPythonString(shape) + fend
        f2 = shapewriter.toPythonString(shapeMinus) + fend

        filename.append((f1,f2, ind))

    #################### EXAMPLE USING THE SEQUENTIAL LAUNCHER #################################
    ### List of filename that contains the simulation to run
    scenefiles = ["scene.pyscn","controller.py", "shape.py", "shapeinv.py"]
    filesandtemplates = []
    for scenefile in scenefiles:
        filesandtemplates.append( (open(basedir+"/"+scenefile).read(), scenefile) )

    for f1, f2, ind in filename:
        runs = []
        for f1,f2,ind in filename:
            runs.append( {"GENERATION": str(pop.id),
                          "INDIVIDUAL": str(ind.id),
                          "SHAPECONTENT": f1, "SHAPEINVCONTENT": f2, "nbIterations":10,
                          "LIBRARYPATH" : os.path.dirname(geneticalgo.__file__),
                          "THICKNESS" : thickness
                          } )

    results = launcher.startSofa(runs, filesandtemplates, launcher=launcher.SerialLauncher())

    for res in results:

               print("Results: ")
               print("    directory: "+res["directory"])
               print("        scene: "+res["scene"])
               print("      logfile: "+res["logfile"])
               print("     duration: "+str(res["duration"])+" sec")

    ### Associate the results to the individuals.
    for i in range(len(filename)):
        f1, f2, ind = filename[i]
        ind.results = results[i]
        data = getJSONFragmentFrom( ind.results["logfile"] )


        if data == None:
            print "SOFA CRASHED DOWN!!  TRY TO LAUNCH MANUALLY scene.pyscn"
            ind.level = - float(sys.maxint)
        else:

            Z0 = data["Z0"]
            Zmax = data["Zmax"]
            V0 = data["V0"]

            level =(Zmax-Z0)


            if level > 20.0:
                print "a scene is ill defined, excessive result"
                ind.level = - float(sys.maxint)
            elif level < 0.0:
                print "becareful, negative result"
                ind.level = level
            else:
                ind.level = level



    ## Wait here the results.
    for f1,f2,ind in filename:
        if ind.level > bestscore:
            bestscore = ind.level

    return bestscore


####################################################################################################
### Select
###
def selectionFunc(pop, params):

    nbInd=params["nbInd"]

    if not isinstance(nbInd, int):
        raise ValueError, "nbInd is not given or not an int"
    pop.sortAccordingToLevel()
    pop.deleteFrom(params["nbInd"])
#    pop.sortAccordingToIndex()

    return pop
