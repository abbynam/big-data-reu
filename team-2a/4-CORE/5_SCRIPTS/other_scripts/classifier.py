import os, sys, array
import rpy2
import rpy2.robjects as robjects
import cProfile

rpy2.rinterface.initr()
robjects.r.library('MASS')


#DATA_FILE = "/home/dsmackin/apps/dca_R/data/444_10K.csv"
#DATA_FILE = "/home/dsmackin/apps/dca_R/data/gamma_data_small.csv"
#DATA_FILE = "/home/dsmackin/apps/dca_R/data/beam_LT_1_GT_3_weka.csv"
DATA_FILE = "/y_drive/CCData/Co60/Co60_df_small.csv"

MIN_EXPLANATORY_VARS = 2
MAX_EXPLANATORY_VARS = 5

RESPONSE_VARIABLE = "class"
RESPONSE_VARIABLE = "I(dca < 20)"
# SEPARATOR="\\t"
SEPARATOR=","


def getSubsets( vals ):

  subsets = [ [vals[0]] ]
  newSubsets = []
  if len(vals) > 1:
      newSubsets = getSubsets( vals[1:] )

  for set in newSubsets:
      subsets.append( set )
      subsets.append( set + [ vals[0] ] )

  subsets.sort()
  return subsets


def loadDataInR( dataFile, selection ):

   R_command = '''
dataAll<-read.csv("%s", header=TRUE, sep='%s')
dataAll <- dataAll[dataAll$E > 1.0 & dataAll$E < 1.5, ]
fieldlist = c(names(dataAll))
attach(dataAll)

# rbind(dataAll)

   ''' %( dataFile, SEPARATOR)


   print R_command, "\n"
   print robjects.r( R_command  )
   print ""


def getVarlist( varArray ):
   subsets = getSubsets(varArray)
   varlists = []
   for set in subsets:
      varString = ""
      delim = ""

      if len(set) < MIN_EXPLANATORY_VARS or len(set) > MAX_EXPLANATORY_VARS:
         continue
      for i in range( len(set) ):
         varString = "%s%s%s" % ( varString, delim, set[i])
         delim = " + "
         #delim = "*"
      varlists.append( varString )
   return varlists


def testExplanatoryVars( reponseVar, varlists, preselectedVariables ):
    glmResults = []
    glmBestResult = [ 99999, "", [], [] ]
    for fields in varlists:
        #glm_muTrigs = robjects.r.glm("formula = muTrig ~ I(muEta^8)+I(muEta^6)+I(muEta^4)+muInHole")
        glmString = "dataAll$%s ~ %s + %s" % ( reponseVar, fields, preselectedVariables )
        glmString = "%s ~ %s + %s" % ( reponseVar, fields, preselectedVariables )


        #glmString = "(dataAll$distance_closest_approach < 5.0) ~ %s + %s" % ( fields, preselectedVariables )
        #print "\n\n%s" % glmString


        #glm_muTrigs = robjects.r.glm(formula=glmString, family=robjects.r.quasibinomial())
	#glm_muTrigs = robjects.r.glm(formula="dataAll$class ~ dataAll$E1", family=robjects.r.binomial("logit"))

        print glmString
        glm_muTrigs = robjects.r.glm(formula=glmString, family=robjects.r.binomial("logit"))

        glm_summary = "%s" % ( robjects.r.summary(glm_muTrigs) )
        glm_summary = glm_summary.split("\n")

        for i in range(len(glm_summary)):
            string = "%s" % glm_summary[i]
            #print string
            if string.find('AIC') != -1:
#                print string
               aic = float(string.split(" ")[-1].strip())
               figOfMerit = aic
#            if string.find('Residual deviance:') != -1:
#                try:
#                    parts = string.split(" ")
#                    floatArray = []
#                    for part in parts:
#                        try:
#                            floatArray.append( float(part.strip()))
#                        except ValueError:
#                            1 # do nothing
#                    #figOfMerit = floatArray[0]/floatArray[1]
#                    figOfMerit = floatArray[0]
#                except ValueError:
#                    print "WARNING: NO DEVIANCE:\n%s\n" % string
#                    print "aic=%s" % (string.split(" ")[2].strip())
#                    figOfMerit= 12345
                #print aic


        coeffs = robjects.r.coefficients(glm_muTrigs)
        coeffsArray = array.array('d',coeffs)
        covmat = robjects.r.vcov(glm_muTrigs)
        covmatArray = array.array('d',covmat)
        record = [figOfMerit, fields, coeffsArray, covmatArray]
        glmResults.append( record )
        if figOfMerit < glmBestResult[0]:
            glmBestResult = record

    return glmResults, glmBestResult




def testVars():

   loadDataInR( DATA_FILE , " 1==1" )
   varlists = getVarlist([ 'E', 'E1','E2', \
                # 'x1', 'x2', 'y1', 'y2', 'z1', 'z2', \
                'dE', \
                'I(E1-E2)', \
                #'dE_scaled', \
			    'theta1', \
                #'delta_theta', 'cos_delta_theta', 'delta_cos_theta',\
                'sin_delta_theta',\
                # 'cos_theta_known', 'theta_known',   \
                  'alpha', 'phi',\
                           ])

   preselectedVariables = '1'

   glmResults_mu, best_mu =  testExplanatoryVars(RESPONSE_VARIABLE, varlists, preselectedVariables)
   glmResults_mu.sort()
   glmResults_mu.reverse()
   for rec in glmResults_mu:
     print "%.5f(%d)  %s" % ( rec[0], len(rec[1].split("+")), rec[1] )

   print "\n\n", best_mu[0], best_mu[1]
   numVars = len( (best_mu[2]) )
   print "--------------\n[",
   for val in best_mu[2]:
     print val, ", ",
   print "]\n\n--------------\n[",
   for i in range(len(best_mu[3])):
     print best_mu[3][i], ",",
     if (i + 1) % numVars == 0:
         print " \n\t\t",
   print "]\n\n",


#----------------------------------------------------------------------
#------------------------------- MAIN ---------------------------------
#----------------------------------------------------------------------

# testVars()
cProfile.run("testVars()", sort='time')


