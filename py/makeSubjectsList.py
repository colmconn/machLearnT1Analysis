#!/bin/env python

import os.path, re, string, glob
import pandas as pd
import numpy as np

root_dir="/data"
data_dir=os.path.join(root_dir, "sanDiego", "rsfcGraphAnalysis", "data")
admin_data_dir=os.path.join(data_dir, "admin")
# config_data_dir=os.path.join(data_dir, "config")

freeSurfer_root=os.path.join(root_dir, "sanDiego", "freesurferAnalysis")
freeSurfer_root_data_dir=os.path.join(freeSurfer_root, "data")

def makeListOfSubjectFolders():
    contents=glob.glob(os.path.join(data_dir, "*_[ABCDE]"))
    return sorted(map(os.path.basename, contents))

def pickSubjectsByTimepoint(inSubjectList, inTimepoint):
    rr=re.compile("^[0-9]{1,4}_" + inTimepoint + "$")
    return filter(rr.match, inSubjectList)

def removeTimepointFromSubjectId(inIdsList):
    return map(lambda xx: re.sub("_[A:Z]", "", xx), inIdsList)

def readDemographicsFile(inDemographicsFilename):
    print("*** Reading %s" % inDemographicsFilename)
    demographics=pd.read_csv(inDemographicsFilename, comment="#",
                             na_values = ["NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""])
    print "*** Read demographics data for %s unique subjects" %  len(np.unique(demographics['ID']))

    return(demographics)

def makeListOfBrainMaskFiles(inSubjectList):
    return map(lambda xx: os.path.join(freeSurfer_root_data_dir, xx, "mri", "brainmask.mgz"), inSubjectList)

def restrictToExistingBrainMasks(inDf):
    return inDf[inDf.apply(lambda xx: os.path.isfile(xx["filename"]), axis=1)]

def splitSubjectIdIntoIdAndTimepoint(inDf):
    newDf = inDf["subjectId"].apply(lambda xx: pd.Series(xx.split('_')))
    newDf=pd.concat([newDf, inDf], axis=1)
    newDf.rename(columns={0:"ID", 1:"timepoint"}, inplace=True)
    return (newDf)

#def fixSubjectId(inDf, inId="ID"):

def fixDates(inData):
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inData['DOB'] = pd.to_datetime(inData['DOB'], coerce=True, unit="d")
    inData['MRI'] = pd.to_datetime(inData['MRI'], coerce=True, unit="d")

    return(inData)


def computeAge(inData):

    # use this convoluted apply methodology to get rid of the units in the timeseriesdelta
    # see http://stackoverflow.com/questions/17414130/pandas-datetime-calculate-number-of-weeks-between-dates-in-two-columns
    age_in_days=(inData['MRI'] - inData['DOB']).apply(lambda x: x/np.timedelta64(1,'D'))
    age_in_weeks=age_in_days/7
    inData['age_in_years']=age_in_weeks/52

    return(inData)

def convertMgzFile(inMgzFile, inOutputFile, inOutType="nii"):

    command="mri_convert --in_type mgz --out_type %s %s %s" % (inOutType, inOutputFile, outputFile)
    print command

def makeNiftiFilenames(inDataFrame, infix, filenameColumn):
    inDataFrame[filenameColumn]=inDataFrame.apply(lambda g, s, i: os.path.join(data_dir, "%s.%s.%s.nii.gz" % (g, s, i), inDataFrame["Grp"], inDataFrame["subjectId"], infix), axis=0)
    return inDataFrame

# def convertMgzFiles(inDataFrame):

#     inDataFrame[
    
def main():

    demographicsFilename=os.path.join(admin_data_dir, "0-data_entry_current_2014.csv")

    demographics=readDemographicsFile(demographicsFilename)

    #print demographics.iloc[0:5,0:5]

    ## subjectsAtTimepointA=removeTimepointFromSubjectId(pickSubjectsByTimepoint(makeListOfSubjectFolders(), "A"))
    subjectsAtTimepointA=pickSubjectsByTimepoint(makeListOfSubjectFolders(), "A")
    print "*** Found folders for %s subjects at timepoint A" % len(subjectsAtTimepointA)

    subjectIdAndBrainMaskFileList=pd.DataFrame(zip(subjectsAtTimepointA, makeListOfBrainMaskFiles(subjectsAtTimepointA)), columns=["subjectId", "filename"])
    subjectIdAndBrainMaskFileList=splitSubjectIdIntoIdAndTimepoint(restrictToExistingBrainMasks(subjectIdAndBrainMaskFileList))

    #print subjectIdAndBrainMaskFileList["ID"]
    #print subjectIdAndBrainMaskFileList    
    #print subjectIdAndBrainMaskFileList[subjectIdAndBrainMaskFileList["ID"] == "169/300"]
    #print listOfExistingBrainMasks(subjectIdAndBrainMaskFileList)

    mriFilesAndDemographics=pd.merge(demographics, subjectIdAndBrainMaskFileList, left_on="ID", right_on="ID", sort=False)
    #print mriFilesAndDemographics

   
    mriFilesAndDemographics = fixDates(mriFilesAndDemographics)
    print mriFilesAndDemographics
    mriFilesAndDemographics = computeAge(mriFilesAndDemographics)
    print mriFilesAndDemographics

    ## remove any rows for which the age in years couldn't be calculated
    ## print mriFilesAndDemographics[ mriFilesAndDemographics['age_in_years'].isnull()]
    mriFilesAndDemographics=mriFilesAndDemographics[mriFilesAndDemographics.age_in_years.notnull()]
    print mriFilesAndDemographics.age_in_years.isnull().any()

    mriFilesAndDemographics = makeNiftiFilenames(mriFilesAndDemographics, "anat_struc_brain", "skullstrippedFilename")
    print mriFilesAndDemographics.columns
    #print mriFilesAndDemographics[["subjectId", "Grp", "skullstrippedFilename"]]

    


####################################################################################################
            
if __name__ == "__main__":
    main()
