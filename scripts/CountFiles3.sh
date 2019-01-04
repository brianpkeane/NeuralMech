#!/bin/bash
# Written by brian keane, 12/11/18
# Identify subjects with inappropriate number of files

# Count files for a given experiment
listOfSubjects="sub-B04 sub-B06 sub-B10 sub-B14 sub-B21 sub-B23 sub-B27 sub-B29 sub-C02 sub-C05 sub-C10 sub-C13 sub-C14 sub-C15 sub-C22 sub-C23 sub-C24 sub-C25 sub-C26 sub-C27 sub-C29 sub-C32 sub-C33 sub-S06 sub-S10 sub-S12 sub-S16 sub-S20 sub-S25 sub-S31 sub-S32"
listOfSubjectsOldProtocol="sub-B04 sub-B06 sub-B10 sub-B14 sub-B21 sub-B23 sub-C02 sub-C05 sub-C10 sub-C13 sub-C14 sub-C15 sub-C22 sub-S06 sub-S10 sub-S12 sub-S16 sub-S20 sub-S25" 
listOfScanTypes="_T1w _T2w dir-AP dir-PA viscomp retino contour ebbloc ebbreg rest"
listOfScanNumCorrect="4 4 4 4 16 12 12 4 12 4"
cd $nm_dirname/data
nT1FilesCorrect=4;
nT2FilesCorrect=4;
# AP and PA field maps, 4 each (two from 1st half of session and two from 2nd)
nAPFilesCorrect=4;
nPAFilesCorrect=4
nViscompFilesCorrect=16;
nRetinoFilesCorrect=12;
nContourFilesCorrect=12;
nEbblocFilesCorrect=4
nEbbregFilesCorrect=12;
nRestFilesCorrect=4;

for subjNum in $listOfSubjects
do
    cd $nm_dirname/data/${subjNum}/
    
    nT1FilesThisSubj=$( find ./anat/ -type f -name *${subjNum}*_T1w* -print| wc -l)
    if (("${nT1FilesThisSubj}" != "${nT1FilesCorrect}" )); then
	echo "The number of files is incorrect for T1 anatomicals for ${subjNum}, ${nT1FilesThisSubj}"
    fi

    nT2FilesThisSubj=$( find ./anat/ -type f -name *${subjNum}*_T2w* -print| wc -l)
    if (("${nT2FilesThisSubj}" != "${nT2FilesCorrect}" )); then
	echo "The number of files is incorrect for T2 anatomicals for ${subjNum}, ${nT2FilesThisSubj}"
    fi

    nAPFilesThisSubj=$( find ./fmap/ -type f -name *${subjNum}*dir-AP* -print| wc -l)
    if (("${nAPFilesThisSubj}" != "${nAPFilesCorrect}" )); then
	echo "The number of files is incorrect for AP fieldmaps for ${subjNum}, ${nAPFilesThisSubj}"
    fi

    nPAFilesThisSubj=$( find ./fmap/ -type f -name *${subjNum}*dir-PA* -print| wc -l)
    if (("${nPAFilesThisSubj}" != "${nPAFilesCorrect}" )); then
	echo "The number of files is incorrect for PA fieldmaps for ${subjNum}, ${nPAFilesThisSubj}"
    fi

    #echo "Now counting viscomp files for subject ${subjNum}:" 
    nViscompFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*viscomp* -print| wc -l)
    if (("${nViscompFilesThisSubj}" != "${nViscompFilesCorrect}" )); then
	echo "The number of files is incorrect for viscomp for ${subjNum}, ${nViscompFilesThisSubj}"
    fi

    #echo "Now counting retino files for subject ${subjNum}:" 
    nRetinoFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*retino* -print| wc -l)
    if (("${nRetinoFilesThisSubj}" != "${nRetinoFilesCorrect}" )); then
	echo "The number of files is incorrect for retino for ${subjNum}, ${nRetinoFilesThisSubj}"
    fi

    #echo "Now counting contour files for subject ${subjNum}:"
    nContourFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*contour* -print| wc -l)
    if (("${nContourFilesThisSubj}" != "${nContourFilesCorrect}" )); then
	echo "The number of files is incorrect for contour for ${subjNum}, ${nContourFilesThisSubj}"
    fi
    
    #echo "Now counting ebbloc files for subject ${subjNum}:" 
    nEbblocFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*ebbloc* -print| wc -l)
    if (("${nEbblocFilesThisSubj}" != "${nEbblocFilesCorrect}" )); then
	echo "The number of files is incorrect for ebbloc for ${subjNum}, ${nEbblocFilesThisSubj}"
    fi

    #echo "Now counting ebbreg files for subject ${subjNum}:" 
    nEbbregFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*ebbreg* -print| wc -l)
    if (("${nEbbregFilesThisSubj}" != "${nEbbregFilesCorrect}" )); then
	echo "The number of files is incorrect for ebbreg for ${subjNum}, ${nEbbregFilesThisSubj}"
    fi

    nRestFilesThisSubj=$( find ./func/ -type f -name *${subjNum}*rest* -print| wc -l)
    if (("${nRestFilesThisSubj}" != "${nRestFilesCorrect}" )); then
	echo "The number of files is incorrect for rest for ${subjNum}, ${nRestFilesThisSubj}"
    fi


#if [[ $listOfSubjectsOldProtocol != *"${subjNum}"* ]]; then  
#    echo "Subject ${subjNum} uses the new protocol."
#else
#    echo "Subject ${subjNum} uses the old protocol."
#fi

done
