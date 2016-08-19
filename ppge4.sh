#!/bin/bash
# set -x 
function cleanup {
    if [[ -d $TMPDIR ]] ; then 
	rm -fr $TMPDIR
    fi
    exit
}
trap cleanup  SIGHUP SIGINT SIGTERM    

#########################################
###  bash version of epidewarp4.ucsd  ###
###     may csh rest in peace...      ### 
#########################################
#
# ppge4
#
# Pre-process GE DICOM files and then calls epidwarp4.ucsd to perform unwarping
# 
#
# Requires: FSL tools from fMRIB 
#           AFNI tools from NIH
#           dicomrx from cfmriweb.ucsd.edu/fmap/dicomrx
#
# Calls:  epidewarp4.ucsd and dicomrx
#
# What this script does:
#  1) Creates NIFTI from DICOM files. 
#  2) Determines echo times from DICOM headers
#  3) Converts field map DICOM to NIFTI
#  4) Generates complex AVW volumes and magnitudes for masking
#  5) Determine EPI dwell time from DICOM header
#  6) Converts EPI DICOM files to NIFTI
#  7) Calls epidewarp4.ucsd to unwarp the files
#
# New capabilities:
#  1. multiple coils
#
# Version History
#  1.0 040727 TTL initial version
#      040926 TTL adding comments, checking into CVS. 
#      040927 TTL adding option not to merge  EPI AVW volume on input
#  2.0 040928 TTL adding (1) AFNI BRIK output option
#                        (2) option to disable EPI AVW output volume creation
#  3.0 040929 TTL motion correction options added.
#  3.1 041001 TTL added nocleanup option
#  3.2 041015 TTL added phswap option
#  3.3 050413 GTB added fmap option
#  3.4 050505 GTB added briktype option
#  3.5 061605 GTB added NIFTI support
#  3.6 092705 KL  added Rx info to the output hdr
#  3.7 011006 GTB added nogzip option, and a flexible postfix specification 
#  3.8 081506 KL  corrected mismatch of nifti header information by using 
#                 avwcpgeom in multiple places. Ohterwise dewarping won't be correct. 
#  3.9 040507 KL  added support for multiple coils. added option for ncoil and nphase. 
# 3.10 042607 KL  added matrix size check. In case of size conflicts, fm is resampled to match
#                 the size of the data. 
# 3.11 042707 KL  removed -TR option. now extract TR from dicom header (0018,0080). 
# 3.12 080207 CBK  Changed all calls from dicom2 (which is not compiled for the mac) to the afni 
#               command line program dicom_hdr (which is).  
# 3.13 090707 KL  updated all avw... commands to fsl... to work with FSL4.0.
# 3.14 031508 KL  Added ASSET detection and handling. 
# 3.15 06/17/2010 KL added -n option to allow users to supply their own epivol in .nii.gz format.
# 4.00 08/31/2012 MD ported to bash
# Send Comments/Questions to kunlu@ucsd.edu or ttliu@ucsd.edu
#

# ==================================================
Run() {
    cmd=$1
    echo $cmd | tee -a $LF
    $cmd
    [[ $? -ne 0 ]] && { echo -e "\n \
**************************************************\n \
Command $cmd failed.\n PWD = $PWD\n Exiting...\n \
**************************************************\n" ; exit 1; } >&2
}
# --------------------------------------------------
ParseArgs() {
    set -- $*

    echo 1 $1 2 $2
    while [[ $# -gt 0 ]] ; do
	case $1 in 
	    -d1)
		d1=$2
		shift 2;;
	    -d2)
		d2=$2
		shift 2;;
	    -i)
		epidir=$2
		shift 2;;
	    -n)
		epinifti=$2
		usrnifti=1
		shift 2;;
	    -o)
		outstem=$2
		shift 2;;
	    -TR)
		echo "ERROR: flag $flag is obsolete!"
		exit 1;;
	    -tmpdir)
		TMPDIR=$2
		shift 2;;
	    -tediff)
		tediff=$2
		shift 2;;
	    -nomask)
		NOMASK=1  #"--outmask"
		shift;;
	    -nobrik)
		dobrik=0
		shift;;
	    -nomoco)
		NOMOCOIFSET=$1
		shift;;
	    -nocleanup)
		docleanup=0
		shift;;
	    -fmap)
		FMAPIFSET="--fmap=$2"
		shift 2;;
	    -refnum)
		REFNUMIFSET="--refnum=$2"
		shift 2;;
	    -unwarpdir)
		UNWARPDIRIFSET="--unwarpdir $2"
		shift 2;;
            -ncoil)
		ncl=$2
		shift 2;;
	    -nphase)
		nphase=$2
		shift 2;;
	    -betf)
		betf=$2
		shift 2;;
	    *)
		echo "ERROR: Flag $1 not recognized. "
		exit 1
	esac
    done

    echo d1 $d1 d2 $d2 
}
##################################################
UsageExit() {
    echo "Name"
    echo "     ppge4  - preprocesses GE DICOM FILES for unwarping and then"
    echo "              calls epidewarp4.ucsd"
    echo ""
    echo "System requirements"
    echo "     AFNI - AFNI_2006 and up, 32bit version"
    echo "     FSL  - FSL3.2 and up, FSLOUTPUTTYPE = NIFTI_GZ (or NIFTI if use -nogzip option)"
    echo ""
    echo "Synopsis"
    echo "     ppge4 -d1 <TE1 data dir> -d2 <TE2 data dir> -i<EPI dir> -o <outstem> [<options>]"
    echo ""
    echo "Required Arguments"
    echo "     -d1 <TE1 data dir>"
    echo "     -d2 <TE2 data dir>"
    echo "     -i <input EPI data dir>"
    echo "     -o <outstem>"
    echo "     Environment variable FSLOUTPUTTYPE need to be set to NIFTI_GZ or NIFTI "
    echo "     (use -nogzip option if NIFTI)"
    echo ""
    echo "Optional Arguments"
    echo "     -n <epi.nii.gz>:    users supplied epivol data (must be nii.gz format) "
    echo "     -tmpdir <tmpdir>:   temporary file directory; default: <current directory>/tmp "
    echo "     -tediff <tediff>:   TE difference [us]; default: auto"
    echo "     -unwarpdir <dir>:   unwarping direction = x / y / z / x- / y- / z-, default = y"
    echo "     -avwvol         :   enables AVW output volume creation"
    echo "     -fmap           :   enables field map output volume creation"
    echo "     -nomask         :   disables brain masking of output EPI volume" 
    echo "     -nomoco         :   disables motion correction of EPI prior to unwarping"
    echo "     -refnum <refnum>:   Reference image number [0,nframes-1] for motion correction, default is middle image in a series"
    echo "     -nocleanup      :   disables removal of temporary files"
    echo "     -ncoil          :   number of receivers used to acquire fieldmap (needed only when using coils other than BODY and 8HRBRAIN)" 
    echo "     -nphase         :   number of phase split in phase unwrapping (default =8)" 
    echo "     -betf         :   fractional intensity threshold for BET (0->1) (default =0.5)" 
    echo "     -nogzip         :   don't gzip NIFTI files" 
    echo "" 
    echo "Outputs"
    echo "     <outstem> - unwarped volume filename stem" 
    echo ""
    echo "Version"
    echo "     "$VERSION
    echo ""
    echo "Credits"
    echo "     FSL library" 
    echo ""
    echo "Reporting Bugs"
    echo "     Report bugs to kunlu@ucsd.edu"
    echo ""

    exit 1
}
# ==================================================
ExitIfHaveNiftis() {
    if ls *.nii.gz &>/dev/null ; then 
	echo "Some nifti files are already present in $PWD. Remove them first." >&2
	exit 1
    fi
}
# ==================================================
TMPDIR=tmp.$$
merge_in=1
domask=1
dobrik=0
docleanup=1
ncl=-1
nphase=8
betf=0.5
cleanup=1
dopmask=1
dofmask=1
dounwrap=1
tediff=auto
AFNIDIR=${AFNIDIR:-/usr/local/afni}


[[ $# -lt 1 ]] && UsageExit

ParseArgs $*


## Get/Create tmp directory ##
mkdir -p $TMPDIR


curdir=`pwd`
echo -e "\nGrabbing TE1 and TE2 Data\n"

# determine # field map files
nfiles=`ls -1 $curdir/$d1/i* | wc -l`
echo nfiles $nfiles

#determine number of recerivers used
#Changed from dicom2 to dicom_hdr (afni) for the rest of us.
# -cbk
if [[ x$ncl = x-1 ]] ; then
    f1=`ls $d1/i*.1`
    coil=`$AFNIDIR/dicom_hdr $f1 | grep "0018 1250" | awk -F// '{print $3}'`

    case $coil in
	8HRBRAIN)
	    ncl=8;;
	BODY)
	    ncl=1;;
	*)
	    echo "ERROR: COIL $coil is not recognized. "
	    echo "ERROR: Please specify -ncoil "
	    exit 4
    esac

    echo "INFO: $coil coil detected - use $ncl receiver(s)" 
else
    echo "INFO: Using $ncl receivers as user specified" 
fi

[[ $ncl -ne 1 ]] && { echo "This script doesn't handle more than one receiver. Use original csh script or adjust this accordingly. Exiting... " ; exit 1; } >&2

# assuming mag,ph,re,im
mag_te1=$curdir/$TMPDIR/mag_te1
re_te1=$curdir/$TMPDIR/re_te1 #_cl$cnt
im_te1=$curdir/$TMPDIR/im_te1
re_te2=$curdir/$TMPDIR/re_te2
im_te2=$curdir/$TMPDIR/im_te2

# -------------------- first fieldmap --------------------
## cd $d1
echo `pwd`
ExitIfHaveNiftis

    #convert dicoms:
##dcm2niix -o ./ $(ls i* | head -1)
dcm2niix -z y $d1
F=$(ls -tr $d1/*.nii.gz | tail -1)     # grab the last nifti file
Run "fslorient -setsformcode 0 $F"
Run "fslroi $F $re_te1 2 1 "  #re
Run "fslroi $F $im_te1 3 1 "  #im

    # make combined magnitude images of te1 for masking purpose later.
Run "fslroi $F $mag_te1 0 1 "  #mag
mag_te1=$TMPDIR/mag_te1 #get rid of absolute path for mag; since this is passed on to epidewarp

rm -f $F

# -------------------- 2nd fieldmap --------------------
#cd $d2
echo `pwd`
ExitIfHaveNiftis

# convert dcm --> nii.gz
## dcm2niix -o ./ $(ls i* | head -1)
dcm2niix -z y $d2
F=$(ls -tr $d2/*.nii.gz | tail -1)
Run "fslorient -setsformcode 0 $F"
Run "fslroi $F $re_te2 2 1 "
Run "fslroi $F $im_te2 3 1 "

rm -f $F

cd $curdir

# determine echo times
# Changed dicom2 to dicom_hdr
# -cbk
if [[ _$tediff = _auto ]] ; then
    f1=`ls $d1/i*.1`
    te1=`$AFNIDIR/dicom_hdr $f1 | grep "0018 0081" | awk -F// '{print $3}'`
    
    f1=`ls $d2/i*.1`
    te2=`$AFNIDIR/dicom_hdr $f1 | grep "0018 0081" | awk -F// '{print $3}'`
    
    tediff=$(echo "1000.0*($te2 - $te1)" | bc -l) 
    echo "INFO: tediff = $tediff usec" 
fi

# weighted mean method
echo "Calculate complex volume using conjugate product"
Run "fslmaths $re_te2 -mul $re_te1 $TMPDIR/grot1 -odt float"
Run "fslmaths $im_te2 -mul $im_te1 $TMPDIR/grot2 -odt float"
Run "fslmaths  $TMPDIR/grot1 -add $TMPDIR/grot2 -div 1000 $TMPDIR/re_dph -odt float"

Run "fslmaths $re_te2 -mul -1 -mul $im_te1 $TMPDIR/grot1 -odt float"
Run "fslmaths $im_te2 -mul $re_te1 $TMPDIR/grot2 -odt float"
Run "fslmaths  $TMPDIR/grot1 -add $TMPDIR/grot2 -div 1000 $TMPDIR/im_dph -odt float"

rm -f $TMPDIR/grot*

# summation
#echo "Combine $ncl coil(s) ... "

#fslmaths $TMPDIR/re_dph_cl1 -mul 0 $TMPDIR/re_dph
#fslmaths $TMPDIR/im_dph_cl1 -mul 0 $TMPDIR/im_dph

#for cnt in $(seq 1 $ncl) ; do
#    fslmaths $TMPDIR/re_dph -add $TMPDIR/re_dph_cl$cnt $TMPDIR/re_dph
#    fslmaths $TMPDIR/im_dph -add $TMPDIR/im_dph_cl$cnt $TMPDIR/im_dph
#done

# generate complex AVW volumes of the conjugate product
cp_dph=$TMPDIR/cp_dph
Run "fslcomplex -complex  $TMPDIR/re_dph $TMPDIR/im_dph $cp_dph"

nslices=`fslinfo $mag_te1 | awk '$1 == "dim3"{print $2}'`
gp=$mag_te1

### ----------- get some parameters from dicom header ----------
#determine EPI dwell time
#Changed dicom2 to dicom_hdr
f1=`ls $epidir/*.1 `
dwell=`$AFNIDIR/dicom_hdr $f1 | grep -i "0043 102c" | awk -F// '{print $3}'`
echo -e "\nINFO: EPI dwell time = $dwell usec"

asset=`$AFNIDIR/dicom_hdr $f1 | grep -i "0043 1083" | awk -F// '{print $3}' | cut -d'\' -f1`
if [[ -z $asset ]] ; then
    echo "     no asset" 
else
    echo "      asset factor = $asset"
    dwell=`echo "scale = 2; $dwell*$asset " | bc -l`
    echo "      Adjusted dwell time = $dwell usec"
fi

dratio=`echo "scale = 10; $dwell/$tediff " | bc -l `
echo -e "INFO: dratio = $dratio"

#Changed dicom2 to dicom_hdr
TR=`$AFNIDIR/dicom_hdr $f1 | grep "0018 0080" | awk -F// '{print $3}'`
echo -e "INFO: EPI TR = $TR msec\n"
# --------------------------------------------------

if [[ ! -z $usrnifti ]] ;  then

    echo "Use user-supplied NIFTI for EPI VOL $usrnifti"   
    imcp $epinifti $curdir/$TMPDIR/epivol.nii.gz

else
    
    cd $curdir/$epidir
    echo "Creating NIFTI for EPI DICOM"
    ExitIfHaveNiftis

    dcm2niix i*.1 # the actual conversion 
    [[ $(ls *.nii.gz | wc -l ) -gt 1 ]] && { echo "Ambiguity, have more than 1 niftis generated in $PWD. Exiting"; exit 1; } >&2
    mv *.nii.gz $curdir/$TMPDIR/epivol.nii.gz
    cd $curdir

fi


# Check matrix size and resample fmap data if needed (KL):
matfm=`fslinfo $mag_te1 | awk '$1~/^dim[123]/{print $2}'`
matepi=`fslinfo $TMPDIR/epivol | awk '$1~/^dim[123]/{print $2}'`
if [[ $matfm != $matepi ]] ; then 

    echo "**************************** WARNING ****************************************"
    echo "WARNING: fieldmap data ($matfm) does not match the functional data ($matepi) !"
    echo "Resampling the fieldmaps to ($matepi) ... "
    echo "*****************************************************************************"
    
    Run "3dresample -master $TMPDIR/epivol -prefix $TMPDIR/mag_te1_rs -inset $TMPDIR/$mag_te1"
    
    gp=$TMPDIR/mag_te1_rs

    Run "3dresample -master  $gp -prefix $TMPDIR/re_dph_rs -inset $TMPDIR/re_dph"    
    Run "3dresample -master  $gp -prefix $TMPDIR/im_dph_rs -inset $TMPDIR/im_dph"
    
    cp_dph=$TMPDIR/cp_dph_rs
    fslcomplex -complex  $TMPDIR/re_dph_rs $TMPDIR/im_dph_rs $cp_dph

else

    echo "check ... fieldmap data ($matfm) matches the functional data ... OK"

fi

echo ">>> ppge done <<<"
echo "**************************************************"
echo "calling epidewarp4.ucsd.sh"


# call epidewarp now
[[ x$NOMASK = x1 ]] || OUTMASKIFSET="--outmask"
Run "epidewarp4.ucsd.sh --mag $gp --dph $cp_dph --epi $TMPDIR/epivol \
               --tediff $tediff  --esp $dwell  --vsm $TMPDIR/vsm --exfdw $TMPDIR/ex \
               --epidw $outstem --cleanup --nphase $nphase --betf $betf \
               --tmpdir $TMPDIR $OUTMASKIFSET $NOMOCOIFSET $REFNUMIFSET $UNWARPDIRIFSET $FMAPIFSET"

# --postfix .nii.gz taken out. nii.gz default
echo "Done epidewarp from ppge4"
