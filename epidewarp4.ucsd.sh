#!/bin/bash
#########################################
###  bash version of epidewarp4.ucsd  ###
###     may csh rest in peace...      ### 
#########################################

# This script uses FSL prelude and fugue tools to unwarp EPI images.
# Original version is called epidewarp.fsl authored by Doug Greve at MGH for use
# by the fBIRN consortium.
# 
# This version is called epidewarp3.ucsd and is adapted from version 1.1 (2004/7/27) of
# epidewarp.fsl. The intent is to faciliate unwarping from GE DICOM images at UCSD.
# As much as possible, this script will be updated to reflect improvements
# in epidewarp.fsl as these are made available.
#
# For processing at UCSD, this script is intended to be called by the script:  ppge3
#
# Version History
#  1.0  TTL  040728  Added gemode input to support different processing for GE DICOM files at UCSD
#       TTL  040926  Checked into CVS repository  
#       TTL  040927  Modifying to support multiple EPI frames
#  3.0  TTL  040929  Adding motion correction features
#  3.1  TTL  041001  Redefine refnum to go from 0 to nframes-1, in accordance with standard AVW indexing
#  3.2  TTL  041014  Added unwarpdir option
#  3.3  GTB  050413  Added write fieldmap option
#  3.4  GTB  050505  Added -plots flag on mcflirt
#  3.5  GTB  061505  Modified for NIFTI format
#  3.6  KL   081506  Corrected mismatch of nifti header information by using 
#                 avwcpgeom in multiple places. Ohterwise dewarping won't be correct. 
#  3.7  KL   040507  added number of phase splits to use in prelude.
#  3.8  KL   042707  removed option dph2 (dph1 is now the phase difference between echo1 and 2, #                    Therefore dph2 is no longer needed)
#  3.9  KL  090707  updated all avw... commands to fsl... to work with FSL4.0.
#  4.0  KL  06172010: added a step to align brain mask to epi space.
#  5.0  MD  08302012 ported to bash
# Send Comments/Questions to kunlu@ucsd.edu
# --------------------------------------------------
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
ExistsOrExit() {
    [[ $(imtest "$1") = 0 ]] && { echo "File $1 doesn't exist. Exiting..."; exit 2; } >&2
}
# --------------------------------------------------
ParseArgs() {
    set -- $*
    while [[ $# -gt 0 ]] ; do
	case $1 in 
	    --mag) 
		mag=$2
		ExistsOrExit $mag
		shift 2;;
	    --dph)
		dph=$2
		ExistsOrExit $dph
		shift 2;;
	    --epi)
		epi=$2
		ExistsOrExit $epi
		shift 2;;
	    --sepi)
		epi=$2
		epimerged=0
		ExistOrExit ${epi}0001
		shift 2;;
	    --tediff)
		tediff=$2
		shift 2;;
	    --esp)
		esp=$2
		shift 2;;
	    --epidw)
		epidw=$2
		shift 2;;
	    --exfdw)
		exfdw=$2
		shift 2;;
	    --vsm)
		vsm=$2
		shift 2;;
	    --fmap)
		SAVEFMAPIFSET="--savefmap=$2"
		shift 2;;
	    --unwarpdir)
		UNWARPDIRIFSET="--unwarpdir=$2"
		shift 2;;
	    --tmpdir)
		TMPDIR=$2
		cleanup=0
		shift 2;;
	    --nphase)
		nphase=$2;
		shift 2;;
	    --betf)
		betf=$2
		shift 2;;
	    --refnum)
		refnum=$2;
		shift 2;;
	    --nomoco)
		domoco=0
		shift;;
	    --nocleanup)
		cleanup=0
		shift;;
	    --outmask)
		do_outmask=1
		shift;;
	    --cleanup)
		cleanup_forced=1
		shift;;
	    --debug)
		verbose=1
		set -x
		shift;;
	    *)
		echo "ERROR: Option $1 not recognized"
		exit 1;;
	esac
    done
}
# --------------------------------------------------
# --------------------------------------------------
UsageExit() {
    echo ""
    echo "USAGE: epidewarp4.ucsd"
    echo ""
    echo "Inputs"
    echo "  --mag volid     : B0 magnitude volume"
    echo "  --dph volid     : B0 phase difference volume "
    echo "  --epi volid     : epi volume  (note: specify this for merged AVW data)"
    echo "  --tediff tediff : difference in B0 field map TEs"
    echo "  --esp esp       : EPI echo spacing"
    echo "  "
    echo "Outputs"
    echo "  --vsm volid   : voxel shift map (required)"
    echo "  --exfdw volid : dewarped example func volume"
    echo "  --epidw volid : dewarped epi volume"
    echo "  --fmap volid  : fieldmap volume"
    echo "  "
    echo "Options  "
    echo "  --unwarpdir<dir>:   unwarping direction = x / y / z / x- / y- / z-, default = y"
    echo "  --nomoco      : no motion correction"
    echo "  --refnum <refnum>: reference image number for motion correction"
    echo "  --tmpdir dir  : save intermediate results here"
    echo "  --nocleanup   : do not delete tmpdir"
    echo "  --cleanup     : force deletion of tmpdir"
    echo "  --nphase      : number of phase splits to use in prelude"
    echo "  --betf        : fractional intensity threshold (0->1)for BET;"
    echo "  --debug"
    echo "  --help"
    echo ""
    echo "Version"
    echo "     "$VERSION

    awk '/^BEGINHELP/{flag=1; next}/^__EOHELP/{exit(0)}flag' $0
#  sed -n '/^BEGINHELP/,/__EOHELP/p' $0

    exit 1;
}
# ==================================================
[[ $# -lt 1 ]] && UsageExit

ParseArgs $*

LF=${0%%.*}.log

# Some temp files
exf=$TMPDIR/exf
brain=$TMPDIR/brain
head=$TMPDIR/head

# Extract the middle time point for the example func (exf) if not explicitly defined $refnum
if [[ -z $refnum ]] ; then
    nframes=`fslinfo $epi | awk '$1 == "dim4"{print $2}'`
    [[ $nframes -eq 1 ]] && refnum=0 || refnum=`echo "$nframes/2-1" | bc `
    echo "nframes = $nframes, refnum = $refnum" | tee -a $LF
fi

Run "fslroi $epi $exf $refnum 1"

    # See if there's a .mat file to propogate
epimat=$(basename $epi | sed 's/nii.gz$/mat/')
[[ -e $epimat ]] || epimat=""

# Keep only the first frame from mag
Run "fslroi $mag $TMPDIR/mag 0 1"
mag=$TMPDIR/mag


# Create brain mask from the mag
Run "bet $mag $brain -m -f $betf"
brainmask=${brain}_mask
magbrain=$brain

# Create head mask by dilating the brain mask 3 times
Run "fslmaths $brainmask -dilM -dilM -dilM $head"



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# now DIFFERENT for multi-channel  KL
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Do the phase unwrapping of phase difference (-f for 3D, -v for verbose)

nphase=${nphase:-8}
Run "prelude -c $dph -o $TMPDIR/dph -n $nphase -f -v" #-m $head);
dph=$TMPDIR/dph

    # FUGUE wants a phase image for each echo, but we only have
    # phase difference between echoes. So create an image of 0s 
    # and merging with the phase diff.
ph1=$TMPDIR/ph1
Run "fslmaths $dph -mul 0 $ph1"

    #currrent version (FSL3.3) prelude does not preserve header information
    #thus we need to add header to ph1 and dph, copy from brain mask. KL
    #fslcpgeom $brainmask $ph1
    #fslcpgeom $brainmask $dph
    # current version of FSL (4.1.9) does

    # Merge, baby, merge
ph2=$TMPDIR/ph2
Run "fslmerge -t $ph2 $ph1 $dph"

    # Create the voxel shift map (VSM) in the mag/phase space. Use mag as 
    # input to assure that VSM is same dimension as mag. The input only affects
    # the output dimension. The content of the input has no effect
    # on the VSM. The dewarped mag volume is meaningless and will be thrown away. 
vsmmag=$TMPDIR/vsmmag
magdw=$TMPDIR/magdw # To be thrown away
Run "fugue -i $mag -u $magdw -p $ph2 --dwell=$esp --asym=$tediff --mask=$brainmask  --saveshift=$vsmmag $UNWARPDIRIFSET $SAVEFMAPIFSET" 

# Forward warp the mag in order to reg with func
# What does mask do here?
magfw=$TMPDIR/magfw
Run "fugue -i $mag -w $magfw --loadshift=$vsmmag  --mask=$brainmask $UNWARPDIRIFSET"


### ==================================================
### MD:  $magfw is no good to use in registration. 
### Use mag brain instead:

# Register magfw to example func. There are some parameters here
# that may need to be tweeked.
#Run "flirt -in $magfw -ref $exf \
#  -out $TMPDIR/magfw-in-exf \
#  -omat $TMPDIR/magfw-in-exf.fsl.mat \
#  -dof 6"

# --------------------------------------------------
#  -bins 256 -cost corratio -interp trilinear \  # <-- this is default
#  -nosearch \                 # <-- nosearch and searchr[xyz] are conflicting options. Just let it search
#  -searchrx -10 10 \
#  -searchry -10 10 \
#  -searchrz -10 10 \
### ==================================================

mag2exf=$TMPDIR/magfw-in-exf.fsl.mat
Run "flirt -in $magbrain -ref $exf -omat $mag2exf -dof 6 "


# Now resample VSM into epi space. This will take care of any 
# differences in in-plane voxel size.
Run "flirt -in $vsmmag -ref $exf -out $vsm \
  -init $mag2exf -applyxfm -paddingsize 1"


# KL added: Now align brain mask to epi space as well.
Run "flirt -in $brainmask -ref $exf -out $TMPDIR/brain1\
  -init $mag2exf -applyxfm -paddingsize 1"

brainmask=$TMPDIR/brain1

# Propogate mat file:
[[ -z $epimat ]] || cp $epimat $(basename $vsm | sed 's/nii.gz/mat/')

# Check whether we can stop at this point
if [[ -z $exfdw ]] && [[ -z $epidw ]] ; then
    echo "Done0. "
    exit 0
fi


# UNWARP ONLY THE EXAMPLE FRAME HERE
if [[ ! -z $exfdw ]] ;  then
  # Now apply the VSM to the exf (um=unmasked)
    exfdwum=$TMPDIR/exfdwum
    Run "fugue -i $exf -u $exfdwum --loadshift=$vsm  --mask=$brainmask $UNWARPDIRIFSET"

  # Now mask the dewarped exf using the mag brain mask.
  # This is necessary to prevent voxels in the original EPI
  # that are out of the brain from simply being copied into
  # the dewarped image.
    Run "fslmaths $exfdwum -mas $brainmask $exfdw"
fi

 # Propogate mat file:
[[ -z $epimat ]] || cp $epimat $(basename $exfdw | sed 's/.nii.gz/mat')


 # UNWARP ALL EPI FRAMES HERE 
if [[ ! -z $epidw ]] ; then
    if [[ x$domoco = x1 ]] && [[ x$epimerged = x1 ]] ; then
	mepi=$TMPDIR/mepi
	Run "mcflirt -in $epi -out $mepi -refvol $nmidframe -plots"
	epi=$mepi
    fi
fi

for frame in $(seq 1 $nframes) ; do 
    inimg=$TMPDIR/inframe
    roiframe=$(($frame - 1))
    Run "fslroi $epi $inimg $roiframe 1"
    outimg=`printf "$TMPDIR/uavwepi%04d" $frame`
    echo inimg $inimg outimg $outimg

    Run "fugue -i $inimg -u $outimg --loadshift=$vsm --mask=$brainmask $UNWARPDIRIFSET"

  #OPTIONAL MASKING OF THE DATA
    if [[ x$do_outmask = x1 ]] ; then
	outmaskname=m  #for final merge filename 
	moutimg=`printf $TMPDIR/muavwepi%04d $frame`
	Run "fslmaths $outimg -mas $brainmask $moutimg"
    fi

  # Need to add propagation of mat file! 
done

## and merge the results
Run "fslmerge -t $epidw $TMPDIR/${outmaskname}uavwepi* "

[[ x$cleanup_forced = x1 ]] && rm -r $TMPDIR

echo "Done."
exit 0;
#####################################

#---- Everything below here is printed out as part of help -----#
cat <<__EOHELP

BEGINHELP

SUMMARY

Front end for FSLs PRELUDE and FUGUE programs to correct for B0
distortion in functional EPI scans. The programs use a B0 fieldmap.
This is assumed to be two conventional GREs collected at two different
echo times (TEs). The field map should be acquired with the same
slice prescription, slice thickness, slice skip, and field-of-view
as the EPI. The field map can have a higher in-plane resolution.

For the stock Siemens field map, two field maps are required, one to
get the phase and another to get the magnitude. The volume that is
stored as the phase is actually the phase difference and is scaled
between 0 to 4095 for -pi to pi. The magnitude volumes for both 
TEs are saved, but only one is needed.

All volumes are assumed to be in NIFTI 4D format. All volumes should
be referred to as volid.$postfix. If the EPI volume has a .mat file, this
will be propagated to the outputs.

ALGORITHM OVERVIEW

1. Create a brain mask from the mag volume (BET)
2. Create a head mask by dilating the brain mask
3. Rescale the phase image to -pi to pi
4. Unwrap the phase (PRELUDE)
5. Create the voxel shift map (VSM) (FUGUE)
6. Forward warp the mag volume (FUGUE)
7. Register the forward warped mag with the example func (FLIRT)
8. Resample the VSM into EPI space (FLIRT)
9. Dewarp the EPI and/or Example Func (FUGUE).
10. Mask out-of-brain voxels using the mag brain mask.

The EPI and Example Func should be in register with the mag volume.

ARGUMENTS:

--mag magvolid

Magnitude volume from the B0 fieldmap. If more than one frame is
present, then only the first is used.

--dph phasediffvolid

Phase difference volume (Echo2-Echo1). These are assumed to be scaled
between 0 to 4095 for -pi to pi. Eg, dph.nii.gz

--epi epivolid

EPI volume to dewarp (e.g. epivol.nii.gz)

--sepi episeries

EPI series to dewarp. (e.g. avwepi) 

--tediff tediff

Difference in the echo times of the B0 field map in ms. This number is
set when acquiring the field map. This is often set to 2.44 ms at 3T,
which is the amount of time it takes for the fat and water to
rephase. Scales with field strength.

--esp echospacing

Time (ms) between the start of the readout of two successive lines in
k-space during the EPI acquisition.  This can also be thought of as
the time between rows.  This parameter is not available in the Siemens
DICOM header. If one saves the raw kspace data, then the echo spacing
can be found next to the m_lEchoSpacing tag in the meas.asc and is in
usec. Note that the echo spacing is not the same as the time to read
out a line of k-space and it cannot be computed from the bandwidth
because neither of these methods takes into account the dead time
between lines when no data are being read out.

--vsm vsmvolid

Voxel shift map. This is a volume the same size as the EPI. The value
at each voxel indicates the amount the voxel has shifted in the phase
encode direction due to B0 distortion. This argument is required.

--fmap fmpvolid

fieldmap volume. This is a volume the same size as the EPI. The value
at each voxel indicates the B0 distortion (rad/s). This argument is optional.

--exfdw exfdwvolid

This is the middle time point from the EPI time series (the "example
functional") with B0 distortion removed. Mainly good for checking how
good the dewarping is without having to dewarp the entire time
series. NOTE: the voxels outside the brain are set to 0.

--epidw epidwvolid

This is the EPI time series with B0 distortion removed. NOTE: the
voxels outside the brain are set to 0.

--tmpdir tmpdir

Location to put the directory for storing temporary files. By default,
this will be in a directory called tmp-epidewarp.fsl under the
directory to hold the VSM volume. When --tmpdir is used or --nocleanup
is specified, then this directory will not be deleted. Otherwise or
with --cleanup it will automatically be deleted.

--nomoco

Disables motion correction of functionals prior to fugue'ing.

--refnum <refnum>

Defines reference image number for motion correction. Default is middle image in a series. 

--nocleanup

Do not delete the tmp dir. --tmpdir automatically implies --nocleanup.

--cleanup

Forces deleting of the tmp dir regardless of whether --tmpdir or 
--nocleanup have been specified.

--nphase

Number of the phase splits to use in prelude. 

--debug

Prints copious amounts to the screen.


BUGS

1. Currently does not work when the B0 field map has a different
in-plane resolution than the EPI.

AUTHORS

Doug Greve and Dave Tuch (help@nmr.mgh.harvard.edu) with generous 
help from the FSL crew.

UCSD specific mods made by Tom Liu (ttliu@ucsd.edu) and Giedrius Buracas (giedrius@salk.edu) and Kun Lu (kunlu@ucsd.edu)
__EOHELP
