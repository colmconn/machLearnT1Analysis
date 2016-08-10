#! /bin/bash
## reviewscans.bash LISTFILE ["overwrite"]
trap exit SIGHUP SIGINT SIGTERM
## This would have been a lot easier in python, in retrospect.

# David Perlman 2012 Jan 8: reviewscans.bash
# uses afni to save 2d images and then puts them into an html table for quick review
# uses xvfb, which is installed by default on our OS X servers, but not on the linux boxes.
# Site-specific note: actually this runs faster on guero than jaloro, even though
# it has to open all those windows on the screen.
# So consider the xvfb option to only serve the purpose of running this through screen, say.

# If the outdir you give it already exists, it will make something like outdir_01

programName=`basename $0`

GETOPT=$( which getopt )

GETOPT_OPTIONS=$( $GETOPT  -o "onr:" --longoptions "overwrite,nosubdirs,resize:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi
resizeTo=300
noSubDirs=0
# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-o|--overwrite)
	    do_overwrite=1; shift ;;
	-n|--nosubdirs)
	    noSubDirs=1; shift ;;
	-r|--resize)
	    do_resize=1; resizeTo=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

# set the interesting stuff here
# For debugging:
always_show_viewer=0

# the command to use for afni, set here if you need a specific version
AFNI="env AFNI_LAYOUT_FILE=ElvisIsAliveOnPlanetZork afni -noplugins -no_detach"
# AFNI="/apps/x86_64_sci6/AFNI_2011_12_21_1014_12Mar20/afni"

# name to call the output images
nameroot=rendered_image

# the input specification file
specification_file="$1"

# the log file so I can spy on everyone using this program
#logfile=/home/perlman/bin/${0%.*}_log.txt

# for the record later
original_cmd_line="$0 $*"

# the BigX data set to sub for missing data sets
# later make this into a relative path from the script file
#bigxorig="/home/perlman/bin/BigX/BigXb+orig.nii"
#bigxtlrc="/home/perlman/bin/BigX/BigXb+tlrc.nii"

# put a default value in every keyword, in case it's blank it doesn't throw off the count
#underlay="MNI_avg152T1+tlrc"
underlay=None
colorscale=Spectrum:yellow_to_cyan
opacity=9
range=0
view=tlrc

# # We can tell it to overwrite the output instead of incrementing a number
# if [[ "$2" == "overwrite" ]]; then
#   do_overwrite=1
# else
#   do_overwrite=0
# fi


function add_brik_extensions () {
  # This is to workaround a longstanding bug in afni viewer direct file input
  local output_file_list
  local thefile
  local theext
  while [[ ! -z "$1" ]]; do
    theext="${1##*.}" # just the extension
    if [ "$theext" == "gz" ] || [ "$theext" == "bz2" ] || [ "$theext" == "Z" ] ; then
	## remove the compression extension
	tf="${1%.*}"
	theext="${tf##*.}" # just the extension
    fi

    #theext="${theext^^}"  # all capital letters
    theext=$( echo $theext | tr "[:lower:]" "[:upper:]" )
    # test if it's one of the extensions we like
    # if it's not, assume we need to add BRIK
    # if it has some other extension, then we're just screwed
    valid_exts="BRIK HEAD NII HDR IMG"
    if [[ "$valid_exts" =~ "$theext" ]]; then
      # do nothing
      thefile=$1
    else
      # add the extension
      thefile=${1}.BRIK
    fi
    output_file_list="$output_file_list $thefile"
    shift
  done
  echo $output_file_list
}
    


get_output_dir_name() {
  # If I wanted to be really clever, I could first strip off any _## that was already in the name
  # but that's probably overkill
  if [[ -e "$1" ]]; then
    local num=1
    local name=${1}_$(printf "%02d" $num)
    while [[ -e "$name" ]]; do
      (( num++ ))
      name=${1}_$(printf "%02d" $num)
    done
    outdir="$name"
  else
    outdir="$1"
  fi
}

function setup_output_dir() {
  if (( do_overwrite != 0 )); then
    echo Overwriting-DELETING!-output directory $1
    rm -Rf $1
  fi
  # This sets the variable outdir
  get_output_dir_name "$1"
  # Make the output directory
  #local output_dir="$1"
  echo Making output directory $outdir
  mkdir "$outdir"
  #cd "$output_dir"
  # Keep a copy of the source file for posterity
  rm -f ${outdir}/specification_file.txt.old
  mv ${outdir}/specification_file.txt ${outdir}/specification_file.old.txt 2>/dev/null
  #echo cp "$specification_file" ${outdir}/specification_file.txt
  specification_file_copy="${outdir}/specification_file.txt"
  cp "$specification_file" ${specification_file_copy}
  # set up the output html file
  htfile="${outdir}/index.html"
}



function begin_html_file() {
  # $1 is title for the output
  # no wait, we'll use global reviewTitle instead
  # htfile is global
  # specification_file is global
  #local reviewTitle=$1
  echo Creating HTML file $htfile
  echo "<head>" >>${htfile}
  echo "<title>${reviewTitle}</title>" >>${htfile}
  echo "</head>" >>${htfile}
  echo >>${htfile}
  echo "<body>" >>${htfile}
  #echo "<h2>${reviewTitle} ${source_file} $(date)</h2>" >>${htfile}
  echo "<br>" >>${htfile}
  echo "<table cellspacing=10 border=1>" >>${htfile}
  echo "<caption><big><bold>${reviewTitle} <br>${source_file} $(date)</big></bold><br>" >>${htfile}
  echo "<small><a href=specification_file.txt>Specification file that generated this review</a></caption>" >>$htfile
  echo >>${htfile}
}



function write_html_file_headers () {
  # renderimage is a global here
  # it is an array.  each element of the array looks like this:
  # underlay overlay colorscale opacity range threshtype view "RENDERIMAGE" column_title cbrik tbrik sliceaxis thresh x y z
  # 16 items for each renderimage item.
  local numrenders
  numrenders=${#renderimage[@]}
  
  echo "<tr>" >>${htfile}
  echo "<th> Subject </th>" >>$htfile

  #local rownumj=$(printf "%03d" $1)
  #shift
      
  local imcount=1
  
  # This is actually a lot simpler now that I made the title a separate field
  # We just need the title, 
  while [[ ! -z "$1" ]]; do
    local underlay=$1
    shift
    local overlay=$1
    shift
    local colorscalename=$1
    shift
    local opacity=$1
    shift
    local range=$1
    shift
    local threshtype=$1
    shift
    local view=$1
    shift
    # item 8 is just the word RENDERIMAGE
    shift
    local column_label=$1
    shift
    local cbrik=$1
    shift
    local tbrik=$1
    shift
    local sliceaxis=$1
    shift
    local thresh=$1
    shift
    local x=$1
    shift
    local y=$1
    shift
    local z=$1
    shift
    
    # Now assemble this column's header
    
#     if [[ "$overlay" == "NO_OVERLAY" ]]; then
#       overlay="No Overlay"
#     elif [[ "${overlay:0:10}" == "NO_OVERLAY" ]]; then
#       if [[ "${overlay:10:1}" =~ [[:punct:]] ]]; then
#         overlay=${overlay:11}
#       else
#         overlay=${overlay:10}
#       fi
#     ####################################
#     # Put something in here to let us use ::: as a title indicator
#     # Will also need to get rid of it for the actual filename
#     elif [[ "${overlay}" =~ ":::" ]]; then
#       overlay=${overlay#*:::}
#     fi
    
    local theHeader
    #theHeader="<center><small>${overlay}</small> <br> ColorBrik $cbrik ThreshBrik $tbrik <br> $threshtype thresh=$thresh x=$x y=$y z=$z</center>"
    theHeader="<center><small>${column_label}</small> <br> ColorBrik $cbrik ThreshBrik $tbrik <br> $threshtype thresh=$thresh x=$x y=$y z=$z</center>"
    echo "<td> $theHeader </td>" >>$htfile
    (( imcount++ ))
  done
  
  echo "</tr>" >>${htfile}
  echo         >>${htfile}
  
}
  

function write_html_file_row() {
  # make_images_outnames is the current list of filenames to go in the row
  echo "<tr>" >>$htfile
  echo "<td> $* </td>" >>$htfile
  local anImage
  for anImage in $make_images_outnames; do
    echo "<td><img src=\"${anImage}\"></td>" >>$htfile
  done
  echo >>$htfile
}
  


function close_html_file () {  
  # close out the html file
  echo "</table>" >>${htfile}
  echo "</body>" >>${htfile}
}



function dataset_id () {
  #local dsetview=$(3dinfo -av_space $1 2>/dev/null)
  { 3dinfo $1 2>/dev/null || 3dinfo $bigxtlrc 2>/dev/null; } | grep 'Identifier Code:' | awk '{print $3}'
}
#{ 3dinfo WB102_m_T1High+tlrc.niig 2>/dev/null || 3dinfo WB102_m_T1High+tlrc. 2>/dev/null; } | grep 'Identifier Code:' | awk '{print $3}';


#function replace_missing_dataset () {
  # Check if the given dataset $1 exists;
  # if not, replace it with the BigX data set
  


function make_images() {
  # We will read from global variables: outimages, displayoption, nameroot
  # Limitation to add later: discrete colorscale? xhairs?
  # colorscalename should start out as Spectrum:yellow_to_cyan
  # overall command:
  # underlay overlay colorscale opacity range threshtype view "RENDERIMAGE" column_title cbrik tbrik sliceaxis thresh x y z
  # 16 parameters per cycle; rownum only appears once at the beginning
  
  local rownumj=$(printf "%03d" $1)
  shift
  
  local commandstring=" \
    -com \"CLOSE_WINDOW axialimage\" \
    -com \"CLOSE_WINDOW sagittalimage\" \
    -com \"CLOSE_WINDOW coronalimage\" \
    -com \"CHDIR ${outdir}\" \
    -com \"SET_XHAIRS OFF \" "
    
  local imc=1
  make_images_outnames=""
  local inputfilelist=""
  
  while [[ ! -z "$1" ]]; do
    #echo Remaining arguments: $*
    imcount=$(printf "%03d" $imc)
    local outname=${nameroot}_${rownumj}_${imcount}   ### This won't work, numbers need to increment across runs of this function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # add the output file to the list, so we can put it into the html file
    # we're using png here
    make_images_outnames="$make_images_outnames ${outname}.png"
    
    local underlay_path=$1
    inputfilelist="$inputfilelist $underlay_path"
    echo underlay $underlay_path $(dataset_id $underlay_path)
    local underlay_id=$(dataset_id $underlay_path)
    #local underlay=$(basename $underlay_path)
    shift
    local overlay_path=$1
    inputfilelist="$inputfilelist $overlay_path"
    echo overlay $overlay_path $(dataset_id $overlay_path)
    local overlay_id=$(dataset_id $overlay_path)
    #local overlay=$(basename $overlay_path)
    shift
    local colorscalename=$1
    shift
    local opacity=$1
    shift
    local range=$1
    shift
    local threshtype=$1
    if [[ "$threshtype" == "raw" ]]; then
      threshtype=""
    fi
    shift
    local view=$1
    shift
    # this item is just the word "RENDERIMAGE"
    shift
    local column_label=$1
    echo Column Label $column_label
    shift
    local cbrik=$1
    shift
    local tbrik=$1
    shift
    local sliceaxis=$1
    shift
    local thresh=$1
    shift
    local x=$1
    shift
    local y=$1
    shift
    local z=$1
    shift
    
    # Now assemble the command string
    if [[ "${overlay:0:10}" == "NO_OVERLAY" ]]; then
      commandstring="$commandstring \
        -com \"SET_VIEW $view\" \
        -com \"SET_ANATOMY $underlay_id\" \
        -com \"SEE_OVERLAY -\" \
        -com \"SET_DICOM_XYZ $x $y $z\" \
        -com \"OPEN_WINDOW $sliceaxis opacity=$opacity\" \
        -com \"SAVE_PNG $sliceaxis ${outname}\" "
      #-com \"SAVE_FILTERED $sliceaxis 'cat | convert ppm:- -resize 300% ${outname}.png'\" "
      
    else
      commandstring="$commandstring \
        -com \"SET_VIEW $view\" \
        -com \"SET_ANATOMY $underlay_id\" \
        -com \"SET_PBAR_ALL -99 1.0 $colorscalename\" \
        -com \"SET_FUNC_RANGE $range\" \
        -com \"SEE_OVERLAY +\" \
        -com \"SET_FUNCTION $overlay_id $cbrik $tbrik\" \
        -com \"SET_DICOM_XYZ $x $y $z\" \
        -com \"SET_THRESHNEW $thresh $threshtype\" \
        -com \"OPEN_WINDOW $sliceaxis opacity=$opacity\" \
        -com \"SAVE_PNG $sliceaxis ${outname}\" "
	#-com \"SAVE_FILTERED $sliceaxis 'cat > ${outname}.ppm'\" "
	#-com \"SAVE_FILTERED $sliceaxis 'cat | convert ppm:- -resize 300% ${outname}.png'\" "

    fi
    (( imc++ ))
  done
  
  # Terminate the command string
  commandstring="$commandstring -com QUIT "
    
  #########################################################################
  # New stuff 2010-01-31:
  # need to figure out what directory we're in so this makes sense
  # need to make sure $inputfilelist has ${subj} expanded and so on
  # need to change command so afni reads the inputfilelist files
  # need to change command so failed loads leave the X mark
  
  # Make sure the big X files are loaded
  inputfilelist="$inputfilelist $bigxorig $bigxtlrc"
  # Now get rid of duplicates in inputfilelist
  inputfilelist=$(echo $inputfilelist | tr ' ' '\n' | sort | uniq | tr '\n' ' ')
  
  # Now add BRIK extensions if necessary, to work around stupid bug in afni viewer
  inputfilelist=$(add_brik_extensions $inputfilelist)
  
  # Now run the command.  I don't know if this is a reasonable thing to do,
  # or a horrible hack!
  # echo "afni $displayoption $commandstring"
  echo Running afni command with input file list:
  #echo $inputfilelist
  fullcommand="$AFNI $displayoption $commandstring $inputfilelist"
  #echo "afni $displayoption $commandstring $inputfilelist"
  #eval "afni $displayoption $commandstring $inputfilelist"
  echo $fullcommand
  eval $fullcommand
  
  # Now the png images should exist.
  
}





# End of function definitions.  Now run the process.




set it up for afni to be able to run headless
if [[ "$always_show_viewer" == 1 ]]; then
  echo Forcing live X display
else
  if which xvfb 2> /dev/null; then
    xvfb :1 -screen 0 800x800x24 &
    displayproc=$!
    # that saves the PID of the xvfb process so we can kill it later
    displayoption="-display :1.0"
  elif which Xvfb 2> /dev/null; then
    Xvfb :1 -screen 0 800x800x24 &
    displayproc=$!
    # that saves the PID of the xvfb process so we can kill it later
    displayoption="-display :1.0"
  elif which /usr/bin/xvfb 2> /dev/null; then
    /usr/bin/xvfb :1 -screen 0 800x800x24 &
    displayproc=$!
    # that saves the PID of the xvfb process so we can kill it later
    displayoption="-display :1.0"
  elif which /usr/bin/Xvfb 2> /dev/null; then
    /usr/bin/Xvfb :1 -screen 0 800x800x24 &
    displayproc=$!
    # that saves the PID of the xvfb process so we can kill it later
    displayoption="-display :1.0"
  else
    echo "headless operation not possible on this machine."
    displayoption=""
    displayproc=""
  fi
fi


echo Using afni viewer to render images based on list $1
rownum=0
need_headers=1
while read theline; do
  linearray=( ${theline} )    # get it into an array for convenience
  linecode=${linearray[0]}
  linerestarray="${linearray[@]:1}"
  linerest="${linerestarray[*]}"
  
  # skip comment lines
  if [[ "${linecode:0:1}" == "#" ]]; then continue; fi
  
  # skip blank lines
  if [[ -z "${linecode}" ]]; then continue; fi
  
  # read in keyword info
  case $linecode in 
    # Put in defaults here in case people give blank keyword values???
    TITLE )
      reviewTitle=$linerest
      echo Review Title will be $reviewTitle
      continue ;;
    DATADIR )
      datadir="${linerest}"
      curr_data_dir=${datadir}/
      echo Data source directory is $datadir
      continue ;;
    SUBDIR )
      subdir="${linerest}"
      echo "Data source subdirectory (after the subject code) is $subdir"
      continue ;;
    OUTDIR )
      # Set variable outdir
      # Set variable htfile
      setup_output_dir "$linerest"
      echo "Output directory is $outdir"
      echo "HTML file is $htfile"
      begin_html_file
      continue ;;
    UNDERLAY )
      underlay="$linerest"
      #underlay=${subdir}/${underlay}/
      echo Underlay file is $underlay
      # now assemble the underlay name as a filename so afni can load it
      #inputfiles="$inputfiles $underlay"
      continue ;;
    OVERLAY )
      overlay="$linerest"
      #overlay=${subdir}/${overlay}/
      echo Overlay file is $overlay
      if [[ ! $do_resize ]] ; then
	  echo "*** You might want to use the resize option with an actual overlay as the images can be quite small"
      fi
      # now assemble the overlay name as a filename so afni can load it
      # inputfiles="$inputfiles $overlay"
      continue ;;
    COLORSCALE )
      colorscale="$linerest"
      echo Colorscale to use is $colorscale
      continue ;;
    OPACITY )
      opacity="$linerest"
      echo Opacity of functional display is $opacity
      continue ;;
    THRESHTYPE )
      threshtype="$linerest"
      echo Threshold type is $threshtype
      continue ;;
    RANGE )
      range="$linerest"
      echo Functional display range is $range
      continue ;;
    VIEW )
      view="$linerest"
      echo Display view is $view
      continue ;;
    RESETIMAGE )
      renderimage=""
      #inputfiles=""
      need_headers=1
      echo "Resetting render parameters, hopefully we'll get some new ones soon"
      continue ;;
    RENDERIMAGE )
      # When we get here, we expect the following paths to be valid already:
      # underlay, overlay, datadir, subdir
      # but we can't assemble the current files because we don't know $SUBJ yet
      # but we can cleverly put it in as a variable
	# Assemble the current files
	if [[ $noSubDirs -eq 1 ]] ; then
	    current_underlay_path=${datadir}/${underlay}
	    current_overlay_path=${datadir}/${overlay}
	else
	    current_underlay_path=${datadir}/'${SUBJ}'/${subdir}/${underlay}
	    current_overlay_path=${datadir}/'${SUBJ}'/${subdir}/${overlay}
	fi
      # Now those variables will get filled in at the same time as any that the user put in the template
      # Now we need to be tricky to keep track of two things:
      # 1. the list of input file full paths, which provide the inputs to afni command line
      # 2. the input file codes, which go in the actual renderimage spec?
      #renderimage="$renderimage $linerest"
      #linearray[2]=${datadir}/${linearray[2]}/${subdir}
      current_renderimage="$current_underlay_path $current_overlay_path $colorscale $opacity $range $threshtype $view $theline"
      # underlay overlay colorscale opacity range threshtype view "RENDERIMAGE" column_title cbrik tbrik sliceaxis thresh x y z
      # 16 items for each renderimage item.
      # renderimage="$renderimage $underlay $colorscale $opacity $range $threshtype $view ${linearray[@]:1}"
      # build renderimage array of individual renderimage specs
      # renderimage=( "${renderimage[@]}" "$current_renderimage" )
      renderimage="$renderimage $current_renderimage"
      #echo "Adding image render parameters ${linearray[@]:1}"
      #overlayspec=${linearray[2]}
      # inputfiles="$inputfiles $overlayspec"
      continue ;;
    EMAIL )
      finished_email_address="$linerest"
      continue ;;
    BREAK )
      break ;;
      # hopefully this will just jump us out of this loop, and close out the web page
#     * )
#       continue ;;
  esac
  
  # Now we've processed the keywords.  If we get here, it means it's a subject.
  # This means it is a subject to process, if we get here.
  (( rownum++ ))
  echo row $rownum processing ${linecode}
  # Need to put together the command string
  # What can be different from one to the next?
  # underlay colorscale threshtype range view 
  #image_command="$renderimage $underlay $colorscale $threshtype $range $view $linerest"
  #make_images $nameroot $underlay $parameters
  #echo make_images $rownum $renderimage
  # Note that this means all images have to be in the same directory for each row
  
  # We probably don't need to do this, because I am going to calculate all full paths
  # cd ${datadir}/${linecode}/${subdir}/
  
  # OK.  To render each row of images, we need:
  # rownum, subj, theline (subj plus whatever other title), full renderimage
  
  # Here we need to take into account the ${subj} notation in the input file
  SUBJ=$linecode
  renderimage_withsubj=$(eval "echo $renderimage")
  
  # If we just started rendering, then write headers.  Otherwise just go ahead with the rendering images.
  if (( need_headers )); then
    write_html_file_headers  $renderimage_withsubj # $rownum not needed? $renderimage will be accessed as a global, since I don't know about passing an array to a function in bash
    need_headers=0
  fi

  #echo make_images $rownum $renderimage_withsubj
  make_images $rownum $renderimage_withsubj
  write_html_file_row $theline

done <$specification_file

close_html_file

cd $outdir

if (( do_resize != 0 )) ; then
    echo "Resizing images"
    for f in rendered* ; do
	convert $f -resize ${resizeTo}% new-$f
	mv -f new-$f $f
    done
    echo "Resizing images. Done."
fi

# kill the headless display server
if [[ -n "$displayproc" ]]; then
  kill $displayproc
fi


# Send an email if necessary

if [[ -n "$finished_email_address" ]]; then
  echo "$original_cmd_line done" | mail -s "$original_cmd_line" $finished_email_address
fi

