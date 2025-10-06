#!/bin/bash

# This script runs on Linux and MaxOS and downloads all the selected files to the current working directory in up to 5 parallel download streams.
# Should a download be aborted just run the entire script again, as partial downloads will be resumed. 
# We've finally got our server-side software sorted out and now positively encourage the use of parallel downloads
# of single files from version 2020-OCT.

# connect / read timeout for wget / curl
export TIMEOUT_SECS=100
# how many times do we want to automatically resume an interrupted download?
export MAX_RETRIES=5
# after a timeout, before we retry, wait a bit. Maybe the servers were overloaded, or there was some scheduled downtime.
# with the default settings we have 15 minutes to bring the dataportal service back up.
export WAIT_SECS_BEFORE_RETRY=30
# the files to be downloaded
LIST=("
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1d5_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1dc_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1de_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1e3_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1e5_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ea_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ec_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1f1_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1f8_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1fa_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ff_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X206_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X20d_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X214_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X216_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X218_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X21b_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X222_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X229_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X22b_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X230_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X232_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X237_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X23e_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X240_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_X97_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xac_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xc1_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xc8_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xcf_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xd6_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xdd_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xe4_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xe6_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xeb_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xf2_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xf9_001_of_001.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1d5_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1dc_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1de_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1e3_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1e5_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ea_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ec_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1f1_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1f8_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1fa_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X1ff_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X206_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X20d_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X214_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X216_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X218_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X21b_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X222_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X229_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X22b_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X230_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X232_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X237_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X23e_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X1467_X240_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_X97_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xac_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xc1_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xc8_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xcf_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xd6_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xdd_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xe4_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xe6_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xeb_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xf2_auxiliary.tar
https://almascience.eso.org/dataPortal/2019.1.00195.L_uid___A001_X146c_Xf9_auxiliary.tar
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1d5.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1dc.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1de.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1e3.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1e5.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1ea.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1ec.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1f1.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1f8.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1fa.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X1ff.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X206.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X20d.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X214.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X216.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X218.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X21b.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X222.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X229.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X22b.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X230.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X232.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X237.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X23e.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X1467_X240.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_X97.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xac.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xc1.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xc8.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xcf.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xd6.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xdd.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xe4.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xe6.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xeb.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xf2.README.txt
https://almascience.eso.org/dataPortal/member.uid___A001_X146c_Xf9.README.txt
")
# If we terminate the script using CTRL-C during parallel downloads, the remainder of the script is executed, asking if
# the user wants to unpack tar files. Not very nice. Exit the whole script when the user hits CTRL-C.
trap "end_session && exit" INT

# quickly log in to the server (if required) and create a cookies file. This is because some users are using shared machines for long-running
# downloads where the credentials can be seen in plain-text in the task list.
function start_session {
  local username=anonymous
  export authentication_status=0
  if [ "${username}" != "anonymous" ]; then
    # only prompt for the password again if we haven't already. What could be going on here? The servers
    # have been restarted and the local cookie file is no longer valid.
    if [ -z ${password} ]; then
      echo ""
      echo -n "Please enter the password for ALMA account ${username}: "
      read -s password
      echo ""
      export password
    fi

    # ICT-17647 curl v7.29.0 is very slow when the --speed-limit option is used. Use wget by default until ALMA 
    # uses a more recent version of RHEL.
    if command -v "wget" > /dev/null 2>&1; then
      login_command=(wget --quiet --delete-after --no-check-certificate --auth-no-challenge --keep-session-cookies --save-cookies alma-rh-cookie.txt "--http-user=${username}" "--http-password=${password}")
    elif command -v "curl" > /dev/null 2>&1; then
      login_command=(curl -s -k -o /dev/null -c alma-rh-cookie.txt --fail "-u" "${username}:${password}")
    fi
    # echo "${login_command[@]}" "https://almascience.eso.org/dataPortal/login"
    $("${login_command[@]}" "https://almascience.eso.org/dataPortal/login")
    authentication_status=$?
    if [ $authentication_status -eq 0 ]; then
      echo "	    OK: credentials accepted."
    else
      echo "	    ERROR: login failed. Error code is ${authentication_status}"
    fi
  fi
}
export -f start_session

# clean up the cookies file after downloading is complete
function end_session {
	# Included unset variable: password, for ticket : ICT-16332
	unset password
	rm -fr alma-rh-cookie.txt
}
export -f end_session

export failed_downloads=0

# download a single file.
# attempt the download up to N times
function dl {
  local file=$1
  local filename=$(basename $file)
  # the nth attempt to download a single file
  local attempt_num=0

  # wait for some time before starting - this is to stagger the load on the server (download start-up is relatively expensive)
  sleep $[ ( $RANDOM % 10 ) + 2 ]s

  if command -v "wget" > /dev/null 2>&1; then
    local tool_name="wget"
    local download_command=(wget -c -q -nv --no-check-certificate --auth-no-challenge --load-cookies alma-rh-cookie.txt --timeout=${TIMEOUT_SECS} --tries=1)
  elif command -v "curl" > /dev/null 2>&1; then
    local tool_name="curl"
    local download_command=(curl -C - -S -s -k -O -f -b alma-rh-cookie.txt --speed-limit 1 --speed-time ${TIMEOUT_SECS})
  fi

  # manually retry downloads. 
  # I know wget and curl can both do this, but curl (as of 10.04.2018) will not allow retry and resume. I want consistent behaviour so 
  # we implement the retry mechanism ourselves.
  echo "starting download of $filename"
  until [ ${attempt_num} -ge ${MAX_RETRIES} ]
  do
    # echo "${download_command[@]}" "$file"
    $("${download_command[@]}" "$file")
    status=$?
    # echo "status ${status}"
    if [ ${status} -eq 0 ]; then
      echo "	    successfully downloaded $filename"
      break
    else
      failed_downloads=1
      # users requested a string instead of a simple exit code - attempt to make it look like the curl output
      if [ ${status} -eq 8 ] && [ ${tool_name} = "wget" ]; then 
        echo "wget: (8) Server issued an error response."
      elif [ ${status} -eq 4 ] && [ ${tool_name} = "wget" ]; then
        echo "wget: (4) Network failure."
      elif { [ ${status} -eq 6 ] && [ ${tool_name} = "wget" ]; } || { [ ${status} -eq 22 ] && [ ${tool_name} = "curl" ]; }; then
        # if we restart the servers then our cookies are invalidated
        end_session
        echo "authentication error - retrying the login"
        start_session
        if [ $authentication_status -eq 0 ]; then
          # our retry attempt was hampered by a credentials issue. We shouldn't count this as an attempt. Give a bonus try...
          echo "		resuming download of $filename, still attempt $[${attempt_num}+1]"
          $("${download_command[@]}" "$file")
        fi
      fi
      #echo "${login_command[@]}" "https://almascience.eso.org/dataPortal/login"
      echo $("${download_command[@]}" "$file/client/status/$status")

      echo "		download $filename was interrupted with error code ${tool_name}/${status}"
      attempt_num=$[${attempt_num}+1]
      if [ ${attempt_num} -ge ${MAX_RETRIES} ]; then
        echo "	  ERROR giving up on downloading $filename after ${MAX_RETRIES} attempts  - rerun the script manually to retry."
      else
        echo "		download $filename will automatically resume after ${WAIT_SECS_BEFORE_RETRY} seconds"
        sleep ${WAIT_SECS_BEFORE_RETRY}
        echo "		resuming download of $filename, attempt $[${attempt_num}+1]"
      fi
    fi
  done
}
export -f dl

tars_exist=0
invalid_tars_exist=0
function check_integrity {
	for nextfile in ${LIST}; do
		tarname=$( basename ${nextfile} )
		if [[ ${nextfile} =~ ^.*\.tar$ ]] && [ -f ${tarname} ]; then
			tars_exist=1
			tar tf ${tarname} >/dev/null 2>&1
			tar_corrupted=$?
			if [[ ${tar_corrupted} -eq 1 ]]; then
				echo "tar invalid: ${tarname}"
				invalid_tars_exist=1
				break
			fi
		fi
	done
	if [[ ${invalid_tars_exist} -eq 1 ]]; then
		echo ""
		echo "Some of the downloads unfortunately failed. Please re-run the download script again in a few hours and, if this still fails, submit a ticket to the ALMA Helpdesk"
	fi
}

function unpack_tars {
	echo ""
	# override TIMEOUT so we never time-out. These files may take days to download. It's annoying if the user comes back
	# to the terminal after a weekend and finds that they weren't quick enough to automatically untar the files.
	# That's a pretty big timeout. One year in seconds.
	echo "All files have been downloaded successfully."
	echo ""
	for nextfile in ${LIST}; do
		if [[ ${nextfile} =~ ^.*\.tar$ ]]; then
			tarname=$( basename ${nextfile} )
			echo "	    unpacking ${tarname}"
			tar xf ${tarname}
			rm ${tarname}
		fi
	done
}

# temporary workaround for ICT-13558: "xargs -I {}" fails on macos with variable substitution where the length of the variable
# is greater than 255 characters. For the moment we download these long filenames in serial. At some point I'll address this issue
# properly, allowing parallel downloads.
# Array of filenames for download where the filename > 251 characters
# 251? Yes. The argument passed to bash is "dl FILENAME;" In total it cannot exceed 255. So FILENAME can only be 251
export long_files=()
# arrayf of filenames with length <= 255 characters - can be downloaded in parallel.
export ok_files=()
function split_files_list {
	for nextfile in ${LIST}; do
		length=${#nextfile}
		if [[ $length -ge 251 ]]; then
			long_files+=($nextfile)
		else
			ok_files+=($nextfile)
		fi
	done
}

# Main body
# ---------

# check that we have one of the required download tools
if ! (command -v "wget" > /dev/null 2>&1 || command -v "curl" > /dev/null 2>&1); then
   echo "ERROR: neither 'wget' nor 'curl' are available on your computer. Please install one of them.";
   exit 1
fi

echo "Downloading the following files in up to 5 parallel streams. Total size is 537GiB."
echo "${LIST}"
echo "In case of errors each download will be automatically resumed up to 3 times after a 5 minute delay"
echo "To manually resume interrupted downloads just re-run the script."
# tr converts spaces into newlines. Written legibly (spaces replaced by '_') we have: tr "\_"_"\\n"
# IMPORTANT. From the 2020-OCT release we have finally got our server-side software fixed, please now feel free to
# split downloads of a single file into multiple parallel pieces. We support range requests now (although we still
# don't support a single HTTP request with multiple range-requests)

# Included for retrials of passwords, for ticket : ICT-16332
export RETRY_ATTMPTS=5

while [ $RETRY_ATTMPTS -gt 0 ];
do
	start_session
	if [ $authentication_status -eq 0 ]; then
		RETRY_ATTMPTS=0
    else
		end_session
		RETRY_ATTMPTS=$(( $RETRY_ATTMPTS - 1 ))
    fi
done

start_session
if [ $authentication_status -eq 0 ]; then
	echo "your downloads will start shortly...."
	split_files_list
	# "dl" is a function for downloading. I abbreviated the name to leave more space for the filename
	echo ${ok_files[@]} | tr \  \\n | xargs -P5 -n1 -I '{}' bash -c 'dl {}'
	for next_file in ${long_files[@]}; do
		dl ${next_file}
	done
fi
end_session
check_integrity
if [ $failed_downloads -eq 0 ] && [ $tars_exist -eq 1 ] && [ $invalid_tars_exist -eq 0 ]; then
	unpack_tars
fi
echo "Done."

