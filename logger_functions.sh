#!/bin/bash

green='\e[0;37;42m'
yellow='\e[0;37;43m'
red='\e[0;37;41m'
endColor='\e[0m'

prefix_tag="***"

function message {
    echo "${prefix_tag} ${message}"
}

function info_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${green}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag}${message}"
    fi
}


function warn_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${yellow}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag}${message}"
    fi
}

function error_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${red}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag}${message}"
    fi
}

#info_message  "This is a test info message"
#warn_message  "This is a test warning message"
#error_message "This is a test error message"
