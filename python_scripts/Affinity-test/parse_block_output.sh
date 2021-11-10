#!/bin/bash

################################################################################
# Parse block outputs from a log file that ASPECT generates
#
# Dependencies:
#
# Example Usage:
#    ./bash_scripts/parse_block_output.sh analyze_affinity_test_results 
# /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/results/spherical_shell_expensive_solver/peloton-ii-32tasks-core-openmpi-4.0.1/output_16_2_1 
# temp
################################################################################


parse_block_output(){
    ##
    # Parse value in output, looking for the block output of aspect
    # Inputs:
    #   $1(str): logfile
    # Outputs:
    #   ??: 
    #       entries:
    #           0: Total wallclock time elapsed since start
    #           1: Assemble Stokes system
    #           2: Assemble composition system
    #           3: Assemble temperature system
    #           4: Build Stokes preconditioner
    #           5: Build composition preconditioner
    #           6: Build temperature preconditioner
    #           7: Initialization
    #           8: Postprocessing
    #           9: Setup dof systems
    #           10: Setup initial conditions
    #           11: Setup matrices
    #           12: Solve Stokes system
    #           13: Solve composition system
    #           14: Solve temperature system
    local logfile="$1"
    local key="$2"
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile doesn't exist"
    # read file from the end
    local parse_results=$(eval "awk '/${key}/{print}' ${logfile} | sed 's/${key}//g' |  awk '{print \$5}' | sed ':a;N;\$!ba;s/\n/ /g'")
    return_values=("${parse_results}")
}

parse_block_output_wallclock(){
    ##
    # Parse value in output, looking for Total wallclock time
    # It turns output the previous one doesn't work for wallclock time
    local logfile="$1"
    local key="Total wallclock time elapsed since start"
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile doesn't exist"
    # read file from the end
    local parse_results=$(eval "awk '/${key}/{print}' ${logfile} | sed 's/${key}//g' |  awk '{print \$3}' | sed ':a;N;\$!ba;s/\n/ /g'")
    return_values=("${parse_results}")
}


parse_block_output_to_file(){
    ##
    # Parse value in output, looking for the block output of aspect
    # Inputs:
    #   $1(str): logfile
    #   $2(str): ofile
    #   $3-: keys
    local logfile="$1"
    local ofile="$2"
    # checkfile exist
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile(${logfile}) doesn't exist"

    local key="Total wallclock time elapsed since start";
    local header="# ${key}\n"
    parse_block_output_wallclock "${log_file}"
    output=$(echo ${return_values[@]} | sed -E "s/[^0-9.e+]/ /g")  # could be wrong when it is scientific expression
    local contents="${output}"


    # loop for key
    local i=3
    key=${!i}
    while [[ -n ${key} ]]; do
        header="${header}# ${key}\n"
        parse_block_output "${log_file}" "${key}"
        output=$(echo ${return_values[@]} | sed -E "s/[^0-9.e+]/ /g")  # could be wrong when it is scientific expression
        [[ -n ${contents} ]] && contents="${contents}\n${output}" || contents=${output}
        ((i++))
        key=${!i}
    done

    # output
    printf "${header}" > "${ofile}"
    printf "${contents}" >> "${ofile}"
}

main(){
    if [[ "$1" = "parse_block_results" ]]; then
        # this doesn't work for $3 with whitespace in it.
        parse_block_output "$2" "$3"
        printf "${return_values}"
    fi

    if [[ "$1" = "analyze_affinity_test_results" ]]; then
        local log_file="$2"
        local ofile="$3"
        parse_block_output_to_file "${log_file}" "${ofile}" "Assemble Stokes system" "Solve Stokes system"
    fi
}

set +a  # return to default setting

if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi