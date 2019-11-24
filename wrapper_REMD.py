#!/usr/bin/env python
# -*- coding: utf-8 -*-

# python /home/REMD/scripts/launch_REMD/main.py -seq /home/REMD/data/seq-cyclic_peptide_baker.txt -cyclic True -time 200 -temperature "300 318 337.97 358.81 380.85 404.27 429.12 455.50"

def get_args():
    parser = argparse.ArgumentParser(description = 'REMD wrapper')
    parser.add_argument('-seq', dest='seqfile', type=str, help="the file that contains the sequence(s)", required=True)
    parser.add_argument('-time', dest='time', type=str, help="time of the simulation")
    parser.add_argument('-temperature', dest='temperature', type=int, help="the number of parallel tasks (default: the number of UniProt in the input)")
    parser.add_argument('-c', dest="rpbs_cluster", action='store_true', help="run on RPBS dev or prod cluster")
    parser.add_argument('-cyclic', dest="cyclic", action='store_true', help="activate the cyclization")

    args = parser.parse_args()

    return args.seqfile, args.time, args.temperature, args.rpbs_cluster, args.cyclic



def main():

    seqfile, time, temperature, rpbs_cluster, cyclic = get_args()

    if (rpbs_cluster):
        sys.path.append("/service/env")
        import cluster

        cmd = "python /home/REMD/scripts/launch_REMD/main.py -cyclic True -time 200 -temperature '300 318 337.97 358.81 380.85 404.27 429.12 455.50' -seq "
        args = [ '$(sed -n "${TASK_ID}p" ALL_SUBOUT)',]
        cluster.progress("Start REMD")
        cluster.runTasks(cmd, args, tasks = number_of_extractions, docker_img = "remd_project", job_opts = "-c 8")


if __name__ == '__main__':
    import sys, os, argparse, re
    main()

