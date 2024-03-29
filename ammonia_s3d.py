#!/usr/bin/python3

import os
import stat
import sys
import getopt
import re
import tempfile
import shutil
import subprocess

def usage(msg):
    if msg:
        sys.stderr.write(msg + "\n")
    exit(1)

def parse_dim3(s):
    m1 = re.match(r'^[0-9]+$', s)
    if m1:
        return (int(s), int(s), int(s))

    m2 = re.match(r'^([0-9]+)x([0-9]+)x([0-9]+)$', s)
    if m2:
        return (int(m2.group(1)), int(m2.group(2)), int(m2.group(3)))

    usage("Invalid grid size: " + s)

def run_s3d(opts, envvars):
    if "-h" in opts:
        usage()
    basedir = opts.get("-b") or os.environ.get("S3DROOT") or "."
    workdir = opts.get("-w") or os.environ.get("WORK") or "."
    keepdir = opts.get("-k")
    outfile = opts.get("-o") or "-"
    mechanism = opts.get("-m")
    print(opts)
    exe = opts.get("-x") or usage("Must specify executable via -x!")
    lib = opts.get("-l")
    local_grid = parse_dim3(opts.get("-g", "32"))
    proc_grid = parse_dim3(opts.get("-p", "1"))
    timesteps = int(opts.get("-t", "10"))
    spawn = opts.get("-s", "local")
    hosts = opts.get("-H")
    no_lagging = "-f" in opts
    node_file = opts.get("-r")
    vary = parse_dim3(opts.get('-v', "1"))
    title = opts.get("--title", "bomb")
    bomb_file = opts.get("--bomb-file")

    if not mechanism:
        # try to guess mechanism from binary
        with open(exe, "rb") as f:
            data = f.read()
            m = re.search(r'Using chemical mechanism: ([A-Za-z0-9_]+)', str(data))
            if not m:
                print("Not able to guess mechanism yet!")
                exit(1)
            mechanism = m.group(1)
            print("Mechanism appears to be:", mechanism)

    # create a temporary directory for stuff to live in (or use what we're told)
    if keepdir:
        tdir = keepdir
        if not os.path.isdir(keepdir):
            os.mkdir(keepdir)
    else:
        tdir = tempfile.mkdtemp(dir = workdir)
    tdir = os.path.abspath(tdir)

    # need three subdirectories
    idir  = os.path.join(tdir, "input")
    ddir  = os.path.join(tdir, "data")
    rdir  = os.path.join(tdir, "run")
    pdir  = os.path.join(tdir, "post")
    pdirt = os.path.join(tdir, "post/tecplot")

    if not os.path.isdir(idir):
        os.mkdir(idir)
    if not os.path.isdir(ddir):
        os.mkdir(ddir)
    if not os.path.isdir(rdir):
        os.mkdir(rdir)
    if not os.path.isdir(pdir):
        os.mkdir(pdir)
    if not os.path.isdir(pdirt):
        os.mkdir(pdirt)

    # copy over input files
    indep_dir = os.path.join(basedir, "input")
    chemdep_dir = os.path.join(basedir, "input", "chemistry_dependent", mechanism)
    if not os.path.exists(chemdep_dir):
        print("Error: %s does not exist - typo in mechanism?" % chemdep_dir)
        exit(1)

    def copy_files(src, tgt):
        for f in os.listdir(src):
            fullpath = os.path.join(src, f)
            fullpathd = os.path.join(tgt, f)
            if os.path.isfile(fullpath):
    #            shutil.copy(fullpath, tgt)
                shutil.copyfile(fullpath, fullpathd)

    ##copy_files(indep_dir, idir)
    copy_files(chemdep_dir, idir)

    # copy over executable and shared lib, if it exists
    shutil.copy(exe, tdir)
    exe = os.path.join(tdir, os.path.basename(exe))
    if node_file is not None:
        shutil.copy(node_file, tdir)
    #    shutil.copyfile(node_file, tdir)

    if lib:
        # library always gets copied in as librhst.so
        print("copying following file:", lib)
        tlib = os.path.join(tdir, "librhsf.so")
        shutil.copy(lib, tlib)
        print("to tlib:", tlib)
    #    shutil.copyfile(lib, tlib)
        lib = tlib
        #print os.environ["LD_LIBRARY_PATH"]
        os.environ["LD_LIBRARY_PATH"] = tdir + ":" + os.environ["LD_LIBRARY_PATH"]
        print(os.environ["LD_LIBRARY_PATH"])

    handoff = 10
    xmax = 1.0
    delx = xmax/(local_grid[0] * proc_grid[0])
    ymax = (local_grid[1] * proc_grid[1]) * delx
    zmax = (local_grid[2] * proc_grid[2]) * delx
    tstep = .25 * delx * 1e-7
    # now generate the right s3d.in
    with open(os.path.join(idir, "s3d.in"), "w") as f:
        f.write("==========================================================================================\n")
        f.write("MODE\n")
        f.write("==========================================================================================\n")
        f.write("0                   - 0 = run DNS code, 1 = post-process DNS results     (mode)\n")
        f.write("==========================================================================================\n")
        f.write("GRID DIMENSION PARAMETERS\n")
        f.write("==========================================================================================\n")
        f.write("%-6d                 - global number of grid points in the x-direction    (nx_g)\n" % (local_grid[0] * proc_grid[0]))
        f.write("%-6d                 - global number of grid points in the y-direction    (ny_g)\n" % (local_grid[1] * proc_grid[1]))
        f.write("%-6d                 - global number of grid points in the z-direction    (nz_g)\n" % (local_grid[2] * proc_grid[2]))
        f.write("%-6d                 - number of processors in x-direction                (npx)\n" % proc_grid[0])
        f.write("%-6d                 - number of processors in y-direction                (npy)\n" % proc_grid[1])
        f.write("%-6d                 - number of processors in z-direction                (npz)\n" % proc_grid[2])
        f.write("==========================================================================================\n")
        f.write("RUN-TIME PARAMETERS\n")
        f.write("==========================================================================================\n")
        f.write("0                   - 0 for write output to screen, 1 for write to file  (i_write)\n")
        f.write("0                   - 0 for new run, 1 for restart                       (i_restart)\n")
        f.write("%-6d                - ending time step                                   (i_time_end)\n" % timesteps)
        f.write("%d                  - frequency to save fields in restart files          (i_time_save)\n" % handoff)
        f.write("1.0e+5                  - time period to save fields in restart files        (time_save_inc)\n")
        f.write("==========================================================================================\n")
        f.write("GEOMETRY PARAMETERS\n")
        f.write("==========================================================================================\n")
        f.write("%s                  - title of run, sets initialiation of flow field     (run_title)\n" % title)
        f.write("%-6d                - 0 for no x-direction dependence, 1 for so          (vary_in_x)\n" % vary[0])
        f.write("%-6d                - 0 for no y-direction dependence, 1 for so          (vary_in_y)\n" % vary[1])
        f.write("%-6d                - 0 for no z-direction dependence, 1 for so          (vary_in_z)\n" % vary[2])
        f.write("0                   - 0 for non-periodic in x-direction, 1 for periodic  (periodic_x)\n")
        f.write("1                   - 0 for non-periodic in y-direction, 1 for periodic  (periodic_y)\n")
        f.write("1                   - 0 for non-periodic in z-direction, 1 for periodic  (periodic_z)\n")
        f.write("1                   - 0 for stretched edges in x-dir, 1 for uniform      (unif_grid_x)\n")
        f.write("1                   - 0 for stretched edges in y-dir, 1 for uniform      (unif_grid_y)\n")
        f.write("1                   - 0 for stretched edges in z-dir, 1 for uniform      (unif_grid_z)\n")
        f.write("40.0e-6             - minimum grid spacing for streching in x direction  (min_grid_x)\n")
        f.write("40.0e-6             - minimum grid spacing for streching in y direction  (min_grid_y)\n")
        f.write("40.0e-6             - minimum grid spacing for streching in z direction  (min_grid_z)\n")
        f.write("0                   - 0 for no turbulence, 1 for isotropic turbulence    (i_turbulence)\n")
        f.write("-7                  - BC at x=0 boundary; 1 nonreflecting, 0 hard inflow (nrf_x0)\n")
        f.write("-7                  - BC at x=L boundary; 1 nonreflecting, 0 hard inflow (nrf_xl)\n")
        f.write("1                   - BC at y=0 boundary; 1 nonreflecting, 0 hard inflow (nrf_y0)\n")
        f.write("1                   - BC at y=L boundary; 1 nonreflecting, 0 hard inflow (nrf_yl)\n")
        f.write("1                   - BC at z=0 boundary; 1 nonreflecting, 0 hard inflow (nrf_z0)\n")
        f.write("1                   - BC at z=L boundary; 1 nonreflecting, 0 hard inflow (nrf_zl)\n")
        f.write("0.2                 - fix factor for pressure drift                      (relax_ct)\n")
        f.write("1.0                 - fix factor for pressure drift                      (pout)\n")
        f.write("==========================================================================================\n")
        f.write("PHYSICAL PARAMETERS\n")
        f.write("==========================================================================================\n")
        f.write("0.0                 - minimum value of grid in x-direction in cm         (xmin)\n")
        f.write("0.0                 - minimum value of grid in y-direction in cm         (ymin)\n")
        f.write("0.0                 - minimum value of grid in z-direction in cm         (zmin)\n")
        f.write("%-10.6f             - maximum value of grid in x-direction in cm         (xmax)\n" % xmax)
        f.write("%-10.6f             - maximum value of grid in y-direction in cm         (ymax)\n" % ymax)
        f.write("%-10.6f             - minimum value of grid in z-direction in cm         (zmax)\n" % zmax)
        f.write("0.001               - Mach number where re_real/mach_no = re_acoustic    (mach_no)\n")
        f.write("100.0               - real convective Reynolds number                    (re_real)\n")
        f.write("0.708               - Prandtl number                                     (pr)\n")
        f.write("==========================================================================================\n")
        f.write("NUMERICS PARAMETERS\n")
        f.write("==========================================================================================\n")
        f.write("1                   - 0 for no reaction, 1 for reaction                    (i_react)\n")
        f.write("8                   - order of spatial derivatives: 6th or 8th only        (iorder)\n")
        f.write("%d                  - frequency to monitor min/max and active              (i_time_mon)\n" % handoff)
        f.write("-1                  - frequency to check spatial resolution                (i_time_res)\n")
        f.write("-1                 - frequency to write tecplot file                      (i_time_tec)\n")
        f.write("10                  - order of spatial filter                              (i_filter)\n")
        f.write("%d                   - frequency to filter solution vector                  (i_time_fil)\n" %handoff)
        f.write("==========================================================================================\n")
        f.write("REQUIRED REFERENCE VALUES\n")
        f.write("==========================================================================================\n")
        f.write("1.364               - reference ratio of specific heats                    (g_ref)\n")
        f.write("528.0               - reference speed of sound (m/s)                       (a_ref)\n")
        f.write("801.9               - freestream temperature (K)                           (to)\n")
        f.write("12.23               - reference density (kg/m^3)                           (rho_ref)\n")
        f.write("57.69e-3            - reference thermal conductivity (W/m-s)               (lambda_ref) \n")
        f.write("==========================================================================================\n")
        f.write("flag to enable/disable tracer\n")
        f.write("==========================================================================================\n")
        f.write("0                   - tracer control                                       (tracer_ctrl)\n")
        f.write("==========================================================================================\n")
        f.write("flag to enable/disable MPI I/O\n")
        f.write("==========================================================================================\n")
        f.write("1                   - I/O method: 0:Fortran I/O, 1:MPI-IO, 2:PnetCDF, 3:HDF5\n")

    if bomb_file == "0d":
        with open(os.path.join(idir, "bomb.in"), "w") as f:
            f.write("! input parameters for bomb test case \n")
            f.write("30.0        !p_mean\n")
            f.write("800.0      !T_mean \n")
            f.write("0.0         !T_rms \n")
            f.write("1.0         !phi\n")
    else:
        with open(os.path.join(idir, "bomb.in"), "w") as f:
            f.write("! input parameters for bomb test case \n")
            f.write("10.0        !p_mean\n")
            f.write("1250.0      !T_mean \n")
            f.write("300.0        !T_rms \n")
            f.write("1.0         !phi\n")

    with open(os.path.join(idir, "erk.in"), "w") as f:
        f.write("-64         explicit Runge-Kutta (ERK) time integration method          (rk_method)\n")
        f.write("0           0 = no controller (set fixed timestep on next line), 1 = on (cont_switch)\n")
        f.write("%.0E        initial timestep (sec)                                      (tstep_init)\n" % tstep)
        f.write("%.0E        minimum timestep (sec)                                      (tstep_min)\n" % tstep)
        f.write("%.0E        maximum timestep (sec)                                      (tstep_max)\n" % tstep)
        f.write("1.0e-03     relative Runge-Kutta error tolerance                        (rk_rtol)\n")
        f.write("1.0e-15     absolute Runge-Kutta error tolerance                        (rk_atol)\n")
        f.write("0.90        controller saftey factor                                    (cont_safety)\n")
        f.write("0.20        controller integral gain                                    (k_I)\n")
        f.write("0.20        controller proportional gain                                (k_P)\n")
        f.write("0.05        controller derivative gain                                  (k_D)\n")
        f.write("0.00        controller second derivative gain                           (k_D2)\n")
        f.write("3.0         overal maximum change in timestep                           (cont_factor)\n")
        f.write("1           acoustic cfl check, 0=off, 1=on                             (cfl_switch)\n")
        f.write("20          timestep frequency to write controller info to ts.dat file  (i_time_cont)\n")
        f.write("0           1 for johnmc, 0 for original\n\n\n\n\n")
        f.write("Notes About Settings\n")
        f.write("--------------------\n")
        f.write(" 1. EKR scheme should be set to -64.\n")
        f.write(" 2. Initial timestep should be set 10 to 100 smaller than the anticipated CFL limit.\n")
        f.write(" 3. Minimum timestep should be used with CAUTION! It should be set at least 100 times\n")
        f.write("    less than the anticipated CFL limit.\n")
        f.write(" 4. Relative Runge-Kutta error tolerance should be set between 1.0e-3 and 1.0e-4.\n")
        f.write(" 5. Absolute Runge-Kutta error tolerance should be set between 1.0e-12 and 1.0e-14.\n")
        f.write("    Keep this value very small, near machine precision!\n")
        f.write(" 6. Controller saftey factor should be set between 0.8 and 0.9.\n")
        f.write(" 7. Controller gains should always be set to 0.2, 0.2, 0.05, and 0.0, respectively.\n")
        f.write(" 8. Overall maximum change in timestep should be set between 3.0 and 5.0.\n")


    with open(os.path.join(idir, "grid.in"), "w") as f:
        f.write("!  input file for grid stretching module\n")
        f.write("0.00031        ! lfine_x (m)\n")
        f.write("0.000        ! lfine_y (m)\n")
        f.write("0.000        ! lfine_z (m)\n")
        f.write("0            ! imap_x \n")
        f.write("0            ! imap_y \n")
        f.write("0            ! imap_z \n")
        f.write("1.0          ! b_x \n")
        f.write("1.0          ! b_y \n")
        f.write("1.0          ! b_z \n")

    with open(os.path.join(idir, "mixavg.in"), "w") as f:
        f.write("\n")
        f.write("F           F = no baro-diffusion,     T = on        (baro_switch)\n")
        f.write("F           F = no thermal diffusion,  T = on        (thermDiff_switch)\n")
        f.write("\n")
        f.write("T           F = no lagging of coefficients,  T = on  (lagging_switch)\n")
        f.write("1           number of steps lagged                   (lag_steps)\n")
        f.write("1           0 for default diffflux calc, other for johnmc routine\n")
        f.write("1           0 for default heatflux calc, other for johnmc routine\n")
        f.write("1           0 for default mcavis calc, other for johnmc routine\n")

    with open(os.path.join(idir, "streams.in"), "w") as f:
        f.write("MOLE\n")
        f.write("FUEL\n")
        f.write("NH3 0.40\n")
        f.write("N2  0.15\n")
        f.write("END\n")
        f.write("/\n")
        f.write("OXID\n")
        f.write("O2       0.21\n")
        f.write("N2       0.79\n")
        f.write("END\n")

    # now we can actually run it!
    n_ranks = proc_grid[0] * proc_grid[1] * proc_grid[2]
    if spawn == "local":
        assert n_ranks == 1
        if outfile == "-":
            outf = sys.stdout
        else:
            outf = open(outfile, "w")
        #print exe, rdir
        # apply environment variables locally
        for e in envvars:
            v = e.split('=', 1)
            if len(v) == 1:
                os.environ[v[0]] = '1'
            else:
                os.environ[v[0]] = v[1]

        subprocess.call([ exe ], shell = False, cwd = rdir, stdin = None, stdout = outf, stderr = outf)
        if not keepdir:
            # nuke temp tree
            shutil.rmtree(tdir)

    elif spawn == "mpi":
        cmdline = [ "mpirun", "-n", str(n_ranks)]
        if hosts:
            cmdline.extend([ "-H", hosts ])

        cmdline.append("--report-bindings")

        if lib:
            cmdline.extend([ "-x", "LD_LIBRARY_PATH", "--pernode", "--bind-to", "none"])
        else:
            cmdline.extend([ "--bind-to-core", "--bysocket" ])

        for e in envvars:
            cmdline.extend([ "-x", e ])

        #if outfile != "-":
        #    cmdline.extend([ "--output-filename", os.path.abspath(outfile) ])

        # cmdline.extend(["valgrind", "--leak-check=full", "--show-reachable=yes", "--error-limit=no",
        #                 # "--suppressions=/scratch2/seshu/legion_s3d_nscbc/Ammonia_Cases/suppression/cond.supp",
        #                 "--log-file=raw.log", "--gen-suppressions=all"])
        # cmdline.extend(["valgrind", "--leak-check=full", "--show-reachable=yes", "--error-limit=no",
        #                 "--log-file=raw.log"])

        # cmdline.extend(["vtune", "-collect", "hotspots", "-r", "r"])
        # cmdline.extend(["nsys", "profile"])

        cmdline.append(exe)

        if outfile == "-":
            outf = sys.stdout
        else:
            outf = open(outfile, "w")
        print(cmdline)
        subprocess.call(cmdline, cwd = rdir, stdin = None, stdout = outf, stderr = outf)
        if not keepdir:
            # nuke temp tree
            shutil.rmtree(tdir)

    elif spawn == "slurm":
        cmdline = ["srun", "-n", str(n_ranks), "--ntasks-per-node", "4", "--gpus-per-task", "1", "--cpus-per-task", "7"]
        if hosts:
            cmdline.extend([ "-H", hosts ])

        env = "--export=ALL,%s" % ','.join(envvars)
        env += f",LD_LIBRARY_PATH={os.environ['LD_LIBRARY_PATH']}"
        cmdline.append(env)

        # launch = os.path.join(rdir, 'launch.sh')
        # with open(launch, 'w') as f:
        #     f.write("#!/bin/bash\n")
        #     # f.write("if [[ $(($SLURM_PROCID % 4)) == 0 ]]; then\n")
        #     # f.write("\texport LD_PRELOAD=/global/homes/s/seshu/mem_trace/mem_trace.so\n")
        #     # f.write("fi\n")
        #     f.write("$@\n")
        # st = os.stat(launch)
        # os.chmod(launch, st.st_mode | stat.S_IEXEC)
        # cmdline.append("launch.sh")

        cmdline.append(exe)
        # cmdline.extend(["launch.sh", "vtune", "-collect", "hotspots", "-r", "r", exe])

        if outfile == "-":
            outf = sys.stdout
        else:
            outf = open(outfile, "w")
        print(cmdline)
        subprocess.call(cmdline, cwd = rdir, stdin = None, stdout = outf, stderr = outf)
        if not keepdir:
            # nuke temp tree
            shutil.rmtree(tdir)

    elif spawn == "jsrun":
        n_cores = int(os.environ["LSB_MAX_NUM_PROCESSORS"]) - 1
        assert(n_cores % 42 == 0)
        n_nodes = n_cores/42
        print("n_ranks:", n_ranks, "n_cores:", n_cores, "n_nodes:", n_nodes)
        if n_ranks == n_nodes * 1:
            cmdline = [ "jsrun", "-n", str(n_ranks), "--rs_per_host", "1", "--tasks_per_rs", "1", "--cpu_per_rs", "42", "--gpu_per_rs", "6", "--bind", "rs" ]
        elif n_ranks == n_nodes * 2:
            cmdline = [ "jsrun", "-n", str(n_ranks), "--rs_per_host", "2", "--tasks_per_rs", "1", "--cpu_per_rs", "21", "--gpu_per_rs", "3", "--bind", "rs" ]
        elif n_ranks == n_nodes * 4:
            cmdline = [ "jsrun", "-n", str(n_ranks), "--rs_per_host", "4", "--tasks_per_rs", "1", "--cpu_per_rs", "10", "--gpu_per_rs", "1", "--bind", "rs" ]
        elif n_ranks == n_nodes * 5:
            cmdline = [ "jsrun", "-n", str(n_ranks), "--rs_per_host", "5", "--tasks_per_rs", "1", "--cpu_per_rs", "8", "--gpu_per_rs", "1", "--bind", "rs" ]
        elif n_ranks == n_nodes * 6:
            cmdline = [ "jsrun", "-n", str(n_ranks), "--rs_per_host", "6", "--tasks_per_rs", "1", "--cpu_per_rs", "7", "--gpu_per_rs", "1", "--bind", "rs" ]
        else:
            raise Exception('ranks do not evenly divide nodes')

        if node_file != None:
            assert False # not supported by jsrun

        # Workaround for crash caused by libpami_cudahook.so CUDA hijack
        cmdline.append("--smpiargs=-x PAMI_DISABLE_CUDA_HOOK=1 -disable_gpu_hooks")

        for e in envvars:
            assert False # environment variales are already passed by default

        # shim_src = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pick_hcas.sh")
        # shim_tmp = os.path.join(tdir, "pick_hcas.sh")
        # shutil.copy(shim_src, shim_tmp)

        # cmdline.append(shim_tmp)

        # if launcher != '':
        #     cmdline.extend(launcher.split(' '))
        cmdline.append(exe)

        if outfile == "-":
            outf = sys.stdout
        else:
            outf = open(outfile, "w")

        print(cmdline)
        rv = subprocess.call(cmdline, cwd = rdir, stdin = None, stdout = outf, stderr= outf)
        # if checksums and not rv:
        #     print_checksums(os.path.join(tdir, 'data'))
        if not keepdir:
            shutil.rmtree(tdir)
        exit(rv)

    elif spawn == "aprun":
        #cmdline = [ "aprun", "-n", str(n_ranks), "-N", "1", "-d", "16", "-cc", "none" ]
        cmdline = [ "aprun", "-n", str(n_ranks), "-N", "1", "-m", "16gh", "-cc", "none" ]
        if node_file is not None:
            node_file = os.path.join(tdir, node_file)
            cmdline.extend([ "-l", node_file ])
        #if lib:
            #print os.environ["LD_LIBRARY_PATH"]
            #cmdline.extend([ "-e", "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" ])

        #tlib = os.path.join(tdir, "libgasnet-gemini-par.so")
        #shutil.copy("/ccs/home/mebauer/gasnet/lib/libgasnet-gemini-par.so", tlib)

        for e in envvars:
            cmdline.extend([ "-e", e ])

        cmdline.append(exe)

        if outfile == "-":
            outf = sys.stdout
        else:
            outf = open(outfile, "w")

        print(cmdline)
        subprocess.call(cmdline, cwd = rdir, stdin = None, stdout = outf, stderr= outf)
        if not keepdir:
            shutil.rmtree(tdir)

    else:
        assert False

def build_local_opts():
    return {"-x": "fortran/dme/s3d.x",
            "-l": "build/dme/librhst.so",
            "-p": "1",
            "-s": "local"}

def setup_environ():
    cwd = os.getcwd()
    legion_lib = os.path.join(cwd, "legion/language/build/lib")
    kernels = os.path.join(cwd, "kernels/dme_kernels")
    assert(os.path.exists(legion_lib))
    assert(os.path.exists(kernels))
    os.environ['LD_LIBRARY_PATH'] = ':'.join([os.environ['LD_LIBRARY_PATH'], legion_lib, kernels])

def setup_cpu_envvars():
    return ["RHSF_MAPMODE=allcpu", "REALM_FREEZE_ON_ERROR=1",
            "RHST_ARGS='-ll:cpu 1 -ll:csize 16384 -ll:rsize 1024 -ll:gsize 0 -hl:safe_mapper 1 -ll:stacksize 16 -ll:cpu 8"]

def run_local_cpu_test_case(keepdir, testcase_options):
    setup_environ()
    envvars = setup_cpu_envvars()
    basedir = os.path.join(os.getcwd(), "s3d")
    assert(os.path.exists(basedir))

    keepdir = os.path.join(os.getcwd(), keepdir)
    outputfile = os.path.join(keepdir, "dme_mixed_20_energy.txt")
    if os.path.exists(keepdir):
        shutil.rmtree(keepdir)

    opts = build_local_opts()
    opts.update(testcase_options)
    opts['-b'] = basedir
    opts['-k'] = keepdir
    opts['-o'] = outputfile

    run_s3d(opts, envvars)

def run_0d_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_0d_local"
    opts = {"-g": "16",
            "-t": "20000",
            "--bomb-file": "0d"}

    run_local_cpu_test_case(keepdir, opts)

def run_1d_x_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_1d_x_local"
    opts = {"-g": "256x16x16",
            "-v": "1x0x0",
            "-t": "3000"}

    run_local_cpu_test_case(keepdir, opts)

def run_1d_y_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_1d_y_local"
    opts = {"-g": "16x256x16",
            "-v": "0x1x0",
            "-t": "3000"}

    run_local_cpu_test_case(keepdir, opts)

def run_1d_z_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_1d_z_local"
    opts = {"-g": "16x16x256",
            "-v": "0x0x1",
            "-t": "3000"}

    run_local_cpu_test_case(keepdir, opts)

def run_2d_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_2d_local"
    opts = {"-g": "128x128x16",
            "--title": "bomb_test_2d_xy",
            "-t": "2000"}

    run_local_cpu_test_case(keepdir, opts)

def run_3d_local_cpu():
    keepdir = "tmp_dme_mixed_20_energy_3d_local"
    opts = {"-g": "64",
            "--title": "bomb_test",
            "-t": "2000"}

    run_local_cpu_test_case(keepdir, opts)

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "w:o:m:x:l:g:p:b:t:k:s:H:e:v:fhr:", ["title=", "bomb-file="])
    envvars = [ x[1] for x in opts if x[0] == '-e' ]
    opts = dict(opts)
    run_s3d(opts, envvars)
