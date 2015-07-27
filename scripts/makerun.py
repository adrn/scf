# coding: utf-8

""" Make an SCF run given initial conditions, timestep, and potential parameters """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os
import sys
import shutil

# Third-party
from astropy.constants import G
from astropy import log as logger
import astropy.units as u

# Project
from scf import project_path
import scf.potential as sp

def main(name, pos, vel, scfpars, potentials, run_path, overwrite=False, submit=False):
    run_path = os.path.abspath(run_path)
    # run_path = os.path.abspath(os.path.join(project_path, "simulations", "runs"))
    template_path = os.path.abspath(os.path.join(project_path, "templates"))
    src_path = os.path.abspath(os.path.join(project_path, "src"))
    logger.debug("Run path: {}".format(run_path))

    path = os.path.join(run_path, name)
    if os.path.exists(path) and overwrite:
        logger.info("nuking directory...")
        shutil.rmtree(path)
    elif os.path.exists(path) and not overwrite:
        logger.warning("Path '{}' already exists -- do you want to overwrite?"
                       .format(path))
        yn = raw_input("[y/N]: ")
        if yn.lower().strip() == 'y':
            logger.info("nuking directory...")
            shutil.rmtree(path)
        else:
            logger.info("aborting...")
            sys.exit(0)

        raise IOError("Path '{}' already exists!".format(path))
        sys.exit(0)

    if not os.path.exists(path):
        os.makedirs(path)

    # parse pos and vel
    x_vals,x_unit = pos.split()
    v_vals,v_unit = vel.split()

    pos = map(float, x_vals.split(',')) * u.Unit(x_unit)
    pos = pos.to(u.kpc).value

    vel = map(float, v_vals.split(',')) * u.Unit(v_unit)
    vel = vel.to(u.km/u.s).value

    # combine all of the different blocks
    allblocks = dict()
    for potential_name in potentials:
        blocks = getattr(sp,potential_name)()
        for k,v in blocks.items():
            if k not in allblocks:
                allblocks[k] = ""

            allblocks[k] += v

            if k != "SCFPOT":
                allblocks[k] += "\n"

    # ------------------------------------------------------------------------
    # SCFPAR
    scfbipath = os.path.join(src_path, "SCFBI")
    with open(os.path.join(template_path, "SCFPAR.tpl"), 'r') as f:
        base_SCFPAR = f.read()

    with open(os.path.join(path, "SCFPAR"), 'w') as f:
        f.write(base_SCFPAR.format(x=pos, v=vel, scfbipath=scfbipath, **scfpars))

    # ------------------------------------------------------------------------
    # Makefile
    with open(os.path.join(template_path, "Makefile.tpl"), 'r') as f:
        base_Makefile = f.read()

    with open(os.path.join(path, "Makefile"), 'w') as f:
        f.write(base_Makefile)

    # cmd = "unexpand {fn} > {fn}2".format(fn=os.path.join(path, "Makefile"))
    # os.system(cmd)
    # cmd = "mv {fn}2 {fn}".format(fn=os.path.join(path, "Makefile"))
    # os.system(cmd)

    # ------------------------------------------------------------------------
    # potential.h
    with open(os.path.join(template_path, "potential.h.tpl"), 'r') as f:
        base_potentialh = f.read()

    with open(os.path.join(path, "potential.h"), 'w') as f:
        f.write(base_potentialh.format(hblock=allblocks['hblock']))

    # ------------------------------------------------------------------------
    # potential.f
    with open(os.path.join(template_path, "potential.f.tpl"), 'r') as f:
        base_potentialf = f.read()

    with open(os.path.join(path, "potential.f"), 'w') as f:
        f.write(base_potentialf.format(**allblocks))

    # ------------------------------------------------------------------------
    # SCFPOT
    with open(os.path.join(path, "SCFPOT"), 'w') as f:
        f.write(allblocks['SCFPOT'])

    if submit:
        with open(os.path.join(template_path, "submit.tpl"), 'r') as f:
            base_submit = f.read()

        with open(os.path.join(path, "submit.sh"), 'w') as f:
            f.write(base_submit.format(path=path, name=name))

    # Copy scf.f and scf.h
    shutil.copyfile(os.path.join(src_path, 'scf.f'),
                    os.path.join(path, 'scf.f'))
    shutil.copyfile(os.path.join(src_path, 'scf.h'),
                    os.path.join(path, 'scf.h'))

if __name__ == '__main__':
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Be chatty! (default = False)")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                        default=False, help="Be quiet! (default = False)")
    parser.add_argument("-o", "--overwrite", action="store_true", dest="overwrite",
                        default=False, help="DESTROY. DESTROY.")

    parser.add_argument("--path", "-p", dest="path", type=str, help="Path.", required=True)
    parser.add_argument("--name", dest="name", type=str, help="Name.", required=True)
    parser.add_argument("--potentials", dest="potentials", type=str,
                        help="Potential names. One of: {0}".format(",".join(sp.__all__)),
                        nargs='+', required=True)
    parser.add_argument("-s", "--submit", action="store_true", dest="submit",
                        default=True, help="Generate a submit.sh file as well.")
    parser.add_argument("--pos", dest="x", required=True,
                        type=str, help="Initial position in 3 space. Must specify a "
                                       "unit, e.g., --x='15.324,123.314,51.134 kpc'")
    parser.add_argument("--vel", dest="v", required=True,
                        type=str, help="Initial position in 3 space. Must specify a "
                                       "unit, e.g., --v='15.324,123.314,51.134 km/s'")

    # Parameters for SCF
    parser.add_argument("--dt", dest="dt", default=1., type=float,
                        help="Time step in Myr.")
    parser.add_argument("--nsteps", dest="nsteps", type=int,
                        help="Number of steps.", required=True)
    parser.add_argument("--ncen", dest="ncen", type=int, default=10,
                        help="Output to SCFCEN every (this number) steps.")
    parser.add_argument("--nsnap", dest="nsnap", type=int, default=1000,
                        help="Output to a SNAP file every (this number) steps.")
    parser.add_argument("--ntide", dest="ntide", default=500, type=int,
                        help="Number of steps over which to turn on the tidal field.")
    parser.add_argument("--mass", dest="mass", type=float, required=True,
                        help="Mass in solar masses within 35 scale radii of the "
                             "cluster/satellite.")
    parser.add_argument("--rscale", dest="rscale", default=None, type=float,
                        help="Scale radius of the cluster/satellite.")
    parser.add_argument("--SCFBI", dest="scfbi_path", default=None, type=str,
                        help="Path to an SCFBI file.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    # Contains user-specified parameters for SCF
    scfpars = dict()

    if args.rscale is None:
        ru = 0.43089*(args.mass/2.5e9)**(1/3.)
    else:
        ru = args.rscale
    scfpars['rscale'] = ru
    scfpars['mass'] = args.mass

    _G = G.decompose(bases=[u.kpc,u.M_sun,u.Myr]).value
    X = (_G / ru**3 * args.mass)**-0.5

    length_unit = u.Unit("{0} kpc".format(ru))
    mass_unit = u.Unit("{0} M_sun".format(args.mass))
    time_unit = u.Unit("{:08f} Myr".format(X))

    scfpars['dt'] = args.dt / (1*time_unit).to(u.Myr).value
    scfpars['nsteps'] = args.nsteps
    scfpars['ncen'] = args.ncen
    scfpars['nsnap'] = args.nsnap
    scfpars['ntide'] = args.ntide
    if args.scfbi_path is None:
        scfpars['SCFBIpath'] = os.path.abspath(os.path.join(project_path, 'src', 'SCFBI'))

    main(name=args.name, pos=args.x, vel=args.v, scfpars=scfpars,
         potentials=args.potentials, overwrite=args.overwrite,
         submit=args.submit, run_path=args.path)
