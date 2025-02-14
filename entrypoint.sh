#!/bin/bash

# start iris
/iris-main "$@" &

/usr/irissys/dev/Cloud/ICM/waitISC.sh

iop --init

iop --migrate /irisdev/app/src/python/opm/settings.py

iop --default Opm.Production

export LD_LIBRARY_PATH=${ISC_PACKAGE_INSTALLDIR}/bin
iop --start