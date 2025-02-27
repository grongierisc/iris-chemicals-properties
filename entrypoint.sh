#!/bin/bash
export LD_LIBRARY_PATH=${ISC_PACKAGE_INSTALLDIR}/bin
# start iris
/iris-main "$@" &

/usr/irissys/dev/Cloud/ICM/waitISC.sh

iop --init

iop --migrate /irisdev/app/src/python/opm/settings.py

iop --default Opm.Production

iop --start --detach

# start the web server
python3 /irisdev/app/src/python/opm/app.py