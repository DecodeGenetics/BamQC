#!/bin/bash

#------------------------------------------------------------------
# Install script for bamqualcheck
#------------------------------------------------------------------
# $Source$
# $Revision$ 
# $Author$
# $Date$

set -e                # Exit on errors

[ $# = 1 ] || { echo "You must supply an existing CVS tag"; exit 1; }
CVSTAG="$1"


mod=bamqualcheck
prefix=/nfs_mount/bioinfo/apps-$(uname -i)
DESTDIR=${prefix}/${mod}/${CVSTAG}
BINDIR=${prefix}/bin
MODULE=tools/${mod}
BINARIES="bamqualcheck"

NEEDMAKE=1

[ -d "${DESTDIR}" ] && { echo "Tag '${CVSTAG}' already installed"; exit 1; }

TMPDIR=$(mktemp -d) || { echo "mktemp failed"; exit 1; }

cd ${TMPDIR}

git archive --format=tar --remote git@lsource2.decode.is:/bioinfo/dam.git ${CVSTAG}:${MODULE} > git.tar
tar xvf git.tar > /dev/null
echo ${CVSTAG} > version.txt

[ -z ${NEEDMAKE} ] || make || { echo "make failed"; exit 1; }

if [ "${CVSTAG}" == "HEAD" ]
then

    DESTDIR=${DESTDIR}_${USER}
    if [ -e ${DESTDIR} ]
    then
      rm -rf ${DESTDIR}
    fi

fi

mkdir -p ${DESTDIR}
install -m 755 ${BINARIES} ${DESTDIR}

if [ "${CVSTAG}" == "HEAD" ]
then
    echo "Not deploying version ${CVSTAG}"
    chmod -R u+w ${DESTDIR}
    chmod o-rx ${DESTDIR}
else
    chmod -R a-w ${DESTDIR}
(
    cd ${BINDIR}
    for x in ${BINARIES}
    do
      ln -fs ../${mod}/${CVSTAG}/${x} .
    done
)
fi

cd /
rm -rf ${TMPDIR}

## Done
