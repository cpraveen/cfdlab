#!/bin/sh
# This requires network access

if [ ! -d ".svn" ] ; then
  echo "ERROR: Not in a Subversion working copy" 1>&2
  exit 1
fi

LOGFILE=./ChangeLog
TRUNK=`svn info | grep "URL" | tail -c 6`
if [ ${TRUNK} = "trunk" ] ; then
  URL=`svn info | grep "Repository Root" | cut -c18-`
  echo "Generating ChangeLog for ${URL} (trunk, branches)"
  svn log -v ${URL} trunk branches > ${LOGFILE}
else
  URL=`svn info | grep "URL" | cut -c6-`
  echo "Generating ChangeLog for ${URL}"
  svn log -v .@HEAD > ${LOGFILE}
fi
