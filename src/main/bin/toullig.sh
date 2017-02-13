#!/usr/bin/env bash

#
# This script set the right classpath and start the application
#
# Author : Aurélien Birer
#

# Function to create lib paths
make_paths() {

	local RESULT=
	for lib in `ls $1`
	do
		if [ -f $1/$lib ]; then
			RESULT=$RESULT:$1/$lib
		fi
	done

	echo $RESULT
}

# Get the path to this script
REAL_SCRIPT_PATH=`readlink -f $0`
BASEDIR=`dirname $REAL_SCRIPT_PATH`

# Set the Toullig libraries path
LIBDIR=$BASEDIR/lib

# Set the memory in MiB needed by Toullig (only Java part, not external tools)
# By Default 2048
if [ -n "$TOULLIG_MEMORY" ]; then
	MEMORY=$TOULLIG_MEMORY
else
	MEMORY=2048
fi

# Additional JVM options
if [ -n "$TOULLIG_JVM_OPTS" ]; then
	JVM_OPTS=$TOULLIG_JVM_OPTS
else
	JVM_OPTS="-server"
fi

# Add here your plugins and dependencies
if [ -n "$TOULLIG_PLUGINS" ]; then
	PLUGINS=$TOULLIG_PLUGINS
else
	PLUGINS=
fi

# Parse options
OPTERR=0
while getopts “j:m:J:p:” OPTION
do
	case $OPTION in
		j)
			JAVA_HOME=$OPTARG
		;;
		m)
			MEMORY=$OPTARG
		;;
		J)
			JVM_OPTS=$OPTARG
		;;
		p)
			PLUGINS=$OPTARG
		;;
		w)
			cd $OPTARG
		;;
	esac
done

# Set the path to java
if [ -n "$TOULLIG_JAVA_HOME" ] ; then
	JAVA_CMD="$TOULLIG_JAVA_HOME/bin/java"
elif [ -n "$JAVA_HOME" ] ; then
	JAVA_CMD="$JAVA_HOME/bin/java"
else
	JAVA_CMD="java"
fi

# Set the temporary directory
TMPDIR_JVM_OPT=""
if [ -n "$TMPDIR" ]; then
	TMPDIR_JVM_OPT="-Djava.io.tmpdir=$TMPDIR"
fi

COMMON_LIBS=$(make_paths $LIBDIR)
PLUGINS_LIBS=$(make_paths $TOULLIG_PLUGINS)
APP_CLASSPATH=$COMMON_LIBS:$PLUGINS:$PLUGINS_LIBS

# Launch Toullig
$JAVA_CMD \
		$TMPDIR_JVM_OPT \
		$JVM_OPTS \
		-Xmx${MEMORY}m \
		-cp $APP_CLASSPATH \
		-Deoulsan.script.path="$0" \
		-Deoulsan.classpath=$APP_CLASSPATH \
		-Deoulsan.memory=$MEMORY \
		-Deoulsan.launch.mode=local \
		-Deoulsan.launch.script.path=$0 \
		-Deoulsan.hadoop.libs=$COMMON_LIBS:$PLUGINS:$PLUGINS_LIB \
		fr.ens.biologie.genomique.toullig.Main "$@"
