FILENAME=$0

RELDIR=${FILENAME%set.bash}

## rsyncit="rsync -Crzut"

if [[ -z $1 ]]; then
    echo "You must specify \"Dynamic\" or \"Static\"."
    exit
fi

#name=$1
name=`echo $1 | awk '{print tolower($0)}'`
#echo $name

if [[ $name == dynamic ]]; then
name=Dynamic
fi
if [[ $name == static ]]; then
name=Static
fi

set -x

# Copy files
cp ${RELDIR}${name}/DESCRIPTION ${RELDIR}../DESCRIPTION
cp ${RELDIR}${name}/NAMESPACE   ${RELDIR}../NAMESPACE
cp ${RELDIR}${name}/Makevars    ${RELDIR}../src/Makevars

set +x

# Remove -Wall if nowall is specified as second argument.
if [[ $2 == nowall ]]; then
sed -i~ -e 's/\-Wall//g' ${RELDIR}../src/Makevars
fi

FILES="man/CUBS.Rd man/AR1.Rd man/AR1Ind.Rd R/CUBS.R"

## Add files for dynamic.
if [[ $name == Dynamic ]]; then
  name=Dynamic  
  for FILE in $FILES; do
	set -x
	cp "${RELDIR}${name}/$FILE" "${RELDIR}../$FILE"
	set +x
  done
fi

## Remove files for static.
if [[ $name == Static ]]; then
  name=Static
  for FILE in $FILES; do
	PTF="${RELDIR}../$FILE"
	if [[ -e "$PTF" ]]; then
	  set -x
	  rm "$PTF"
	  set +x
    fi
  done	
fi
