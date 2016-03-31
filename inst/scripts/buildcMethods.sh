file=$1

functions=$(cat $file | grep -P ^void)

for i in 
do
echo $i | awk -F " " '{print $2}' | awk -F "(" '{print $1}'
done

cat $file | grep -P ^void | awk -F " " '{print $2}' | awk -F "(" '{print $1}'

tform{
    i=$1
    if [ $i -eq int ]
    then
	echo "INTSXP"
    elif [ $i == char ]
    then
	echo "STRSXP"
    elif [ $i -eq double ]
    then
	echo "REALSXP"
    fi
}

for i in $(cat $file | grep -P ^void | awk -F "(" '{print $2}' | awk -F ")" '{print $1}' | sed 's/\ /-/g')
do
echo $(echo $i | sed "s/,/\n/g" | awk -F "-" '{print $1}')
done
