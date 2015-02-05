#base='./2015-02-05_00-12-18_Ez_00006'

base=$1
ncores=$2
if [ -z "$ncores" ]; then
    ncores=1;
fi

width=$(identify -ping -format '%w' ${base}_rank-0_t0.png)
height=$(identify -ping -format '%h' ${base}_rank-0_t0.png)

echo "ncores: $ncores"
echo "$width x $height"

for (( j=0; j < ncores; j++ )); do
    convert -size ${width}x${height} canvas:white ${base}_tmp_${j}.bmp &
done
wait

for (( i=0; i<127; i=i+ncores )); do
    for (( j=0; j < ncores; j++ )); do
        if [ -f "${base}_rank-$((i+j))_t0.png" ]; then
            echo "Merging ${base}_rank-$((i+j))_t0.png onto ${base}_tmp_${j}.bmp"
            convert "${base}_tmp_${j}.bmp" -transparent white "${base}_rank-$((i+j))_t0.png" -transparent white -composite "${base}_tmp_${j}.bmp" &
        fi
    done
    wait
done

convert -size ${width}x${height} canvas:white ${base}_tmp.bmp
for (( j=0; j < ncores; j++ )); do
    echo "Merging ${base}_tmp_${j}.bmp onto ${base}_tmp.bmp"
    convert "${base}_tmp.bmp" -transparent white "${base}_tmp_${j}.bmp" -transparent white -composite "${base}_tmp.bmp"
done
convert ${base}_tmp.bmp ${base}.png

for (( j=0; j < ncores; j++ )); do
    rm ${base}_tmp_${j}.bmp
done

