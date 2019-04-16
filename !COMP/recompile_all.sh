cd gen/Release
make all
cd ../..
echo "---------------- gen done --------------------"

cd comp/Release
make all
cd ../..
echo "---------------- comp done --------------------"

cd post_proc/Release
make all
cd ../..
echo "---------------- post_proc done --------------------"

./apply_recompiled.sh