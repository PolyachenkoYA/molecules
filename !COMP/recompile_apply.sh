echo "{---------------- gen begin --------------------"
cd gen/Release
make clean
make all
cd ../..
cp -v ./gen/Release/gen ./!go
echo "---------------- gen done --------------------}"

echo "{---------------- comp begin --------------------"
cd comp/Release
make clean
make all
cd ../..
cp -v ./comp/Release/comp ./!go
echo "---------------- comp done --------------------}"

echo "{---------------- post_proc begin --------------------"
cd post_proc/Release
make clean
make all
cd ../..
cp -v ./post_proc/Release/post_proc ./!go
cp -v ./post_proc/Release/post_proc ../RES
echo "---------------- post_proc done --------------------}"
