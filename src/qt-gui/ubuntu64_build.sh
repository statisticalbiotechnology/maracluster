sudo apt install patchelf

mkdir -p ~/build/ubuntu64/maracluster-gui
cd ~/build/ubuntu64/maracluster-gui

cmake -DCMAKE_INSTALL_PREFIX=../maracluster-gui-zip -DCMAKE_PREFIX_PATH=~/build/ubuntu64/tools/build/Qt-dynamic/ ~/maracluster-gui/
make -j4
sudo make install

sudo patchelf --set-rpath '$ORIGIN/../../lib/maracluster' ../maracluster-gui-zip/bin/platforms/libqxcb.so

