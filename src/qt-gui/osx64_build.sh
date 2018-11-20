mkdir -p ~/build/osx64/maracluster-gui
cd ~/build/osx64/maracluster-gui

cmake -DCMAKE_INSTALL_PREFIX=../maracluster-gui-zip -DCMAKE_PREFIX_PATH=~/build/osx64/tools/build/Qt-dynamic/ ~/maracluster-gui/
make -j4
make package -j4
#sudo make install
