cd ~/build/osx64/tools

if [ ! -f qtbase-everywhere-src-5.11.2.tar.xz ]; then
  wget http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/qtbase-everywhere-src-5.11.2.tar.xz
fi

tar xf qtbase-everywhere-src-5.11.2.tar.xz
cd qtbase-everywhere-src-5.11.2

./configure -prefix ../build/Qt-dynamic -opensource -confirm-license -nomake tools -nomake examples -nomake tests

make -j4
