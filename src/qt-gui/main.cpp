/******************************************************************************  
  Copyright 2019 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
 ******************************************************************************/

#include <QtWidgets>
#include "mainwidget.h"

int main(int argc, char *argv[])
{
    // Creates an instance of QApplication
    QApplication a(argc, argv);

    // This is our MainWidget class containing our GUI and functionality
    MainWidget w;
    w.show(); // Show main window

    // run the application and return execs() return value/code
    return a.exec();
}
