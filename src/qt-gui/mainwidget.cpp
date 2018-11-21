#include <QtWidgets>
#include "mainwidget.h"

// Constructor for main window
MainWidget::MainWidget(QWidget *parent) :
    QWidget(parent)
{
  // input files
  selectFileButton_ = new QPushButton(tr("Add ms2 files"));
  deleteFileButton_ = new QPushButton(tr("Remove selected files"));
  fileListWidget_ = new QListWidget();
  fileListWidget_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  fileDialog_ = new QFileDialog();
  
  // parameters
  QLabel* pvalThreshLabel = new QLabel(tr("log10(p-value) threshold:"));
  pvalThreshInput_ = new QDoubleSpinBox();
  pvalThreshInput_->setRange(-100.0, 0.0);
  pvalThreshInput_->setValue(-10.0);
  
  QLabel* precTolLabel = new QLabel(tr("Precursor tolerance:"));
  precTolInput_ = new QDoubleSpinBox();
  precTolInput_->setRange(0.0, 1000.0);
  precTolInput_->setValue(20.0);
  precTolUnitInput_ = new QComboBox();
  precTolUnitInput_->insertItem(0,tr("ppm"));
  precTolUnitInput_->insertItem(1,tr("Da"));
  
  QLabel* outputFolderLabel = new QLabel(tr("Output folder:"));
  outputFolderInput_ = new QLineEdit();
  outputFolderDialog_ = new QFileDialog();
  outputFolderButton_ = new QPushButton(tr("Browse..."));
  
  QLabel* otherParamsLabel = new QLabel(tr("Extra command line params:"));
  otherParamsInput_ = new QLineEdit();
  otherParamsInput_->setPlaceholderText(tr("e.g. \"--verbatim 1 --prefix ABC\""));
  
  // run
  runButton_ = new QPushButton(tr("Run MaRaCluster"));
  QLabel* stdoutLabel = new QLabel(tr("Command line output:"));
  textBrowser_ = new QTextBrowser();
  textBrowser_->setLineWrapMode(QTextEdit::NoWrap);
  clearStdoutButton_ = new QPushButton(tr("Clear log"));
  saveStdoutButton_ = new QPushButton(tr("Save log file"));
  
  QLabel* copyrightText = new QLabel(QChar(0x00A9) + tr(" 2018 Matthew The"));
  QLabel* citeText = new QLabel(tr("Please cite: \"MaRaCluster: A Fragment Rarity Metric for Clustering Fragment Spectra"));
  QLabel* linkText = new QLabel(tr("in Shotgun Proteomics\", <a href=\"https://doi.org/10.1021/acs.jproteome.5b00749\">https://doi.org/10.1021/acs.jproteome.5b00749</a>"));
  linkText->setTextFormat(Qt::RichText);
  linkText->setTextInteractionFlags(Qt::TextBrowserInteraction);
  linkText->setOpenExternalLinks(true);
  
  mainLayout_ = new QGridLayout;
  
  // left panel
  mainLayout_->addWidget(selectFileButton_,  0,0);
  mainLayout_->addWidget(deleteFileButton_,  0,1,1,2);
  
  mainLayout_->addWidget(fileListWidget_,    1,0,1,3);
  
  mainLayout_->addWidget(pvalThreshLabel,    2,0);
  mainLayout_->addWidget(pvalThreshInput_,   2,1,1,2);
  
  mainLayout_->addWidget(precTolLabel,       3,0);
  mainLayout_->addWidget(precTolInput_,      3,1);
  mainLayout_->addWidget(precTolUnitInput_,  3,2);
  
  mainLayout_->addWidget(outputFolderLabel,  4,0);
  mainLayout_->addWidget(outputFolderInput_, 4,1);
  mainLayout_->addWidget(outputFolderButton_,4,2);
  
  mainLayout_->addWidget(otherParamsLabel,   5,0);
  mainLayout_->addWidget(otherParamsInput_,  5,1,1,2);
  mainLayout_->addWidget(runButton_,         6,0,1,3);
  
  mainLayout_->addWidget(copyrightText,      7,0,1,3);
  mainLayout_->addWidget(citeText,           7,3,1,3);
  mainLayout_->addWidget(linkText,           8,3,1,3);
  
  // right panel
  mainLayout_->addWidget(stdoutLabel,        0,3);
  mainLayout_->addWidget(saveStdoutButton_,  0,4);
  mainLayout_->addWidget(clearStdoutButton_, 0,5);
  mainLayout_->addWidget(textBrowser_,       1,3,6,3);
  
  setLayout(mainLayout_);
  setWindowTitle(tr("MaRaCluster"));

  connect(selectFileButton_, SIGNAL(released()), this, SLOT(onSelectFilesButtonReleased()));
  connect(deleteFileButton_, SIGNAL(released()), this, SLOT(onDeleteFilesButtonReleased()));
  connect(outputFolderButton_, SIGNAL(released()), this, SLOT(onOutputFolderButtonReleased()));
  connect(runButton_, SIGNAL(released()), this, SLOT(onRunButtonReleased()));
  connect(saveStdoutButton_, SIGNAL(released()), this, SLOT(onSaveStdoutButtonReleased()));
  connect(clearStdoutButton_, SIGNAL(released()), this, SLOT(onClearStdoutButtonReleased()));
  connect(&process_, SIGNAL(readyReadStandardOutput()), this, SLOT(onCaptureProcessOutput()));
  connect(&process_, QOverload<QProcess::ProcessError>::of(&QProcess::errorOccurred), 
    [=](QProcess::ProcessError error){ 
      if (error == QProcess::FailedToStart) {
        textBrowser_->insertPlainText(tr("MaRaCluster execution failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.\n"));
      } else {
        textBrowser_->insertPlainText(tr("MaRaCluster execution failed after starting successfully.\n"));
      }
      textBrowser_->ensureCursorVisible();
      runButton_->setText("Run MaRaCluster");
    });
  connect(&process_, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
    [=](int exitCode, QProcess::ExitStatus exitStatus){ 
      if (exitStatus == QProcess::NormalExit) {
        textBrowser_->insertPlainText(tr("MaRaCluster execution finished successfully\n"));
      } else {
        textBrowser_->insertPlainText(tr("MaRaCluster execution failed (exit code: ") + QString::number(exitCode) + tr(")\n"));
      }
      textBrowser_->ensureCursorVisible();
      runButton_->setText("Run MaRaCluster");
    });
}

// Destructor
MainWidget::~MainWidget()
{
  delete selectFileButton_;
  delete deleteFileButton_;
  delete outputFolderButton_;
  delete saveStdoutButton_;
  delete clearStdoutButton_;
  delete runButton_;
  
  delete pvalThreshInput_;
  delete precTolInput_;
  delete precTolUnitInput_;
  delete otherParamsInput_;
  delete outputFolderInput_;
  
  delete fileDialog_;
  delete outputFolderDialog_;
  delete fileListWidget_;
  delete textBrowser_;
  delete mainLayout_;
}

// Handler for button click
void MainWidget::onSelectFilesButtonReleased()
{  
  fileDialog_->setDirectory(QDir::homePath());
  fileDialog_->setFileMode(QFileDialog::ExistingFiles);
  fileDialog_->setNameFilter(tr("ms2 files *.mzML, *.ms2, *.mgf (*.mzML *.ms2 *.mgf)"));
  
  QStringList fileNames;
  if (fileDialog_->exec())
    fileNames = fileDialog_->selectedFiles();
  
  for (const auto& file : fileNames) {
    new QListWidgetItem(file, fileListWidget_);
  }
}

void MainWidget::onDeleteFilesButtonReleased()
{
  qDeleteAll(fileListWidget_->selectedItems());
}

void MainWidget::onOutputFolderButtonReleased()
{  
  outputFolderDialog_->setDirectory(QDir::homePath());
  outputFolderDialog_->setFileMode(QFileDialog::Directory);
  
  QStringList fileNames;
  if (outputFolderDialog_->exec())
    fileNames = outputFolderDialog_->selectedFiles();
    if (fileNames.size() > 0)
      outputFolderInput_->setText(fileNames.front());
}

void MainWidget::onRunButtonReleased()
{
  if (runButton_->text() == "Stop MaRaCluster") {
    textBrowser_->insertPlainText(tr("Attempting to stop MaRaCluster run.\n"));
    textBrowser_->ensureCursorVisible();
    process_.kill();
  } else {
    if (fileListWidget_->count() == 0) {
      textBrowser_->insertPlainText(tr("No input files detected.\n"));
      textBrowser_->ensureCursorVisible();
      return;
    }
    
    QString outputFolder(outputFolderInput_->text());
    if (outputFolder.size() == 0) {
      textBrowser_->insertPlainText(tr("No output folder set.\n"));
      textBrowser_->ensureCursorVisible();
      return;
    }
    
    QDir outputDir(outputFolder);
    if (!outputDir.exists()) {
      bool success = outputDir.mkpath(outputFolder);
      if (!success) {
        textBrowser_->insertPlainText(tr("Could not create output directory: ") + outputFolder + tr("\n"));
        textBrowser_->ensureCursorVisible();
        return;
      }
    }
    
    QFile file(outputFolder + QDir::separator() + "file_list.txt");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      textBrowser_->insertPlainText(tr("Could not create output directory: ") + outputFolder + tr("\n"));
      textBrowser_->ensureCursorVisible();
      return;
    }

    QTextStream out(&file);
    for (int i = 0; i < fileListWidget_->count(); ++i)
    {
      QListWidgetItem* item = fileListWidget_->item(i);
      out << item->text() << endl;
    }
    file.close();
    
    
    process_.setProcessChannelMode(QProcess::MergedChannels);
    
    // Set up our process to write to stdout and run our command
    process_.setCurrentWriteChannel(QProcess::StandardOutput); // Set the write channel
    QString cmd;
    cmd += "\"" + QCoreApplication::applicationDirPath() + "/maracluster\" batch --batch " + file.fileName();
    cmd += " --output-folder " + outputFolder;
    cmd += " --pvalThreshold " + QString::number(pvalThreshInput_->value());
    cmd += " --precursorTolerance " + QString::number(precTolInput_->value());
    cmd += " " + otherParamsInput_->text();
    
    textBrowser_->insertPlainText(tr("Command: ") + cmd + '\n');
    textBrowser_->ensureCursorVisible();
    
    runButton_->setText("Stop MaRaCluster");
    
    process_.start(cmd); // Start a terminal command
  }
}

void MainWidget::onClearStdoutButtonReleased()
{
  textBrowser_->clear();
}

void MainWidget::onSaveStdoutButtonReleased()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                            "maracluster.log",
                            tr("Log file (*.txt *.log)"));
  if (fileName.size() > 0) {
    QFile file(fileName);
    
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      textBrowser_->insertPlainText(tr("Could not write to selected file\n."));
      textBrowser_->ensureCursorVisible();
      return;
    }
    
    QTextStream out(&file);
    for (int i = 0; i < fileListWidget_->count(); ++i)
    {
      out << textBrowser_->toPlainText() << endl;
    }
    file.close();
  }
}

// This is called whenever the QProcess::readyReadStandardOutput() signal is received
void MainWidget::onCaptureProcessOutput()
{
  // Determine whether the object that sent the signal was a pointer to a process
  QProcess* process = qobject_cast<QProcess*>(sender());
  // If so, append the output to the textbrowser
  if (process) {
    textBrowser_->insertPlainText(process->readAllStandardOutput());
    textBrowser_->ensureCursorVisible();
  }
}
