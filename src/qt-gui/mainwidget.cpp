#include <QtWidgets>
#include "mainwidget.h"

// Constructor for main window
MainWidget::MainWidget(QWidget *parent) :
    QMainWindow(parent)
{
  // input files
  QPushButton* selectFileButton = new QPushButton(tr("Add ms2 files"), this);
  QPushButton* deleteFileButton = new QPushButton(tr("Remove selected files"), this);
  fileListWidget_ = new QListWidget(this);
  fileListWidget_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  fileDialog_ = new QFileDialog(this);
  
  // parameters
  QLabel* pvalThreshLabel = new QLabel(tr("log10(p-value) threshold:"), this);
  pvalThreshInput_ = new QDoubleSpinBox(this);
  pvalThreshInput_->setRange(-100.0, 0.0);
  pvalThreshInput_->setValue(-10.0);
  
  QLabel* precTolLabel = new QLabel(tr("Precursor tolerance:"), this);
  precTolInput_ = new QDoubleSpinBox(this);
  precTolInput_->setRange(0.0, 1000.0);
  precTolInput_->setValue(20.0);
  precTolUnitInput_ = new QComboBox(this);
  precTolUnitInput_->insertItem(0, tr("ppm"));
  precTolUnitInput_->insertItem(1, tr("Da"));
  
  QLabel* outputFolderLabel = new QLabel(tr("Output folder:"), this);
  outputFolderInput_ = new QLineEdit(this);
  outputFolderDialog_ = new QFileDialog(this);
  QPushButton* outputFolderButton = new QPushButton(tr("Browse..."), this);
  
  QLabel* consensusFileLabel = new QLabel(tr("Consensus spectra file (optional):"), this);
  consensusFileInput_ = new QLineEdit(this);
  consensusFileDialog_ = new QFileDialog(this);
  QPushButton* consensusFileButton = new QPushButton(tr("Browse..."), this);
  
  QLabel* otherParamsLabel = new QLabel(tr("Extra command line params:"), this);
  otherParamsInput_ = new QLineEdit(this);
  otherParamsInput_->setPlaceholderText(tr("e.g. \"--verbatim 1 --prefix ABC\""));
  
  // run
  runButton_ = new QPushButton(tr("Run MaRaCluster"), this);
  
  QLabel* stdoutLabel = new QLabel(tr("Command line output:"), this);
  textBrowser_ = new QTextBrowser(this);
  textBrowser_->setLineWrapMode(QTextEdit::NoWrap);
  QPushButton* clearStdoutButton = new QPushButton(tr("Clear log"), this);
  QPushButton* saveStdoutButton = new QPushButton(tr("Save log file"), this);
  
  QLabel* copyrightText = new QLabel(QChar(0x00A9) + tr(" 2018 Matthew The"), this);
  QLabel* citeText = new QLabel(tr("Please cite: \"MaRaCluster: A Fragment Rarity Metric for Clustering Fragment Spectra"), this);
  QLabel* linkText = new QLabel(tr("in Shotgun Proteomics\", <a href=\"https://doi.org/10.1021/acs.jproteome.5b00749\">https://doi.org/10.1021/acs.jproteome.5b00749</a>"), this);
  linkText->setTextFormat(Qt::RichText);
  linkText->setTextInteractionFlags(Qt::TextBrowserInteraction);
  linkText->setOpenExternalLinks(true);
  
  // Menu Bar
  QMenuBar* menuBar = new QMenuBar;

  QMenu* fileMenu = menuBar->addMenu(tr("&File"));
  QMenu* helpMenu = menuBar->addMenu(tr("&Help"));
  
  QAction* loadAct = new QAction(tr("&Open project"), this);
  loadAct->setShortcuts(QKeySequence::Open);
  connect(loadAct, &QAction::triggered, this, &MainWidget::loadProject);
  
  QAction* saveAct = new QAction(tr("&Save project"), this);
  saveAct->setShortcuts(QKeySequence::Save);
  connect(saveAct, &QAction::triggered, this, &MainWidget::saveProject);
  
  QAction* exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcuts(QKeySequence::Quit);
  connect(exitAct, &QAction::triggered, this, &QWidget::close);
  
  fileMenu->addAction(loadAct);
  fileMenu->addAction(saveAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);
  
  QAction* helpAct = new QAction(tr("&Help"), this);
  connect(helpAct, &QAction::triggered, this, &MainWidget::help);
  
  helpMenu->addAction(helpAct);
  
  setMenuBar(menuBar);
  
  // Widgets
  mainLayout_ = new QGridLayout;
  
  // left panel
  mainLayout_->addWidget(selectFileButton,      0,0);
  mainLayout_->addWidget(deleteFileButton,      0,1,1,2);
  
  mainLayout_->addWidget(fileListWidget_,       1,0,1,3);
  
  mainLayout_->addWidget(pvalThreshLabel,       2,0);
  mainLayout_->addWidget(pvalThreshInput_,      2,1,1,2);
  
  mainLayout_->addWidget(precTolLabel,          3,0);
  mainLayout_->addWidget(precTolInput_,         3,1);
  mainLayout_->addWidget(precTolUnitInput_,     3,2);
  
  mainLayout_->addWidget(outputFolderLabel,     4,0);
  mainLayout_->addWidget(outputFolderInput_,    4,1);
  mainLayout_->addWidget(outputFolderButton,    4,2);
  
  mainLayout_->addWidget(consensusFileLabel,    5,0);
  mainLayout_->addWidget(consensusFileInput_,   5,1);
  mainLayout_->addWidget(consensusFileButton,   5,2);
  
  mainLayout_->addWidget(otherParamsLabel,      6,0);
  mainLayout_->addWidget(otherParamsInput_,     6,1,1,2);
  mainLayout_->addWidget(runButton_,            7,0,1,3);
  
  mainLayout_->addWidget(copyrightText,         9,0,1,3);
  
  // right panel
  mainLayout_->addWidget(stdoutLabel,           0,3);
  mainLayout_->addWidget(saveStdoutButton,      0,4);
  mainLayout_->addWidget(clearStdoutButton,     0,5);
  mainLayout_->addWidget(textBrowser_,          1,3,7,3);
  mainLayout_->addWidget(citeText,              8,3,1,3);
  mainLayout_->addWidget(linkText,              9,3,1,3);
  
  setCentralWidget(new QWidget);
  centralWidget()->setLayout(mainLayout_);
  
  setWindowTitle(tr("MaRaCluster"));

  connect(selectFileButton, SIGNAL(released()), this, SLOT(onSelectFilesButtonReleased()));
  connect(deleteFileButton, SIGNAL(released()), this, SLOT(onDeleteFilesButtonReleased()));
  connect(outputFolderButton, SIGNAL(released()), this, SLOT(onOutputFolderButtonReleased()));
  connect(consensusFileButton, SIGNAL(released()), this, SLOT(onConsensusFileButtonReleased()));
  connect(runButton_, SIGNAL(released()), this, SLOT(onRunButtonReleased()));
  
  connect(saveStdoutButton, SIGNAL(released()), this, SLOT(onSaveStdoutButtonReleased()));
  connect(clearStdoutButton, SIGNAL(released()), this, SLOT(onClearStdoutButtonReleased()));
  
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

// Destruction of QObjects is done by parent-child relationships
MainWidget::~MainWidget() {}

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

void MainWidget::onConsensusFileButtonReleased()
{  
  consensusFileDialog_->setDirectory(QDir::homePath());
  consensusFileDialog_->setFileMode(QFileDialog::AnyFile);
  fileDialog_->setNameFilter(tr("ms2 files *.mzML, *.ms2, *.mgf (*.mzML *.ms2 *.mgf)"));
  
  QStringList fileNames;
  if (consensusFileDialog_->exec())
    fileNames = consensusFileDialog_->selectedFiles();
    if (fileNames.size() > 0)
      consensusFileInput_->setText(fileNames.front());
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
    
    
    QString consensusFileString(consensusFileInput_->text());
    if (!consensusFileString.isEmpty()) {
      QFile consensusFile(consensusFileString);
      if (!consensusFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        textBrowser_->insertPlainText(tr("Could not create consensus file: ") + consensusFileString + tr("\n"));
        textBrowser_->ensureCursorVisible();
        return;
      }
      consensusFile.close();
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
    if (!consensusFileString.isEmpty()) {
      cmd += " --specOut " + consensusFileString;
    }
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

void MainWidget::help() 
{
  QDialog *dialog = new QDialog(this);
  dialog->setWindowTitle(QString("Help"));
  dialog->resize(500, 400);
 
  QTextBrowser *scroll = new QTextBrowser(dialog);
  
  scroll->setOpenExternalLinks(true);
  scroll->setHtml("To cluster fragment spectra, follow these steps: <br/>\
<br/>\
1. Select ms2 files using the \"Add ms2 files\" button; currently, the following formats are supported: mzMl, ms2, mgf. To remove previously added files, simply select the files (using Ctrl/Shift to select multiple files) and use the \"Remove selected files\" button.<br/>\
2. Select a log10(p-val) threshold. The default of -10.0 works for most use cases; larger clusters can be obtained by increasing this value, e.g. -5.0, smaller clusters by lowering this values, e.g. -15.0.<br/>\
3. Select a precursor tolerance. The default of 20ppm works for high-resolution MS1 data. For low-resolution MS1 data, the units can be changed from \"ppm\" to \"Da\", e.g. 5Da.<br/>\
4. Select an output folder. Note that you need to have writing permissions to this folder.<br/>\
5. (optional) Select a (new) file for where to write the consensus spectra. Currently, the following formats are supported: mzML, ms2, mgf.<br/>\
6. (optional) Add extra command line parameters in the \"Extra command line parameters\" box. These advanced parameters can be found by running \"maracluster --help\" on the command line.<br/>\
7. Run MaRaCluster by pressing the \"Run MaRaCluster\" button.<br/>\
<br/>\
Log messages from MaRaCluster will be displayed in the panel on the right. If everything ran successfully you will see the message \"MaRaCluster execution finished successfully\". If you notice a problem with the input parameters or with the execution of MaRaCluster, you can stop the execution by pressing the \"Stop MaRaCluster\" button, which will replace the \"Run MaRaCluster\" button upon starting MaRaCluster.<br/>\
<br/>\
N.B.: MaRaCluster will try to reuse intermediate files present in the output folder produced by previous runs, so to start completely fresh, remember to empty the specified output folder.<br/>\
<br/>\
If you encounter any issues, please report them on <a href=\"https://github.com/statisticalbiotechnology/maracluster/issues\">GitHub</a>");
  
  // Add a layout for QDialog
  QHBoxLayout *dialog_layout = new QHBoxLayout(dialog);
  dialog->setLayout(dialog_layout);
  dialog->layout()->addWidget(scroll); // add scroll to the QDialog's layout
  dialog->show();
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
