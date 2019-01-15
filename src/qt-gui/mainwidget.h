#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>
#include <QMainWindow>
#include <QProcess>

class QGridLayout;
class QPushButton;
class QListWidget;
class QTextBrowser;
class QFileDialog;
class QLineEdit;
class QDoubleSpinBox;
class QComboBox;
class QMessageBox;

// This is the declaration of our MainWidget class
// The definition/implementation is in mainwidget.cpp
class MainWidget : public QMainWindow
{
  Q_OBJECT

 public:
  explicit MainWidget(QWidget *parent = 0); //Constructor
  ~MainWidget(); // Destructor

 private slots:
  void onSelectFilesButtonReleased(); // Handler for "add files" button press
  void onDeleteFilesButtonReleased(); // Handler for "delete files" button press
  void onOutputFolderButtonReleased(); // Handler for "output folder" button press
  void onConsensusFileButtonReleased(); // Handler for "output folder" button press
  void onRunButtonReleased(); // Handler for "run" button press
  void onSaveStdoutButtonReleased();
  void onClearStdoutButtonReleased(); 
  void onCaptureProcessOutput(); // Handler for Process output
  
  void saveProject();
  void loadProject();
  void help();

 private:
  QString prepareCmd();
  
  QGridLayout* mainLayout_;
  
  QFileDialog* fileDialog_;
  QListWidget* fileListWidget_;
  
  QDoubleSpinBox* pvalThreshInput_;
  QDoubleSpinBox* precTolInput_;
  QComboBox* precTolUnitInput_;
  
  QLineEdit* outputFolderInput_;
  QFileDialog* outputFolderDialog_;
  
  QLineEdit* consensusFileInput_;
  
  QLineEdit* otherParamsInput_;
  
  QPushButton* runButton_;
  
  QTextBrowser* textBrowser_;
  QPushButton* saveStdoutButton_;
  QPushButton* clearStdoutButton_;
  
  QProcess process_;   // This is the process that runs MaRaCluster
};

#endif // MAINWIDGET_H
