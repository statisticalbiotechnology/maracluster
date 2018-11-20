#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>
#include <QProcess>

class QGridLayout;
class QPushButton;
class QListWidget;
class QTextBrowser;
class QFileDialog;
class QLineEdit;
class QDoubleSpinBox;
class QComboBox;

// This is the declaration of our MainWidget class
// The definition/implementation is in mainwidget.cpp
class MainWidget : public QWidget
{
  Q_OBJECT

 public:
  explicit MainWidget(QWidget *parent = 0); //Constructor
  ~MainWidget(); // Destructor

 private slots:
  void onSelectFilesButtonReleased(); // Handler for "add files" button press
  void onDeleteFilesButtonReleased(); // Handler for "delete files" button press
  void onOutputFolderButtonReleased(); // Handler for "output folder" button press
  void onRunButtonReleased(); // Handler for "run" button press
  void onSaveStdoutButtonReleased();
  void onClearStdoutButtonReleased(); 
  void onCaptureProcessOutput(); // Handler for Process output

 private:
  QGridLayout* mainLayout_;
  QPushButton* selectFileButton_;
  QPushButton* deleteFileButton_;
  QPushButton* runButton_;
  QPushButton* outputFolderButton_;
  QPushButton* saveStdoutButton_;
  QPushButton* clearStdoutButton_;
  QListWidget* fileListWidget_;
  QDoubleSpinBox* pvalThreshInput_;
  QDoubleSpinBox* precTolInput_;
  QComboBox* precTolUnitInput_;
  QLineEdit* otherParamsInput_;
  QLineEdit* outputFolderInput_;
  QTextBrowser* textBrowser_;
  QFileDialog* fileDialog_;
  QFileDialog* outputFolderDialog_;
  QProcess process_;   // This is the process the button will fire off
};

#endif // MAINWIDGET_H
