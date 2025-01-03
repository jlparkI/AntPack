#ifndef MAIN_WINDOW_HEADER_H
#define MAIN_WINDOW_HEADER_H

// C++ headers
#include <QMainWindow>
#include <QMenuBar>
#include <QBoxLayout>
#include <QTabWidget>
#include <QTableWidget>



namespace GUIApplicationRunner {

class MainWindow : public QMainWindow {
 public:
    MainWindow(QWidget *parent = nullptr);

 private:
    void quit_now();
};


}  // namespace GUIApplicationRunner

#endif
