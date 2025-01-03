// C++ headers
#include <QMainWindow>
#include <QApplication>

// Library headers

// Project headers
#include "main_window.h"


namespace GUIApplicationRunner {

MainWindow::MainWindow(QWidget *parent)
: QMainWindow(parent) {
}


/// @brief Quits the application. 
void MainWindow::quit_now(void) {
    QApplication::quit();
}


}  // namespace GUIApplicationRunner
