// C++ headers
#include <string>
#include <map>

// Library headers
#include <QApplication>
#include <QStyleFactory>
#include <QSplashScreen>
#include <QIcon>

// Project headers
#include "application_runner.h"
#include "main_window.h"

namespace GUIApplicationRunner {

void run_antpack_gui(std::string consensus_filepath,
        std::map<std::string, std::string> icon_filepaths) {
    int narg = 1;
    char *argv[1];
    char arg1[] = {'a', 'p', 'p'};
    argv[1] = arg1;
    QApplication app(narg, argv);

    // QPixmap pixmap("resources/intro_pic.jpg");
    // QSplashScreen splash(pixmap);
    // splash.show();

    // Sleep(3000);
    // splash.close();

    MainWindow w;
    w.setWindowTitle("AntPack");

    QString icon_filepath = QString::fromStdString(
            icon_filepaths.at("app_icon"));
    w.setWindowIcon(QIcon(icon_filepath));
    w.resize(800, 500);
    w.show();
    app.exec();
}

}  // namespace GUIApplicationRunner