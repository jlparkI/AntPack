#ifndef GUI_APPLICATION_RUNNER_HEADER_H
#define GUI_APPLICATION_RUNNER_HEADER_H
// C++ headers
#include <string>
#include <map>


namespace GUIApplicationRunner {

void run_antpack_gui(std::string consensus_filepath,
        std::map<std::string, std::string> icon_filepaths);

}  // namespace GUIApplicationRunner

#endif