#ifndef DETSDWSYSTEMCONFIGFILEHANDLE_H_
#define DETSDWSYSTEMCONFIGFILEHANDLE_H_

#include <memory>
#include <fstream>

struct DetSDW_SystemConfig_FileHandle {
    std::unique_ptr<std::ofstream> phi_output_text;
    std::unique_ptr<std::ofstream> cdwl_output_text; 
    std::unique_ptr<std::ofstream> phi_output_binary;
    std::unique_ptr<std::ofstream> cdwl_output_binary;
};



#endif //DETSDWSYSTEMCONFIGFILEHANDLE_H_
