#ifndef DETSDWSYSTEMCONFIGFILEHANDLE_H_
#define DETSDWSYSTEMCONFIGFILEHANDLE_H_

#include <memory>
#include <fstream>

struct DetSDW_SystemConfig_FileHandle {
    std::unique_ptr<std::ofstream> phi_output_text;
    std::unique_ptr<std::ofstream> cdwl_output_text; 
    std::unique_ptr<std::ofstream> phi_output_binary;
    std::unique_ptr<std::ofstream> cdwl_output_binary;

    void flush() {
        if (phi_output_text) phi_output_text->flush();
        if (cdwl_output_text) cdwl_output_text->flush();        
        if (phi_output_binary) phi_output_binary->flush();
        if (cdwl_output_binary) cdwl_output_binary->flush();        
    }
};



#endif //DETSDWSYSTEMCONFIGFILEHANDLE_H_
