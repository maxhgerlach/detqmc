#ifndef DETSDWSYSTEMCONFIGFILEHANDLE_H_
#define DETSDWSYSTEMCONFIGFILEHANDLE_H_

#include <memory>
#include <fstream>

struct DetSDW_SystemConfig_FileHandle {
    typedef std::shared_ptr<std::ofstream> OfstreamPointer;
    OfstreamPointer phi_output_text;
    OfstreamPointer cdwl_output_text; 
    OfstreamPointer phi_output_binary;
    OfstreamPointer cdwl_output_binary;

    void flush() {
        if (phi_output_text) phi_output_text->flush();
        if (cdwl_output_text) cdwl_output_text->flush();        
        if (phi_output_binary) phi_output_binary->flush();
        if (cdwl_output_binary) cdwl_output_binary->flush();        
    }
};



#endif //DETSDWSYSTEMCONFIGFILEHANDLE_H_
